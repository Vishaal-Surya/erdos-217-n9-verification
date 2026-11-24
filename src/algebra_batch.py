"""Batch driver for src/algebra.py with per-pattern timeouts (killable worker processes).

Usage: python src/algebra_batch.py IN.txt[.gz] OUT.jsonl [--engines 1,2] [--workers 2] [--timeout 60]
Each output line: {"pattern", "e1": {...}, "e2": {...}}.  A timed-out pattern is recorded with
verdict "UNRESOLVED" (never as excluded).  CPU cap: workers default 2 (user thermal limit).
"""
from __future__ import annotations

import gzip
import json
import multiprocessing as mp
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

MAX_WORKERS = min(int(os.environ.get("ALGEBRA_MAX_WORKERS", "2")), 5)  # hard cap 5 (user CPU thermal limit)


def _clean(d):
    return {k: (v if isinstance(v, (int, float, str, bool, type(None), list, dict)) else str(v))
            for k, v in d.items() if not k.startswith("_")}


def run_one(s, engines, e2_limits=tuple(int(x) for x in os.environ.get("ALGEBRA_E2_LIMITS", "10000,100000,100000").split(","))):
    import algebra as A
    out = {"pattern": s}
    if 1 in engines:
        out["e1"] = _clean(A.engine1(s))
    if 2 in engines:
        out["e2"] = _clean(A.engine2(s, limits=e2_limits, eliminants=True))
    return out


def _worker(inq, outq, engines):
    while True:
        item = inq.get()
        if item is None:
            return
        idx, s = item
        try:
            outq.put((idx, run_one(s, engines)))
        except Exception as ex:  # recorded, never treated as exclusion
            outq.put((idx, {"pattern": s, "error": repr(ex)}))


def run_batch(patterns, engines=(1, 2), workers=2, timeout=60.0, progress=None):
    workers = min(workers, MAX_WORKERS)
    ctx = mp.get_context("spawn")
    outq = ctx.Queue()
    slots = []

    def spawn():
        q = ctx.Queue()
        p = ctx.Process(target=_worker, args=(q, outq, tuple(engines)), daemon=True)
        p.start()
        return {"p": p, "q": q, "task": None, "t0": None}

    for _ in range(workers):
        slots.append(spawn())
    results = [None] * len(patterns)
    nxt = 0
    done = 0
    while done < len(patterns):
        for i, sl in enumerate(slots):
            if sl["task"] is None and nxt < len(patterns):
                sl["task"] = nxt
                sl["t0"] = time.time()
                sl["q"].put((nxt, patterns[nxt]))
                nxt += 1
        try:
            idx, res = outq.get(timeout=0.2)
            results[idx] = res
            done += 1
            for sl in slots:
                if sl["task"] == idx:
                    sl["task"] = None
            if progress and done % progress == 0:
                print(f"{done}/{len(patterns)}", flush=True)
        except Exception:
            pass
        for i, sl in enumerate(slots):
            if sl["task"] is not None and time.time() - sl["t0"] > timeout:
                idx = sl["task"]
                sl["p"].kill()
                results[idx] = {"pattern": patterns[idx], "timeout": timeout,
                                "e1": {"verdict": "UNRESOLVED", "reason": "timeout"} if 1 in engines else None,
                                "e2": {"verdict": "UNRESOLVED", "reason": "timeout"} if 2 in engines else None}
                done += 1
                slots[i] = spawn()
    for sl in slots:
        sl["q"].put(None)
    return results


def main(argv):
    inp, outp = argv[0], argv[1]
    engines = (1, 2)
    workers, timeout = 2, 60.0
    if "--engines" in argv:
        engines = tuple(int(x) for x in argv[argv.index("--engines") + 1].split(","))
    if "--workers" in argv:
        workers = int(argv[argv.index("--workers") + 1])
    if "--timeout" in argv:
        timeout = float(argv[argv.index("--timeout") + 1])
    op = gzip.open if inp.endswith(".gz") else open
    with op(inp, "rt") as f:
        pats = [x.strip() for x in f if x.strip()]
    res = run_batch(pats, engines, workers, timeout, progress=500)
    with open(outp, "w") as f:
        for r in res:
            f.write(json.dumps(r) + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
