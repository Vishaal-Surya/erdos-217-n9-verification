"""Full Engine-1 run over data/patterns6_affine.txt.gz with checkpointing (Builder A).

python src/algebra_full.py [--workers 4] [--chunk 5000]
Chunks -> data/algebra_E1_chunks/chunk_NNNN.jsonl.gz (skipped if present => resumable).
Engine 1 is run WITH its own affine tests (Lemma 2.3) re-applied: they are cheap and also
produce the linear relations Az=0 that enter W.  An affine exclusion here would be a
discrepancy with Builder E's file and is recorded as such.
Every quartic exclusion's certificate is independently replayed over Q (verify_certificate).
"""
import gzip
import json
import os
import sys
import time
from multiprocessing import get_context

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CH = os.path.join(ROOT, "data", "algebra_E1_chunks")


def work(s):
    import algebra as A
    try:
        r = A.engine1(s)
    except Exception as ex:
        return {"pattern": s, "verdict": "ERROR", "error": repr(ex)}
    out = {k: r.get(k) for k in ("pattern", "r", "verdict", "stage", "reason", "a", "h", "finite",
                                  "dimW", "modp", "modp_discrepancy", "t")}
    if r.get("stage") == "quartic":
        out["cert"] = r["cert"]
        t0 = time.perf_counter()
        out["cert_replay_ok"] = A.verify_certificate(s, r["cert"])
        out["t_replay"] = time.perf_counter() - t0
    return out


def main(argv):
    workers = int(argv[argv.index("--workers") + 1]) if "--workers" in argv else 4
    chunk = int(argv[argv.index("--chunk") + 1]) if "--chunk" in argv else 5000
    os.makedirs(CH, exist_ok=True)
    pats = gzip.open(os.path.join(ROOT, "data", "patterns6_affine.txt.gz"), "rt").read().split()
    nchunks = (len(pats) + chunk - 1) // chunk
    ctx = get_context("spawn")
    with ctx.Pool(workers) as pool:
        for k in range(nchunks):
            fn = os.path.join(CH, f"chunk_{k:04d}.jsonl.gz")
            if os.path.exists(fn):
                continue
            t0 = time.time()
            part = pats[k * chunk:(k + 1) * chunk]
            res = pool.map(work, part, chunksize=20)
            tmp = fn + ".tmp"
            with gzip.open(tmp, "wt") as f:
                for r in res:
                    f.write(json.dumps(r) + "\n")
            os.replace(tmp, fn)
            ret = sum(r["verdict"] in ("F", "U") for r in res)
            print(f"chunk {k+1}/{nchunks} done in {time.time()-t0:.0f}s, retained {ret}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
