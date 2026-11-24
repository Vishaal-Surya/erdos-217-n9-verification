"""Engine 1 vs Engine 2 comparison + timing on a sample (Builder A).

python src/algebra_compare.py N_RANDOM N_FILE OUT_PREFIX [--workers 2] [--timeout 60]
 - N_RANDOM random label strings passing the local tests (Lemma 2.1/2.2, implemented in algebra.local_ok)
 - N_FILE patterns sampled from data/patterns6_local.txt.gz (if present)
Writes OUT_PREFIX.jsonl (per-pattern results) and OUT_PREFIX_summary.json.
"""
import gzip
import json
import os
import random
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import algebra as A  # noqa: E402
from algebra_batch import run_batch  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def sample_random(n, seed=101):
    rng = random.Random(seed)
    out = set()
    while len(out) < n:
        s = A.normalize([rng.randint(1, 6) for _ in range(15)])
        if A.local_ok(A.parse(s)[0]):
            out.add(s)
    return sorted(out)


def sample_file(n, seed=202):
    p = os.path.join(ROOT, "data", "patterns6_local.txt.gz")
    if not os.path.exists(p) or n == 0:
        return []
    L = gzip.open(p, "rt").read().split()
    return sorted(random.Random(seed).sample(L, n))


def pct(ts, q):
    ts = sorted(ts)
    return ts[min(len(ts) - 1, int(q * len(ts)))] if ts else None


def main(argv):
    nr, nf, prefix = int(argv[0]), int(argv[1]), argv[2]
    workers = int(argv[argv.index("--workers") + 1]) if "--workers" in argv else 2
    timeout = float(argv[argv.index("--timeout") + 1]) if "--timeout" in argv else 60.0
    groups = {"random_local": sample_random(nr), "file_local": sample_file(nf)}
    summary = {}
    with open(prefix + ".jsonl", "w") as f:
        for g, pats in groups.items():
            if not pats:
                continue
            res = run_batch(pats, (1, 2), workers, timeout, progress=250)
            mat = Counter()
            t1, t2 = [], []
            bad = []
            for r in res:
                f.write(json.dumps({"group": g, **r}) + "\n")
                e1 = (r.get("e1") or {})
                e2 = (r.get("e2") or {})
                v1 = e1.get("verdict", "ERR")
                if v1 == "EXCLUDED":
                    v1 = "EXCL-" + (e1.get("stage") or "?")
                v2 = e2.get("verdict", "ERR")
                mat[(v1, v2)] += 1
                if "t" in e1:
                    t1.append(e1["t"])
                if "t" in e2:
                    t2.append(e2["t"])
                if v1.startswith("EXCL") and v2 in ("ZERODIM", "POSDIM"):
                    bad.append(r["pattern"])
                if v1 in ("F", "U") and v2 == "EXCLUDED":
                    pass  # allowed: E2 is stronger (full saturation)
            summary[g] = {
                "n": len(pats),
                "matrix": {f"{a} | {b}": c for (a, b), c in sorted(mat.items())},
                "SERIOUS_e1excl_e2feasible": bad,
                "e1_time": {"median": pct(t1, .5), "p99": pct(t1, .99), "max": max(t1) if t1 else None,
                            "mean": sum(t1) / len(t1) if t1 else None},
                "e2_time": {"median": pct(t2, .5), "p99": pct(t2, .99), "max": max(t2) if t2 else None,
                            "mean": sum(t2) / len(t2) if t2 else None},
            }
            print(g, json.dumps(summary[g], indent=1), flush=True)
    with open(prefix + "_summary.json", "w") as f:
        json.dump(summary, f, indent=1)


if __name__ == "__main__":
    main(sys.argv[1:])
