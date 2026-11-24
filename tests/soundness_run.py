"""Builder E soundness test (Phase 2, mandatory).

Generate REALIZABLE 6-point patterns from exact point sets in general position with <= 6
distinct distances, and check every one's canonical pattern is in
data/patterns6_local.txt.gz AND data/patterns6_affine.txt.gz.

Point sources (all exact integer arithmetic):
  * Lattices Z^2 with Euclidean metric given by a positive definite integral binary quadratic
    form Q(dx,dy) = a dx^2 + b dx dy + c dy^2 (b^2 < 4ac): this is the squared Euclidean
    distance of the lattice spanned by vectors e1,e2 with |e1|^2=a, 2 e1.e2=b, |e2|^2=c, so every
    finite subset is a genuine planar point set. (1,0,1)=square lattice, (1,1,1)=triangular
    lattice (coords (2a+b, b*sqrt3)/2, squared distance a^2+ab+b^2), (3,0,1)=the lattice
    {(u*sqrt3/2, v/2)} containing Palasti's 8 points (scaled by 2), plus others.
    Collinearity of 3 lattice points <=> integer cross product of coordinate differences = 0
    (linear change of basis preserves collinearity). For 4 points with no 3 collinear:
    concyclic <=> det(D_ij)_{4x4} = 0 (Wang eq. (circlefactorization): det = -4 (det B)^2).
  * Palasti's 8-point set: all 6-subsets with <= 6 distances that are in general position.
Translations: WLOG the lex-smallest point is (0,0).

Usage: python tests/soundness_run.py [W] [procs]
"""
import gzip
import itertools
import json
import os
import sys
import time
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
import canon  # noqa: E402

FORMS = [(1, 0, 1), (1, 1, 1), (3, 0, 1), (1, 0, 2), (1, 1, 2), (1, 0, 3), (2, 1, 2),
         (1, 1, 3), (2, 1, 3), (1, 0, 5), (2, 2, 3), (1, 1, 4), (3, 2, 3), (1, 0, 7),
         (2, 1, 4), (3, 3, 4)]

PALASTI_UV = [(0, 2), (2, 0), (4, 0), (5, 5), (3, 9), (1, 7), (3, 7), (2, 4)]  # (2x, 2y), point=(x*sqrt3, y)
PALASTI_FORM = (3, 0, 1)


def qdist(f, p, q):
    a, b, c = f
    dx, dy = p[0] - q[0], p[1] - q[1]
    return a * dx * dx + b * dx * dy + c * dy * dy


def collinear(p, q, r):
    return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]) == 0


def det4(M):
    # exact integer determinant via Laplace (4x4)
    def det3(m):
        return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))
    s = 0
    for j in range(4):
        minor = [row[:j] + row[j + 1:] for row in M[1:]]
        s += (-1) ** j * M[0][j] * det3(minor)
    return s


def concyclic(f, pts):
    M = [[qdist(f, p, q) for q in pts] for p in pts]
    return det4(M) == 0


def general_position(f, pts):
    for t in itertools.combinations(pts, 3):
        if collinear(*t):
            return False
    for q in itertools.combinations(pts, 4):
        if concyclic(f, q):
            return False
    return True


def pattern_of(f, pts):
    labs = [qdist(f, pts[i], pts[j]) for i in range(len(pts)) for j in range(i + 1, len(pts))]
    return labs


def search(args):
    """Backtracking over 6-subsets containing (0,0) as lex-min point, with second point fixed."""
    f, W, second_idx = args
    cands = [(x, y) for x in range(0, W + 1) for y in range(-W, W + 1) if (x, y) > (0, 0)]
    found = {}  # relabelled lex string -> example pts
    nsets = 0
    origin = (0, 0)

    def rec(chosen, start, dists):
        nonlocal nsets
        if len(chosen) == 6:
            nsets += 1
            key = canon.relabel(pattern_of(f, chosen))
            if key not in found:
                found[key] = list(chosen)
            return
        for idx in range(start, len(cands)):
            p = cands[idx]
            nd = set(dists)
            ok = True
            for q in chosen:
                nd.add(qdist(f, p, q))
                if len(nd) > 6:
                    ok = False
                    break
            if not ok:
                continue
            bad = False
            for q, r in itertools.combinations(chosen, 2):
                if collinear(q, r, p):
                    bad = True
                    break
            if bad:
                continue
            for trip in itertools.combinations(chosen, 3):
                if concyclic(f, list(trip) + [p]):
                    bad = True
                    break
            if bad:
                continue
            rec(chosen + [p], idx + 1, nd)

    p2 = cands[second_idx]
    rec([origin, p2], second_idx + 1, {qdist(f, origin, p2)})
    return f, nsets, found


def main():
    W = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    procs = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    forms = FORMS
    if len(sys.argv) > 3 and sys.argv[3] != "all":
        forms = [tuple(int(x) for x in t.split(",")) for t in sys.argv[3].split(";")]
    tag = f"W{W}" + ("" if forms is FORMS else "_" + "_".join("".join(map(str, f)) for f in forms))
    t0 = time.time()
    with gzip.open(os.path.join(ROOT, "data", "patterns6_local.txt.gz"), "rt") as fh:
        local = set(ln.strip() for ln in fh)
    with gzip.open(os.path.join(ROOT, "data", "patterns6_affine.txt.gz"), "rt") as fh:
        affine = set(ln.strip() for ln in fh)
    print(f"loaded local={len(local)} affine={len(affine)}")

    examples = {}  # canonical -> (source, pts)
    rawcount = Counter()

    # Palasti
    P = PALASTI_UV
    dd = Counter(qdist(PALASTI_FORM, P[i], P[j]) for i in range(8) for j in range(i + 1, 8))
    pal8_ok = general_position(PALASTI_FORM, P) and sorted(dd.values()) == list(range(1, 8))
    print("Palasti 8-set (u,v)=2*(x/sqrt3... ) multiplicities:", sorted(dd.values()), "general position+crescent:", pal8_ok)
    for sub in itertools.combinations(P, 6):
        sub = list(sub)
        labs = pattern_of(PALASTI_FORM, sub)
        if len(set(labs)) <= 6 and general_position(PALASTI_FORM, sub):
            rawcount["palasti"] += 1
            examples.setdefault(canon.canonical(labs, 6), ("palasti", sub))

    import multiprocessing as mp
    W_cands = (W + 1) * (2 * W + 1) - (W + 1)  # number of points lex-greater than origin
    jobs = [(f, W, i) for f in forms for i in range(W_cands)]
    with mp.Pool(procs) as pool:
        for f, nsets, found in pool.imap_unordered(search, jobs):
            rawcount[str(f)] += nsets
            for key, pts in found.items():
                c = canon.canonical(key, 6)
                if c not in examples:
                    examples[c] = (str(f), pts)
    print("raw general-position <=6-distance 6-sets by source:", dict(rawcount))
    print("distinct canonical realizable patterns:", len(examples))

    miss_local = [(c, e) for c, e in examples.items() if c not in local]
    miss_aff = [(c, e) for c, e in examples.items() if c not in affine]
    byr = Counter(max(int(ch) for ch in c) for c in examples)
    print("realizable patterns by r:", dict(sorted(byr.items())))
    print("MISSING from local:", len(miss_local))
    print("MISSING from affine:", len(miss_aff))
    for c, e in (miss_local + miss_aff)[:20]:
        print("  !!! MISS", c, e)
    out = {
        "W": W, "forms": forms, "raw_sets": dict(rawcount), "distinct_patterns": len(examples),
        "by_r": dict(sorted(byr.items())), "missing_local": [(c, e) for c, e in miss_local],
        "missing_affine": [(c, e) for c, e in miss_aff], "palasti8_ok": pal8_ok,
        "time_s": round(time.time() - t0, 1),
    }
    with open(os.path.join(ROOT, "data", f"soundness_phase2_{tag}.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    with gzip.open(os.path.join(ROOT, "data", f"soundness_realizable6_{tag}.txt.gz"), "wt", newline="\n") as fh:
        for c in sorted(examples):
            src, pts = examples[c]
            fh.write(f"{c}\t{src}\t{pts}\n")
    print(f"time {time.time() - t0:.1f}s")
    return 0 if not miss_local and not miss_aff else 1


if __name__ == "__main__":
    sys.exit(main())
