"""Verifier: exhaustive search of the radius-3 hexagonal patch (37 points) of the
triangular lattice for n-point crescent configurations.

Lattice point (i, j) = i*(1,0) + j*(1/2, sqrt3/2); squared distance = di^2 + di*dj + dj^2.
Backtracking over index-increasing subsets with exact integer pruning:
  - no 3 collinear (cross product in lattice coords; collinearity is affine invariant),
  - #distinct distances <= n-1 and k-th largest class size <= n-k (necessary for final
    sizes {1..n-1}, since sizes only grow).
Every leaf is re-checked with the exact Q(sqrt3) checker (verify/checker_v.py), which
also does the concyclicity test. Results are reduced modulo the 12 lattice point
symmetries + translation.
Usage: python lattice_search.py n [radius]
"""
import sys
import os
import time
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from checker_v import check, lattice_to_q3, Q3Field


def patch(r):
    return [(i, j) for i in range(-r, r + 1) for j in range(-r, r + 1) if abs(i + j) <= r]


def norm(a, b):
    di, dj = a[0] - b[0], a[1] - b[1]
    return di * di + di * dj + dj * dj


def canon(pts):
    best = None
    cur = list(pts)
    for refl in (False, True):
        q = [(j, i) for i, j in cur] if refl else list(cur)
        for _ in range(6):
            q = [(-j, i + j) for i, j in q]
            s = sorted(q)
            o = s[0]
            t = tuple((a - o[0], b - o[1]) for a, b in s)
            if best is None or t < best:
                best = t
    return best


def search(n, r):
    P = patch(r)
    N = len(P)
    D = [[norm(P[a], P[b]) for b in range(N)] for a in range(N)]
    col = {}

    def collin(a, b, c):
        (x1, y1), (x2, y2), (x3, y3) = P[a], P[b], P[c]
        return (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1) == 0

    m = n - 1
    sol = []
    counts = {}
    chosen = []
    leaves = [0]

    def feasible():
        if len(counts) > m:
            return False
        vs = sorted(counts.values(), reverse=True)
        for k, v in enumerate(vs, 1):
            if v > n - k:
                return False
        return True

    def rec(start):
        if len(chosen) == n:
            leaves[0] += 1
            vs = sorted(counts.values())
            if vs == list(range(1, n)):
                pts = [lattice_to_q3(*P[x]) for x in chosen]
                ok, info = check(pts, Q3Field)
                if ok:
                    sol.append([P[x] for x in chosen])
            return
        for c in range(start, N - (n - len(chosen)) + 1):
            if any(collin(a, b, c) for a, b in combinations(chosen, 2)):
                continue
            added = []
            for a in chosen:
                d = D[a][c]
                counts[d] = counts.get(d, 0) + 1
                added.append(d)
            if feasible():
                chosen.append(c)
                rec(c + 1)
                chosen.pop()
            for d in added:
                counts[d] -= 1
                if counts[d] == 0:
                    del counts[d]

    rec(0)
    return P, sol, leaves[0]


if __name__ == "__main__":
    n = int(sys.argv[1])
    r = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    t = time.time()
    P, sol, leaves = search(n, r)
    classes = {}
    for s in sol:
        classes.setdefault(canon(s), s)
    print(f"n={n} radius={r} patch={len(P)} pts: {len(sol)} raw solutions, "
          f"{len(classes)} classes mod lattice symmetry, {leaves} size-n leaves, "
          f"{time.time() - t:.1f}s")
    for k, (c, s) in enumerate(sorted(classes.items())):
        ok, info = check([lattice_to_q3(*p) for p in s], Q3Field)
        print(f"  class {k}: lattice(i,j) = {list(c)}  sizes={info['class_sizes']}")
