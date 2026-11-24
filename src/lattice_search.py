"""Exhaustive exact search for crescent configurations on a hexagonal patch of the
triangular lattice.

Lattice point (a, b) <-> complex a + b*w, w = exp(i*pi/3); Cartesian coordinates
(a + b/2, b*sqrt(3)/2).  Everything below is integer arithmetic:

* squared distance of (a,b),(a',b'):  p^2 + p q + q^2,  (p,q) = (a-a', b-b').
* collinearity: the Cartesian 2x2 determinant equals (sqrt(3)/2) * (p1 q2 - p2 q1),
  so A,B,C collinear  <=>  integer det of the (a,b)-differences is 0.
* concyclicity: the Cartesian lifted determinant |X Y X^2+Y^2 1| has columns
  X = a + b/2, Y = (sqrt(3)/2) b, X^2+Y^2 = N := a^2+ab+b^2.  Subtracting 1/2 of the
  b-column (a multiple of the Y column) from X gives
      |X Y N 1| = (sqrt(3)/2) * |a b N 1|,
  so 4 points are concyclic (given no 3 collinear)  <=>  integer det |a b N 1| = 0.
  (Both criteria are cross-checked against src/checker.py in the tests.)

Patch of radius R: { (a,b) : max(|a|, |b|, |a+b|) <= R }  (3R^2+3R+1 points).

Usage:  python src/lattice_search.py R k [k ...] [--all] [--out data/known_examples.json]
Without --all it stops at the first solution for each k.
"""
from __future__ import annotations

import itertools
import json
import sys
import time


def patch(R):
    return [(a, b) for a in range(-R, R + 1) for b in range(-R, R + 1) if abs(a + b) <= R]


def norm(p, q):
    return p * p + p * q + q * q


def collinear(A, B, C):
    return (B[0] - A[0]) * (C[1] - A[1]) - (B[1] - A[1]) * (C[0] - A[0]) == 0


def det3(m):
    (a, b, c), (d, e, f), (g, h, i) = m
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def concyclic(A, B, C, D):
    rows = []
    for P in (A, B, C):
        p, q = P[0] - D[0], P[1] - D[1]
        rows.append((p, q, norm(p, q)))
    return det3(rows) == 0


def to_cartesian_strings(pts):
    """Exact sqrt(3) strings (x, y) for the checker's 'sqrt' input mode."""
    from fractions import Fraction
    out = []
    for a, b in pts:
        x = Fraction(2 * a + b, 2)
        y = Fraction(b, 2)
        out.append([str(x), f"({y})*sqrt(3)"])
    return out


# The 12 symmetries of the lattice fixing the origin, acting on (a,b) coordinates.
def _rot(P):  # multiply by w: (a + b w) w = a w + b w^2 = a w + b (w - 1) = -b + (a + b) w
    a, b = P
    return (-b, a + b)


def _refl(P):  # complex conjugation: conj(w) = 1 - w  -> a + b - b w
    a, b = P
    return (a + b, -b)


def canonical(pts):
    best = None
    cur = list(pts)
    for r in range(6):
        for s in (False, True):
            q = [_refl(P) for P in cur] if s else cur
            m0 = min(q)
            t = tuple(sorted((P[0] - m0[0], P[1] - m0[1]) for P in q))
            if best is None or t < best:
                best = t
        cur = [_rot(P) for P in cur]
    return best


def search(R, k, find_all=False, time_limit=None):
    pts = patch(R)
    N = len(pts)
    D = [[norm(pts[i][0] - pts[j][0], pts[i][1] - pts[j][1]) for j in range(N)] for i in range(N)]
    sols = []
    canon = set()
    t0 = time.time()
    nodes = 0
    chosen = []
    counts = {}

    def feasible():
        if len(counts) > k - 1:
            return False
        cs = sorted(counts.values(), reverse=True)
        return all(c <= k - 1 - i for i, c in enumerate(cs))

    class Stop(Exception):
        pass

    def rec(start):
        nonlocal nodes
        nodes += 1
        if len(chosen) == k:
            if sorted(counts.values()) == list(range(1, k)):
                conf = [pts[i] for i in chosen]
                c = canonical(conf)
                if c not in canon:
                    canon.add(c)
                    sols.append(conf)
                    if not find_all:
                        raise Stop
            return
        if time_limit and time.time() - t0 > time_limit:
            raise Stop
        need = k - len(chosen)
        for p in range(start, N - need + 1):
            P = pts[p]
            # general position against the current set
            bad = False
            for i, j in itertools.combinations(chosen, 2):
                if collinear(pts[i], pts[j], P):
                    bad = True
                    break
            if bad:
                continue
            for i, j, l in itertools.combinations(chosen, 3):
                if concyclic(pts[i], pts[j], pts[l], P):
                    bad = True
                    break
            if bad:
                continue
            added = []
            for i in chosen:
                d = D[i][p]
                counts[d] = counts.get(d, 0) + 1
                added.append(d)
            if feasible():
                chosen.append(p)
                rec(p + 1)
                chosen.pop()
            for d in added:
                counts[d] -= 1
                if counts[d] == 0:
                    del counts[d]

    try:
        rec(0)
        complete = True
    except Stop:
        complete = False  # stopped early (first solution found, or time limit)
    return sols, nodes, time.time() - t0, complete


def main(argv):
    find_all = "--all" in argv
    out = None
    if "--out" in argv:
        out = argv[argv.index("--out") + 1]
    nums = [int(a) for a in argv if a.lstrip("-").isdigit()]
    R, ks = nums[0], nums[1:]
    record = {"lattice": "triangular, (a,b) -> (a + b/2, b*sqrt(3)/2)", "patch_radius": R,
              "patch_size": len(patch(R)), "results": {}}
    for k in ks:
        sols, nodes, dt, complete = search(R, k, find_all)
        print(f"R={R} k={k}: {len(sols)} solution(s) up to lattice symmetry"
              f"{' (exhaustive)' if find_all and complete else ''}; nodes={nodes} time={dt:.1f}s")
        for s in sols[:5]:
            print("   ", s)
        record["results"][str(k)] = {
            "exhaustive": bool(find_all and complete), "num_up_to_symmetry": len(sols),
            "nodes": nodes, "seconds": round(dt, 1),
            "examples": [{"lattice_ab": [list(P) for P in s],
                          "checker_input": {"field": {"type": "sqrt", "d": 3},
                                            "points": to_cartesian_strings(s)}}
                         for s in sols[:10]]}
    if out:
        with open(out, "w", encoding="utf-8") as fh:
            json.dump(record, fh, indent=1)
        print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1:])
