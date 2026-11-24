"""Builder C: exact affine post-filter (Wang Lemma 2.3 tests) for complete n-point patterns.

Input: lines "<arrangement> <class-size string>" (as written by CoverageSearch --out).
R rows: e_a+e_d-e_b-e_c for every ordered quadruple (a,b,c,d) with c_ab=c_cd and c_ac=c_bd
(parallelogram a+d=b+c, Lemma 2.2 proof); optional centroid rows (--centroid).
T1: reject if rowspace(R) has a nonzero vector supported on <=3 points (two points equal or three
    collinear: rows have coordinate sum 0, so such a vector is an affine dependency).
T2: Sum_i rho_i (D_ij - D_i1) = 0 for all rows rho and all j  ->  A z = 0 ; reject if z_l = 0 or
    z_l = z_m is forced (e_l or e_l-e_m in rowspace(A)).
Applying the tests to the complete pattern subsumes applying them to every prefix (prefix rows are
rows of the full matrices). Ranks: python-flint, cross-checked with Fractions (--check).
"""
import itertools
import sys
from fractions import Fraction

import flint


def rank(rows, ncols):
    if not rows:
        return 0
    return flint.fmpz_mat(len(rows), ncols, [x for r in rows for x in r]).rank()


def rank_frac(rows, ncols):
    M = [[Fraction(x) for x in r] for r in rows]
    rk = 0
    for col in range(ncols):
        piv = next((i for i in range(rk, len(M)) if M[i][col] != 0), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][col] != 0:
                f = M[i][col] / M[rk][col]
                M[i] = [x - f * y for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk


def test(s, n, centroid=False, rk=rank):
    c = [[-1] * n for _ in range(n)]
    for (i, j), ch in zip(itertools.combinations(range(n), 2), s):
        c[i][j] = c[j][i] = int(ch) - 1
    r = n - 1
    R = set()
    for a, b, cc, d in itertools.permutations(range(n), 4):
        if c[a][b] == c[cc][d] and c[a][cc] == c[b][d]:
            v = [0] * n; v[a] += 1; v[d] += 1; v[b] -= 1; v[cc] -= 1
            R.add(tuple(v))
    if centroid:
        for a, b, cc in itertools.combinations(range(n), 3):
            if c[a][b] == c[a][cc] == c[b][cc]:
                for d in range(n):
                    if d not in (a, b, cc) and c[a][d] == c[b][d] == c[cc][d]:
                        v = [0] * n; v[a] += 1; v[b] += 1; v[cc] += 1; v[d] -= 3
                        R.add(tuple(v))
    R = sorted(R)
    if not R:
        return True, "noR"
    rkR = rk(R, n)
    for I in itertools.combinations(range(n), 3):
        keep = [k for k in range(n) if k not in I]
        if rk([[row[k] for k in keep] for row in R], n - 3) < rkR:
            return False, "T1:" + ",".join(str(x + 1) for x in I)
    A = set()
    for rho in R:
        for j in range(1, n):
            coef = [0] * r
            for i in range(n):
                if rho[i] == 0:
                    continue
                if i != j:
                    coef[c[i][j]] += rho[i]
                if i != 0:
                    coef[c[i][0]] -= rho[i]
            if any(coef):
                A.add(tuple(coef))
    A = sorted(A)
    if not A:
        return True, "noA"
    rkA = rk(A, r)
    for l in range(r):
        v = [0] * r; v[l] = 1
        if rk(A + [tuple(v)], r) == rkA:
            return False, f"T2:z{l+1}=0"
    for l, m in itertools.combinations(range(r), 2):
        v = [0] * r; v[l] = 1; v[m] = -1
        if rk(A + [tuple(v)], r) == rkA:
            return False, f"T2:z{l+1}=z{m+1}"
    return True, "pass"


def main():
    inp, outp = sys.argv[1], sys.argv[2]
    centroid = "--centroid" in sys.argv
    check = "--check" in sys.argv
    rows = [l.split() for l in open(inp) if l.strip()]
    n = 9
    stats = {}
    with open(outp, "w") as f:
        for arr, s in rows:
            ok, why = test(s, n, centroid)
            if check:
                assert (ok, why) == test(s, n, centroid, rank_frac)
            stats.setdefault(arr, [0, 0])[0 if ok else 1] += 1
            f.write(f"{arr} {s} {'PASS' if ok else 'FAIL'} {why}\n")
    tot = [0, 0]
    for a, (p, q) in stats.items():
        print(f"{a:18s} pass={p} fail={q}")
        tot[0] += p; tot[1] += q
    print(f"TOTAL pass={tot[0]} fail={tot[1]} (centroid rows: {centroid})")


if __name__ == "__main__":
    main()
