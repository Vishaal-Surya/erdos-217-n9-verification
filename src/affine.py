"""Builder E, Phase 2: six-point affine tests (Wang, arXiv:2611.03688, Lemma 2.3), exact arithmetic.

Input pattern: string of 15 label digits over pairs (1,2),(1,3),...,(1,6),(2,3),...,(5,6).

Tests (all exact; integer matrices, rank over Q):
  R  = rows e_a+e_d-e_b-e_c for every ordered quadruple (a,b,c,d) of distinct points with
       c_ab=c_cd and c_ac=c_bd  (parallelogram relation a+d=b+c, Lemma 2.2)
     + rows e_a+e_b+e_c-3e_d whenever c_ab=c_ac=c_bc and c_ad=c_bd=c_cd (centroid).
  (T1) reject if rowspace(R) contains a nonzero vector with <=3 nonzero entries, tested as:
       exists triple I with rank R[:, I^c] < rank R.
  (T2) linear distance equations: for each row rho of R and j=2..6,
       sum_i rho_i (D_ij - D_i1) = 0, D_ii=0, D_ij=z_{c_ij}  ->  A z = 0.
       reject if z_l vanishes identically on ker A (zero column of a kernel basis N), or
       z_l - z_m vanishes identically on ker A (two equal columns of N).
       Basis independence: rows of any basis N span ker A, so column l of N is zero iff every
       kernel vector has z_l = 0, and columns l,m agree iff every kernel vector has z_l = z_m.
       We test it in the basis-free form  e_l in rowspace(A)  /  e_l - e_m in rowspace(A),
       i.e. rank([A; v]) == rank(A)  (ker A = rowspace(A)^perp over Q).
Two engines: python-flint fmpz_mat ranks (fast) and pure-Python Fraction elimination
(independent cross-check, used in tests).
"""
from __future__ import annotations

import gzip
import itertools
import sys
import time
from fractions import Fraction

PAIRS = [(i, j) for i in range(6) for j in range(i + 1, 6)]
TRIPLES = list(itertools.combinations(range(6), 3))


def to_matrix(s: str):
    c = [[-1] * 6 for _ in range(6)]
    for (i, j), ch in zip(PAIRS, s):
        c[i][j] = c[j][i] = int(ch) - 1
    return c


def build_R(c):
    rows = set()
    n = 6
    for a, b, cc, d in itertools.permutations(range(n), 4):
        if c[a][b] == c[cc][d] and c[a][cc] == c[b][d]:
            v = [0] * n
            v[a] += 1; v[d] += 1; v[b] -= 1; v[cc] -= 1
            rows.add(tuple(v))
    for tri in TRIPLES:
        a, b, cc = tri
        if c[a][b] == c[a][cc] == c[b][cc]:
            for d in range(n):
                if d in tri:
                    continue
                if c[a][d] == c[b][d] == c[cc][d]:
                    v = [0] * n
                    v[a] += 1; v[b] += 1; v[cc] += 1; v[d] -= 3
                    rows.add(tuple(v))
    return sorted(rows)


def build_A(c, R, r):
    eqs = []
    for rho in R:
        for j in range(1, 6):
            coef = [0] * r
            for i in range(6):
                if rho[i] == 0:
                    continue
                if i != j:
                    coef[c[i][j]] += rho[i]
                if i != 0:
                    coef[c[i][0]] -= rho[i]
            if any(coef):
                eqs.append(tuple(coef))
    return sorted(set(eqs))


# ---------------- rank engines ----------------
def rank_fraction(rows, ncols):
    M = [[Fraction(x) for x in row] for row in rows]
    rank = 0
    for col in range(ncols):
        piv = None
        for i in range(rank, len(M)):
            if M[i][col] != 0:
                piv = i
                break
        if piv is None:
            continue
        M[rank], M[piv] = M[piv], M[rank]
        p = M[rank][col]
        for i in range(len(M)):
            if i != rank and M[i][col] != 0:
                f = M[i][col] / p
                M[i] = [x - f * y for x, y in zip(M[i], M[rank])]
        rank += 1
    return rank


def _rank_flint_factory():
    import flint

    def rank_flint(rows, ncols):
        if not rows:
            return 0
        return flint.fmpz_mat(len(rows), ncols, [x for row in rows for x in row]).rank()

    return rank_flint


def affine_test(s: str, engine: str = "flint"):
    """Return (passes: bool, reason: str)."""
    rank = rank_fraction if engine == "fraction" else _rank_flint_factory()
    c = to_matrix(s)
    r = max(int(ch) for ch in s)
    R = build_R(c)
    if not R:
        return True, "noR"
    rkR = rank(R, 6)
    for I in TRIPLES:
        keep = [k for k in range(6) if k not in I]
        sub = [[row[k] for k in keep] for row in R]
        if rank(sub, 3) < rkR:
            return False, "support<=3:%s" % (",".join(str(x + 1) for x in I))
    A = build_A(c, R, r)
    if not A:
        return True, "noA"
    rkA = rank(A, r)
    for l in range(r):
        v = [0] * r
        v[l] = 1
        if rank(A + [tuple(v)], r) == rkA:
            return False, "z%d=0" % (l + 1)
    for l in range(r):
        for m in range(l + 1, r):
            v = [0] * r
            v[l] = 1; v[m] = -1
            if rank(A + [tuple(v)], r) == rkA:
                return False, "z%d=z%d" % (l + 1, m + 1)
    return True, "pass"


ENGINE = "flint"


def _init(engine):
    global ENGINE
    ENGINE = engine


def _work(chunk):
    out = []
    for s in chunk:
        ok, why = affine_test(s, ENGINE)
        out.append((s, ok, why))
    return out


def main(inp, outp, procs=4, engine="flint"):
    import multiprocessing as mp
    from collections import Counter
    t0 = time.time()
    with gzip.open(inp, "rt") as f:
        pats = [ln.strip() for ln in f if ln.strip()]
    chunks = [pats[i:i + 2000] for i in range(0, len(pats), 2000)]
    surv = []
    reasons = Counter()
    byr_in, byr_out = Counter(), Counter()
    with mp.Pool(procs, initializer=_init, initargs=(engine,)) as pool:
        for k, res in enumerate(pool.imap(_work, chunks)):
            for s, ok, why in res:
                r = max(int(ch) for ch in s)
                byr_in[r] += 1
                reasons[why.split(":")[0] if not ok else "pass"] += 1
                if ok:
                    surv.append(s)
                    byr_out[r] += 1
            if k % 25 == 0:
                print(f"  {min((k + 1) * 2000, len(pats))}/{len(pats)} survivors={len(surv)} t={time.time() - t0:.0f}s", flush=True)
    surv.sort()
    with gzip.open(outp, "wt", newline="\n") as f:
        for s in surv:
            f.write(s + "\n")
    print(f"input={len(pats)} survivors={len(surv)} time={time.time() - t0:.1f}s")
    print("byR in :", dict(sorted(byr_in.items())))
    print("byR out:", dict(sorted(byr_out.items())))
    print("reasons:", dict(reasons))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 4,
         sys.argv[4] if len(sys.argv) > 4 else "flint")
