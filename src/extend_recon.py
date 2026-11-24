"""Builder D: exact reconstruction of the finite (F) six-point cores via Engine 2 (Groebner).

Route (primary, independent of Wang's quartic span W):
  1. Saturated ideal J = (I_minors + linear relations, z_a = 1) : (prod of nonzero factors)^inf
     computed with flint Buchberger exactly as in algebra.engine2 (setup copied here because
     engine2 does not return the basis; helpers are imported read-only from src/algebra.py).
  2. Exact multiplication matrices M_x of Q[z,t]/J on the standard monomials (pure-Fraction
     normal forms; G is a reduced GB in degrevlex).
  3. Primitive element lam = sum alpha_k z_k.  f = charpoly(M_lam); factor over Q.
     For each irreducible factor q, over K = Q[u]/(q) compute ker(M_lam^T - u I).
       * The evaluation functionals ev_p of all points p with lam(p) = u lie in this kernel,
         and the kernel is invariant under every M_x^T (the M commute).
       * If dim = 1: z_k = (v . M_k e)/(v . e).
       * else: R_k = M_k^T restricted to the kernel (invariance verified exactly); accept only
         if (R_k - tau_k I)^d = 0 with tau_k = tr R_k / d for every variable k (one point);
         then z_k = tau_k.  Otherwise lam does not separate: try the next alpha.
     Global check: number of distinct complex points = rank of the Hermite trace form
     [Tr(M_{b_i b_j})] must equal sum of deg q.  (Exact; certifies that every point is found.)
  4. Wang Sec. 6.1 reconstruction in K_q (frame p1=(0,0), p2=(1,0), p3=(0,1), metric G):
     z nonzero & distinct, det G != 0, PD at some real root (exact sign via isolation),
     p4..p6 from the 2x2 systems, all 15 distances verified, general position.
"""
from __future__ import annotations

import itertools
import sys
import os
import time
from fractions import Fraction

import flint
from flint import fmpq_mat, fmpq, fmpq_poly, fmpz_mpoly_ctx, fmpz_mpoly_vec

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import algebra as A  # read-only reuse of pattern/poly helpers
from extend_kfield import KField, kernel_over_K, poly_to_list

PAIRS, TRIPLES, QUADS = A.PAIRS, A.TRIPLES, A.QUADS


# --------------------------------------------------------------------------- 1. saturated GB

def saturated_gb(s, limits=(10000, 100000, 100000)):
    c, r = A.parse(s)
    out = {"pattern": s, "r": r}
    ok, why, Arows = A.affine_tests(c, r)
    if not ok:
        out["verdict"] = "EXCLUDED"; out["reason"] = "affine:" + why
        return out
    free, expr = A._lin_substitution(r, Arows)
    a = free[0]
    others = [l for l in free if l != a]
    varmap = {l: k for k, l in enumerate(others)}
    names = tuple(f"z{l+1}" for l in others) + ("t",)
    ctx = fmpz_mpoly_ctx.get(names, "degrevlex")
    conv = lambda p: A._primitive(A._to_mpoly(p, r, expr, a, ctx, varmap))
    gens = [g for g in (conv(m) for m in A.gram_minors(c, r)) if not g.is_zero()]
    facs = []
    for l in range(r):
        facs.append((f"z{l+1}", conv(A.p_var(l, r))))
    for l, t in itertools.combinations(range(r), 2):
        facs.append((f"z{l+1}-z{t+1}", conv(A.p_add(A.p_var(l, r), A.p_var(t, r), -1))))
    for T, f in A.triangle_factors(c, r):
        facs.append(("tri", conv(f)))
    for Q, f in A.quad_dets(c, r):
        facs.append(("det", conv(f)))
    out.update(a=a, free=free, expr=expr, others=others, names=names, c=c)
    for nm, f in facs:
        if f.is_zero():
            out["verdict"] = "EXCLUDED"; out["reason"] = f"factor {nm} == 0"
            return out
    if not gens:
        # no minors: only possible if every z is fixed by linear relations (no free var but a)
        if len(names) == 1:
            out.update(verdict="ZERODIM", G=[ctx.gens()[-1] * 0 + 0], ctx=ctx, nolin=True)
            # J = <t*P - 1> with P constant: one point
            P = ctx.constant(1)
            for nm, f in facs:
                P = P * f
            out["G"] = [ctx.gens()[-1] * P - 1]
            out["ctx"] = ctx
            return out
        out["verdict"] = "POSDIM"
        return out
    G, okb = A._gb(gens, ctx, limits)
    if A._is_unit_basis(G):
        out["verdict"] = "EXCLUDED"; out["reason"] = "GB(I)=1"; return out
    if not okb:
        out["verdict"] = "UNRESOLVED"; out["reason"] = "GB(I) limit"; return out
    P = ctx.constant(1)
    for nm, f in facs:
        if f.is_constant():
            continue
        P = A._nf(P * f, G, ctx)
        if P.is_zero():
            out["verdict"] = "EXCLUDED"; out["reason"] = "NF0"; return out
    tvar = ctx.gens()[-1]
    G2, ok2 = A._gb(G + [tvar * P - 1], ctx, limits)
    if A._is_unit_basis(G2):
        out["verdict"] = "EXCLUDED"; out["reason"] = "GB(sat)=1"; return out
    if not ok2:
        out["verdict"] = "UNRESOLVED"; out["reason"] = "GB(sat) limit"; return out
    nv = len(names)
    lms = [g.monoms()[0] for g in G2]
    pure = set()
    for m in lms:
        nz = [k for k in range(nv) if m[k]]
        if len(nz) == 1:
            pure.add(nz[0])
    if len(pure) < nv:
        out["verdict"] = "POSDIM"; return out
    out.update(verdict="ZERODIM", G=G2, ctx=ctx)
    return out


# --------------------------------------------------------------------------- 2. mult matrices

def mult_matrices(G, ctx):
    nv = len(ctx.names())
    lms = [tuple(g.monoms()[0]) for g in G]

    def divisible(m):
        return any(all(m[k] >= lm[k] for k in range(nv)) for lm in lms)
    std, seen, frontier = [], set(), [tuple([0] * nv)]
    seen.add(frontier[0])
    while frontier:
        m = frontier.pop()
        if divisible(m):
            continue
        std.append(m)
        for k in range(nv):
            mm = tuple(m[i] + (k == i) for i in range(nv))
            if mm not in seen:
                seen.add(mm); frontier.append(mm)
    std.sort()
    idx = {m: i for i, m in enumerate(std)}
    h = len(std)
    Gd = [([tuple(m) for m in g.monoms()], [Fraction(int(x)) for x in g.coeffs()]) for g in G]
    key = lambda m: (sum(m), tuple(-x for x in reversed(m)))  # degrevlex

    def nf(terms):
        terms = dict(terms)
        rem = {}
        while terms:
            lt = max(terms, key=key)
            cf = terms.pop(lt)
            if cf == 0:
                continue
            for gm, gc in Gd:
                l0 = gm[0]
                if all(lt[k] >= l0[k] for k in range(nv)):
                    q = tuple(lt[k] - l0[k] for k in range(nv))
                    s_ = cf / gc[0]
                    for m2, c2 in zip(gm[1:], gc[1:]):
                        mm = tuple(m2[k] + q[k] for k in range(nv))
                        v = terms.get(mm, 0) - s_ * c2
                        if v == 0:
                            terms.pop(mm, None)
                        else:
                            terms[mm] = v
                    break
            else:
                rem[lt] = rem.get(lt, 0) + cf
        return rem

    mats = []
    for k in range(nv):
        M = [[fmpq(0)] * h for _ in range(h)]
        for j, m in enumerate(std):
            mm = tuple(m[i] + (k == i) for i in range(nv))
            rem = nf({mm: Fraction(1)})
            for m2, cf in rem.items():
                M[idx[m2]][j] = fmpq(cf.numerator, cf.denominator)
        mats.append(fmpq_mat(M))
    return std, mats


# --------------------------------------------------------------------------- 3. points

def _to_K_matrix(K, M, shift=None):
    h = M.nrows()
    out = []
    for i in range(h):
        row = []
        for j in range(h):
            x = fmpq_poly([M[i, j]]) if M[i, j] != 0 else K.zero
            if shift is not None and i == j:
                x = x - shift
            row.append(K.red(x))
        out.append(row)
    return out


def _matvec_K(K, MK, v):
    return [K.red(sum((a * b for a, b in zip(row, v)), K.zero)) for row in MK]


def hermite_rank(std, mats):
    """Rank of the trace form = number of distinct complex points of V(J)."""
    nv = len(mats)
    h = len(std)
    I = fmpq_mat([[fmpq(int(i == j)) for j in range(h)] for i in range(h)])
    cache = {}

    def mon_mat(m):
        if m in cache:
            return cache[m]
        R = I
        for k, e in enumerate(m):
            for _ in range(e):
                R = R * mats[k]
        cache[m] = R
        return R
    Ms = [mon_mat(m) for m in std]
    T = fmpq_mat([[ (Ms[i] * Ms[j]).trace() if hasattr(fmpq_mat, 'trace') else sum((Ms[i]*Ms[j])[k, k] for k in range(h)) for j in range(h)] for i in range(h)])
    return T.rank()


ALPHAS = [None, "1s", "lin", "pow3", "pow7", "r1", "r2", "r3", "r4"]


def _alpha(kind, nz):
    import random
    if kind == "1s":
        return [1] * nz
    if kind == "lin":
        return [k + 1 for k in range(nz)]
    if kind == "pow3":
        return [3 ** k for k in range(nz)]
    if kind == "pow7":
        return [(-7) ** k for k in range(nz)]
    rng = random.Random(kind)
    return [rng.randint(-50, 50) or 1 for _ in range(nz)]


def points_of(std, mats, nz):
    """All points of V(J) in the z-variables (first nz of the variables), as a list of
    (q, [z_k in K_q]).  Returns (points, info)."""
    h = len(std)
    e0 = std.index(tuple([0] * len(std[0])))
    npts = hermite_rank(std, mats)
    info = {"h": h, "npoints": npts}
    if nz == 0:
        # only t: one point iff npts==1 (no z variables -> the normalized tuple is fixed)
        info["alpha"] = []
        return ([(fmpq_poly([0, 1]), [], {"kernel_dim": 1})] if npts else []), info
    for kind in ALPHAS[1:]:
        al = _alpha(kind, nz)
        Ml = fmpq_mat(h, h)
        for k in range(nz):
            Ml = Ml + mats[k] * al[k]
        f = Ml.charpoly()
        f = fmpq_poly(f.coeffs()) if not isinstance(f, fmpq_poly) else f
        cont, facs = f.factor()
        qs = [fmpq_poly(q) for q, m in facs]
        if sum(q.degree() for q in qs) != npts:
            continue  # lam does not separate the points (or would collapse) -> next alpha
        pts = []
        good = True
        for q in qs:
            K = KField(q)
            u = fmpq_poly([0, 1])
            MlT = _to_K_matrix(K, Ml.transpose(), shift=u)
            basis, free = kernel_over_K(K, MlT)
            d = len(basis)
            if d == 0:
                good = False; break
            MkT = [_to_K_matrix(K, mats[k].transpose()) for k in range(nz)]
            if d == 1:
                v = basis[0]
                ve = v[e0]
                if ve == 0:
                    good = False; break
                ive = K.inv(ve)
                zs = []
                for k in range(nz):
                    col = [K.red(fmpq_poly([mats[k][i, e0]])) for i in range(h)]
                    zs.append(K.mul(ive, K.red(sum((a * b for a, b in zip(v, col)), K.zero))))
                pts.append((K.q, zs, {"kernel_dim": 1}))
            else:
                zs = []
                for k in range(nz):
                    W = [_matvec_K(K, MkT[k], b) for b in basis]  # images of basis vectors
                    # coordinates in the basis: read off on free coords (basis is reduced)
                    R = [[W[j][free[i]] for j in range(d)] for i in range(d)]
                    # invariance check: W[j] == sum_i R[i][j] basis[i]
                    for j in range(d):
                        recon = [K.red(sum((R[i][j] * basis[i][x] for i in range(d)), K.zero)) for x in range(h)]
                        if any(K.red(recon[x] - W[j][x]) != 0 for x in range(h)):
                            raise AssertionError("kernel not invariant: matrices do not commute?")
                    tau = K.mul(K.red(sum((R[i][i] for i in range(d)), K.zero)), K.inv(K.const(d)))
                    N = [[K.red(R[i][j] - (tau if i == j else 0)) for j in range(d)] for i in range(d)]
                    P = N
                    for _ in range(d - 1):
                        P = [[K.red(sum((P[i][l] * N[l][j] for l in range(d)), K.zero)) for j in range(d)] for i in range(d)]
                    if any(x != 0 for row in P for x in row):
                        good = False; break
                    zs.append(tau)
                if not good:
                    break
                pts.append((K.q, zs, {"kernel_dim": d}))
        if good:
            info["alpha"] = al
            info["alpha_kind"] = kind
            return pts, info
    info["failed"] = "no separating linear form found"
    return None, info


# --------------------------------------------------------------------------- 4. reconstruction

def pdet(M, K):
    n = len(M)
    if n == 1:
        return K.red(M[0][0])
    s = K.zero
    for j in range(n):
        if M[0][j] == 0:
            continue
        minor = [row[:j] + row[j + 1:] for row in M[1:]]
        t = K.mul(M[0][j], pdet(minor, K))
        s = s + t if j % 2 == 0 else s - t
    return K.red(s)


def reconstruct_core(c, r, K, z):
    """z: list of r elements of K (label values).  Returns record dict with status."""
    rec = {}
    if any(x == 0 for x in z) or len({K.key(x) for x in z}) < r:
        rec["status"] = "zero or equal label values"
        return rec
    D = lambda i, j: K.zero if i == j else z[c[i][j]]
    g11, g22 = D(0, 1), D(0, 2)
    g12 = K.red((D(0, 1) + D(0, 2) - D(1, 2)) / 2)
    detG = K.red(g11 * g22 - g12 * g12)
    if detG == 0:
        rec["status"] = "det G = 0"
        return rec
    nreal = K.nreal()
    rec["nreal"] = nreal
    if nreal == 0:
        rec["status"] = "no real root"
        return rec
    pd = [k for k in range(nreal) if K.sign_at(g11, k) > 0 and K.sign_at(detG, k) > 0]
    if not pd:
        rec["status"] = "G not positive definite at any real root"
        return rec
    pts = [(K.zero, K.zero), (K.one, K.zero), (K.zero, K.one)]
    idet = K.inv(K.red(4 * detG))
    for i in range(3, 6):
        b1 = D(0, i) + D(0, 1) - D(1, i)
        b2 = D(0, i) + D(0, 2) - D(2, i)
        x = K.mul(idet, K.red(2 * g22 * b1 - 2 * g12 * b2))
        y = K.mul(idet, K.red(-2 * g12 * b1 + 2 * g11 * b2))
        pts.append((x, y))
    nrm = lambda v: K.red(v[0] * v[0] * g11 + 2 * v[0] * v[1] * g12 + v[1] * v[1] * g22)
    for i, j in PAIRS:
        dv = (pts[i][0] - pts[j][0], pts[i][1] - pts[j][1])
        if K.red(nrm(dv) - D(i, j)) != 0:
            rec["status"] = "distance verification failed"
            return rec
    for i, j, k in TRIPLES:
        d = K.red((pts[j][0] - pts[i][0]) * (pts[k][1] - pts[i][1]) - (pts[j][1] - pts[i][1]) * (pts[k][0] - pts[i][0]))
        if d == 0:
            rec["status"] = "collinear triple"
            return rec
    for Q in QUADS:
        M = [[K.one, pts[i][0], pts[i][1], nrm(pts[i])] for i in Q]
        if pdet(M, K) == 0:
            rec["status"] = "concyclic quadruple"
            return rec
    rec["status"] = "REALIZED"
    rec["pd_roots"] = pd
    rec["G"] = [g11, g12, g22]
    rec["points"] = pts
    return rec


def full_z(info_gb, zfree, K):
    """Label values z_0..z_{r-1} in K from the free-variable values (z_a = 1)."""
    r, a, others, expr = info_gb["r"], info_gb["a"], info_gb["others"], info_gb["expr"]
    val = {a: K.one}
    for k, l in enumerate(others):
        val[l] = zfree[k]
    z = []
    for l in range(r):
        s = K.zero
        for var, coef in expr[l].items():
            s = s + val[var] * fmpq(coef.numerator, coef.denominator)
        z.append(K.red(s))
    return z


def reconstruct_pattern(s, limits=(10000, 100000, 100000)):
    """Full E2 route for one pattern. Returns dict (JSON-able) with per-factor records."""
    t0 = time.perf_counter()
    g = saturated_gb(s, limits)
    out = {"pattern": s, "r": g["r"], "e2": g["verdict"], "factors": []}
    if g["verdict"] != "ZERODIM":
        out["reason"] = g.get("reason")
        out["t"] = time.perf_counter() - t0
        return out
    c = g["c"]
    out["a"] = g["a"] + 1
    out["others"] = [l + 1 for l in g["others"]]
    std, mats = mult_matrices(g["G"], g["ctx"])
    nz = len(g["names"]) - 1
    pts, info = points_of(std, mats, nz)
    out.update({k: v for k, v in info.items() if k in ("h", "npoints", "alpha", "alpha_kind", "failed")})
    if pts is None:
        out["status"] = "UNRESOLVED (no separating form)"
        out["t"] = time.perf_counter() - t0
        return out
    for q, zf, meta in pts:
        K = KField(q)
        z = full_z(g, zf, K)
        rec = reconstruct_core(c, g["r"], K, z)
        fr = {"q": poly_to_list(K.q), "deg": K.deg, "status": rec["status"],
              "nreal": rec.get("nreal"), "kernel_dim": meta.get("kernel_dim"),
              "z": [poly_to_list(x) for x in z]}
        if rec["status"] == "REALIZED":
            fr["pd_roots"] = rec["pd_roots"]
            fr["roots"] = [K.root_str(k) for k in rec["pd_roots"]]
            fr["G"] = [poly_to_list(x) for x in rec["G"]]
            fr["points"] = [[poly_to_list(p[0]), poly_to_list(p[1])] for p in rec["points"]]
        out["factors"].append(fr)
    out["t"] = time.perf_counter() - t0
    return out


if __name__ == "__main__":
    import json
    for s in sys.argv[1:]:
        print(json.dumps(reconstruct_pattern(s)))
