"""Builder A: algebraic elimination / certificate engine for 6-point distance patterns
(Erdos #217, re-implementation of Wang arXiv:2611.03688 Sec. 4 + an independent Groebner route).

Pattern input: 15 label digits over pairs (1,2),(1,3),...,(1,6),(2,3),...,(5,6); labels 1..r, r<=6.
Internally points are 0..5 and labels 0..r-1; z_l is the squared distance of label l.

ENGINE 1  (engine1(pattern, ...))   -- Wang's quartic-span test, exact.
  * affine_tests()    Lemma 2.3 (rows of R: parallelograms + centroids; row-space support <=3;
                       linear distance equations Az=0 with zero / equal kernel columns).
  * W0 = span{ z_l * (3x3 minors of L) } + span{ (rows of A) * cubic monomials }  in V = Q[z]_4,
    L_ij = D_1i + D_1j - D_ij (2<=i,j<=6), all 100 minors (55 distinct up to transpose).
  * division fixpoint (Lemma 4.1) by f in {z_l, z_l - z_t, 4D_ijD_ik-(D_ij+D_ik-D_jk)^2}.
  * exclusion (Lemma 4.2): W contains a degree-4 monomial or det of the 4x4 D of 4 points.
  * finite criterion V = W + z_a Q[z]_3, F/U split exactly as Procedure 4.4.
  Linear algebra: exact over Q via flint fmpz_mat (fraction-free rref / integer nullspace),
  optionally preceded by a run modulo a 62-bit prime (nmod_mat) used ONLY as a hint:
  every EXCLUDED verdict is re-derived over Q (by replaying the mod-p enlargement trace over
  Q, falling back to the full Q fixpoint), and every F/U verdict is computed over Q.

ENGINE 2  (engine2(pattern, ...))   -- independent route: Groebner bases over Q
  (python-flint fmpz_mpoly_vec.buchberger_naive, degrevlex).  Ideal: 3x3 minors of L,
  linear relations eliminated by exact substitution, one free z normalised to 1,
  saturated (Rabinowitsch, one variable t) by the product of all nonzero factors
  z_l, z_l-z_t, triangle factors, 4-point det(D).  The product is reduced modulo GB(I)
  factor by factor (P <- NF(P*f)), which keeps it small and does not change the ideal
  I + <1 - tP>.  Outcomes: EXCLUDED (GB = {1}, or NF(prod)=0, or a factor is identically 0
  after the linear substitution), ZERODIM (+ eliminants), POSDIM, or UNRESOLVED (limit/timeout).
"""
from __future__ import annotations

import itertools
import time
from fractions import Fraction
from math import comb, gcd

import flint
from flint import fmpz_mat, nmod_mat, fmpz_mpoly_ctx, fmpz_mpoly_vec, fmpq_mat, fmpq

PAIRS = [(i, j) for i in range(6) for j in range(i + 1, 6)]
TRIPLES = list(itertools.combinations(range(6), 3))
QUADS = list(itertools.combinations(range(6), 4))
PRIME = 4611686018427387847  # 2^62 - 57, prime

# ----------------------------------------------------------------------------- patterns


def parse(s):
    """Return (c, r): c 6x6 label matrix (0-based labels, -1 diagonal), r = #labels."""
    digs = [int(ch) for ch in str(s).strip() if ch.isdigit()]
    if len(digs) != 15:
        raise ValueError(f"pattern needs 15 labels: {s!r}")
    c = [[-1] * 6 for _ in range(6)]
    for (i, j), d in zip(PAIRS, digs):
        c[i][j] = c[j][i] = d - 1
    labs = sorted(set(digs))
    if labs != list(range(1, len(labs) + 1)):
        raise ValueError(f"labels must be 1..r: {s!r}")
    return c, len(labs)


def normalize(labels):
    """First-occurrence renaming of a 15-sequence of hashable labels -> string."""
    m = {}
    out = []
    for x in labels:
        if x not in m:
            m[x] = len(m) + 1
        out.append(str(m[x]))
    return "".join(out)


# ----------------------------------------------------------------------------- polynomials
# Sparse polys as dict {exponent tuple (len r): int}.


def p_add(a, b, s=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + s * v
        if out[k] == 0:
            del out[k]
    return out


def p_mul(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = tuple(x + y for x, y in zip(ka, kb))
            out[k] = out.get(k, 0) + va * vb
    return {k: v for k, v in out.items() if v}


def p_var(l, r):
    e = [0] * r
    e[l] = 1
    return {tuple(e): 1}


def p_det(M):
    n = len(M)
    if n == 1:
        return M[0][0]
    out = {}
    for j in range(n):
        if not M[0][j]:
            continue
        minor = [row[:j] + row[j + 1:] for row in M[1:]]
        term = p_mul(M[0][j], p_det(minor))
        out = p_add(out, term, 1 if j % 2 == 0 else -1)
    return out


def monomials(r, d):
    """All exponent tuples of total degree d in r variables, lex-descending."""
    if r == 0:
        return [()] if d == 0 else []
    out = []
    for e0 in range(d, -1, -1):
        for rest in monomials(r - 1, d - e0):
            out.append((e0,) + rest)
    return out


def Dpoly(c, r, i, j):
    return {} if i == j else p_var(c[i][j], r)


def gram_minors(c, r):
    """All 3x3 minors of L (5x5, L_ij = D_1i + D_1j - D_ij), deduplicated, as polys."""
    L = [[p_add(p_add(Dpoly(c, r, 0, i), Dpoly(c, r, 0, j)), Dpoly(c, r, i, j), -1)
          for j in range(1, 6)] for i in range(1, 6)]
    seen = set()
    out = []
    for R in itertools.combinations(range(5), 3):
        for C in itertools.combinations(range(5), 3):
            if C < R:
                continue  # transpose symmetry: minor(R,C) = minor(C,R)
            m = p_det([[L[i][j] for j in C] for i in R])
            key = tuple(sorted(m.items()))
            if m and key not in seen:
                seen.add(key)
                out.append(m)
    return out


def triangle_factors(c, r):
    out, seen = [], set()
    for i, j, k in TRIPLES:
        a, b, d = Dpoly(c, r, i, j), Dpoly(c, r, i, k), Dpoly(c, r, j, k)
        s = p_add(p_add(a, b), d, -1)
        f = p_add(p_mul(p_mul({(0,) * r: 4}, a), b), p_mul(s, s), -1)
        key = tuple(sorted(f.items()))
        if key not in seen:
            seen.add(key)
            out.append(((i, j, k), f))
    return out


def quad_dets(c, r):
    out, seen = [], set()
    for Q in QUADS:
        f = p_det([[Dpoly(c, r, i, j) for j in Q] for i in Q])
        key = tuple(sorted(f.items()))
        if key not in seen:
            seen.add(key)
            out.append((Q, f))
    return out


# ----------------------------------------------------------------------------- local tests


def local_ok(c):
    """Lemma 2.1 (i)-(iii) and Lemma 2.2 on a full 6-point pattern (necessary conditions)."""
    n = len(c)
    for a in range(n):
        cnt = {}
        for x in range(n):
            if x != a:
                cnt[c[a][x]] = cnt.get(c[a][x], 0) + 1
        if max(cnt.values()) > 3:
            return False
    for a in range(n):
        for b in range(a + 1, n):
            cnt = {}
            eq = 0
            for x in range(n):
                if x in (a, b):
                    continue
                k = (c[a][x], c[b][x])
                cnt[k] = cnt.get(k, 0) + 1
                if cnt[k] > 2:
                    return False
                if c[a][x] == c[b][x]:
                    eq += 1
            if eq > 2:
                return False
    for a, b, cc, d in itertools.combinations(range(n), 4):
        if c[a][b] == c[cc][d] and c[a][cc] == c[b][d] and c[a][d] == c[b][cc]:
            return False
    return True


# ----------------------------------------------------------------------------- affine (Lemma 2.3)


def affine_rows(c):
    n = len(c)
    rows = set()
    for a, b, cc, d in itertools.permutations(range(n), 4):
        if c[a][b] == c[cc][d] and c[a][cc] == c[b][d]:
            v = [0] * n
            v[a] += 1; v[d] += 1; v[b] -= 1; v[cc] -= 1
            rows.add(tuple(v))
    for a, b, cc in itertools.combinations(range(n), 3):
        if c[a][b] == c[a][cc] == c[b][cc]:
            for d in range(n):
                if d not in (a, b, cc) and c[a][d] == c[b][d] == c[cc][d]:
                    v = [0] * n
                    v[a] += 1; v[b] += 1; v[cc] += 1; v[d] -= 3
                    rows.add(tuple(v))
    return sorted(rows)


def _rank(rows, ncols):
    if not rows:
        return 0
    return fmpz_mat([list(x) for x in rows]).rank()


def affine_tests(c, r):
    """Returns (ok: bool, reason, A_rows) where A_rows (integer, length r) span the linear
    relations sum_i rho_i (D_ij - D_i1) = 0 (one per R-row and j)."""
    n = len(c)
    R = affine_rows(c)
    rk = _rank(R, n)
    if R:
        for I in itertools.combinations(range(n), 3):
            keep = [k for k in range(n) if k not in I]
            if _rank([[row[k] for k in keep] for row in R], len(keep)) < rk:
                return False, f"affine-support{I}", []
    A = set()
    for rho in R:
        for j in range(1, n):
            v = [0] * r
            for i in range(n):
                if rho[i] == 0:
                    continue
                if i != j:
                    v[c[i][j]] += rho[i]
                if i != 0:
                    v[c[i][0]] -= rho[i]
            if any(v):
                g = 0
                for x in v:
                    g = gcd(g, x)
                v = [x // g for x in v]
                if next(x for x in v if x) < 0:
                    v = [-x for x in v]
                A.add(tuple(v))
    A = sorted(A)
    if A:
        rA = _rank(A, r)
        for l in range(r):
            e = [0] * r; e[l] = 1
            if _rank(A + [tuple(e)], r) == rA:
                return False, f"lin-zero z{l+1}", A
        for l, t in itertools.combinations(range(r), 2):
            e = [0] * r; e[l] = 1; e[t] = -1
            if _rank(A + [tuple(e)], r) == rA:
                return False, f"lin-equal z{l+1}=z{t+1}", A
        # reduce A to a basis (rref rows)
        M, den, rk2 = fmpz_mat([list(x) for x in A]).rref()
        A = [tuple(int(x) for x in M.tolist()[i]) for i in range(rk2)]
    return True, "", A


# ----------------------------------------------------------------------------- linear algebra


class _Field:
    """Exact field ops on integer row lists. kind 'Q' -> fmpz_mat; 'P' -> nmod_mat mod PRIME."""

    def __init__(self, kind):
        self.kind = kind

    def mat(self, rows, ncols):
        if self.kind == "Q":
            return fmpz_mat(rows) if rows else fmpz_mat(0, ncols)
        return nmod_mat(rows, PRIME) if rows else nmod_mat(0, ncols, PRIME)

    def rref_rows(self, rows, ncols):
        """Return (basis rows as int lists, pivot list)."""
        if not rows:
            return [], []
        M = self.mat(rows, ncols)
        res = M.rref()
        R, rk = res[0], res[-1]
        L = R.tolist()
        out, piv = [], []
        for i in range(rk):
            row = [int(x) for x in L[i]]
            if self.kind == "Q":
                g = 0
                for x in row:
                    if x:
                        g = gcd(g, x)
                if g > 1:
                    row = [x // g for x in row]
            out.append(row)
            piv.append(next(k for k, x in enumerate(row) if x))
        return out, piv

    def left_kernel(self, rows, ncols):
        """Basis of {c : c * M = 0} for M = rows (as int lists)."""
        M = self.mat(rows, ncols).transpose()
        X, nul = M.nullspace()
        L = X.tolist()
        return [[int(L[i][k]) for i in range(len(rows))] for k in range(nul)]

    def rank(self, rows, ncols):
        if not rows:
            return 0
        return self.mat(rows, ncols).rank()


QF = _Field("Q")
PF = _Field("P")

# ----------------------------------------------------------------------------- engine 1


class _Setup:
    """Per-pattern data shared by the Q and mod-p runs."""

    def __init__(self, c, r, A):
        self.c, self.r, self.A = c, r, A
        self.m4 = monomials(r, 4)
        self.i4 = {m: k for k, m in enumerate(self.m4)}
        self.n = len(self.m4)
        self.m3 = monomials(r, 3)
        self.m2 = monomials(r, 2)
        self.m1 = monomials(r, 1)
        self.minors = gram_minors(c, r)
        self.tri = triangle_factors(c, r)
        self.dets = quad_dets(c, r)
        # factors (name, deg, poly)
        fac = []
        for l in range(r):
            fac.append((f"z{l+1}", 1, p_var(l, r)))
        for l, t in itertools.combinations(range(r), 2):
            fac.append((f"z{l+1}-z{t+1}", 1, p_add(p_var(l, r), p_var(t, r), -1)))
        for (i, j, k), f in self.tri:
            fac.append((f"tri{i+1}{j+1}{k+1}", 2, f))
        self.factors = fac
        # initial generators
        gens = []
        for m in self.minors:
            for l in range(r):
                gens.append(self.vec(p_mul(m, p_var(l, r))))
        for a in A:
            lin = {}
            for l in range(r):
                if a[l]:
                    lin[tuple(1 if k == l else 0 for k in range(r))] = a[l]
            for mon in self.m3:
                gens.append(self.vec(p_mul(lin, {mon: 1})))
        self.gens = gens
        self.det_vecs = [(Q, self.vec(f)) for Q, f in self.dets]
        # multiplication data: for factor f of degree e, the matrix F with rows f*m_i
        self.Fmats = []
        for name, e, f in fac:
            base = self.m3 if e == 1 else self.m2
            Frows = [self.vec(p_mul(f, {m: 1})) for m in base]
            self.Fmats.append((base, Frows))

    def vec(self, p):
        v = [0] * self.n
        for k, x in p.items():
            v[self.i4[k]] = x
        return v


def _reduce(v, W, piv, field):
    """Reduce row v modulo rref basis W (integer rows; for Q the pivot entry may be != 1)."""
    v = list(v)
    if field.kind == "P":
        v = [x % PRIME for x in v]
        for row, p in zip(W, piv):
            if v[p]:
                f = v[p]
                v = [(x - f * y) % PRIME for x, y in zip(v, row)]
        return v
    for row, p in zip(W, piv):
        if v[p]:
            a, b = row[p], v[p]
            v = [a * x - b * y for x, y in zip(v, row)]
            g = 0
            for x in v:
                if x:
                    g = gcd(g, x)
            if g > 1:
                v = [x // g for x in v]
    return v


def _excluded(S, W, piv, field):
    """Lemma 4.2 on the current W. Returns reason or None."""
    for row in W:
        nz = [k for k, x in enumerate(row) if x]
        if len(nz) == 1:
            return ("monomial", S.m4[nz[0]])
    for Q, dv in S.det_vecs:
        if not any(_reduce(dv, W, piv, field)):
            return ("det4", Q)
    return None


def _division_step(S, fi, W, piv, field):
    """Lemma 4.1 for factor index fi. Returns list of new rows (possibly empty)."""
    name, e, f = S.factors[fi]
    base, Frows = S.Fmats[fi]
    nq = len(base)
    rows = [list(x) for x in Frows] + W
    K = field.left_kernel(rows, S.n)
    Cs = [k[:nq] for k in K if any(k[:nq])]
    if not Cs:
        return []
    Cs, _ = field.rref_rows(Cs, nq)
    mult = S.m1 if e == 1 else S.m2
    new = []
    for cvec in Cs:
        for mp in mult:
            v = [0] * S.n
            for coef, mon in zip(cvec, base):
                if coef:
                    v[S.i4[tuple(x + y for x, y in zip(mon, mp))]] += coef
            new.append(v)
    return new


def quartic_fixpoint(S, field, early_exit=True, replay=None, trace=None):
    """Run the division fixpoint. replay: optional list of factor indices to apply in order
    (only those, then stop).  Returns (status, W, piv, info) with status 'EXCLUDED' or 'RETAINED'
    ('PARTIAL' when a replay ends without exclusion)."""
    W, piv = field.rref_rows(S.gens, S.n)
    ex = _excluded(S, W, piv, field)
    if ex:
        return "EXCLUDED", W, piv, ex
    order = replay
    passes = 0
    while True:
        changed = False
        seq = order if order is not None else range(len(S.factors))
        for fi in seq:
            new = _division_step(S, fi, W, piv, field)
            if not new:
                continue
            W2, piv2 = field.rref_rows(W + new, S.n)
            if len(W2) > len(W):
                if trace is not None:
                    trace.append((fi, len(W), len(W2)))
                W, piv = W2, piv2
                changed = True
                if early_exit:
                    ex = _excluded(S, W, piv, field)
                    if ex:
                        return "EXCLUDED", W, piv, ex
        passes += 1
        if order is not None:
            ex = _excluded(S, W, piv, field)
            if ex:
                return "EXCLUDED", W, piv, ex
            return "PARTIAL", W, piv, None
        if not changed:
            break
    ex = _excluded(S, W, piv, field)
    if ex:
        return "EXCLUDED", W, piv, ex
    return "RETAINED", W, piv, {"passes": passes}


def finite_labels(S, W, field):
    """Labels a (0-based) with V = W + z_a Q[z]_3."""
    out = []
    for a in range(S.r):
        cols = [k for k, m in enumerate(S.m4) if m[a] == 0]
        sub = [[row[k] for k in cols] for row in W]
        if field.rank(sub, len(cols)) == len(cols):
            out.append(a)
    return out


def engine1(s, use_modp=True, check_affine=True, check_local=False):
    """Classify a pattern.  Returns dict with keys:
    verdict: 'EXCLUDED' | 'F' | 'U'
    stage:   'local' | 'affine' | 'quartic' | None
    reason, exact (True: the verdict was derived over Q), cert (for quartic exclusions:
    {'trace': [(factor name, dim before, dim after)...], 'target': ...}), h, finite (labels),
    a (smallest finite label, 1-based), dimW, r, t (seconds)."""
    t0 = time.perf_counter()
    c, r = parse(s)
    res = {"pattern": s, "r": r}
    if check_local and not local_ok(c):
        res.update(verdict="EXCLUDED", stage="local", reason="local", exact=True)
        res["t"] = time.perf_counter() - t0
        return res
    if check_affine:
        ok, why, A = affine_tests(c, r)
        if not ok:
            res.update(verdict="EXCLUDED", stage="affine", reason=why, exact=True)
            res["t"] = time.perf_counter() - t0
            return res
    else:
        A = affine_tests(c, r)[2]
    S = _Setup(c, r, A)
    res["A"] = [list(a) for a in A]
    qres = None
    if use_modp:
        tr = []
        st, W, piv, info = quartic_fixpoint(S, PF, trace=tr)
        res["modp"] = st
        if st == "EXCLUDED":
            tq = []
            st2, Wq, pq, info2 = quartic_fixpoint(S, QF, replay=[x[0] for x in tr], trace=tq)
            if st2 == "EXCLUDED":
                qres = (st2, Wq, pq, info2, tq)
            # else: fall through to a full Q run
    if qres is None:
        tq = []
        st2, Wq, pq, info2 = quartic_fixpoint(S, QF, trace=tq)
        qres = (st2, Wq, pq, info2, tq)
        if use_modp and res["modp"] != st2:
            res["modp_discrepancy"] = True
    st2, Wq, pq, info2, tq = qres
    res["exact"] = True
    res["dimW"] = len(Wq)
    res["dimV"] = S.n
    if st2 == "EXCLUDED":
        res.update(verdict="EXCLUDED", stage="quartic", reason=f"{info2[0]}:{info2[1]}")
        res["cert"] = {"trace": [(S.factors[fi][0], a, b) for fi, a, b in tq],
                       "target": [info2[0], list(info2[1])]}
    else:
        fin = finite_labels(S, Wq, QF)
        res["finite"] = [a + 1 for a in fin]
        res["h"] = S.n - len(Wq)
        if r in (5, 6) and fin:
            res.update(verdict="F", a=fin[0] + 1)
        else:
            res.update(verdict="U")
        res["stage"] = None
        res["_W"] = Wq
        res["_S"] = S
    res["t"] = time.perf_counter() - t0
    return res


def verify_certificate(s, cert):
    """Re-derive a quartic exclusion over Q from its trace (factor names in order)."""
    c, r = parse(s)
    ok, why, A = affine_tests(c, r)
    S = _Setup(c, r, A)
    names = [f[0] for f in S.factors]
    order = [names.index(nm) for nm, _, _ in cert["trace"]]
    st, W, piv, info = quartic_fixpoint(S, QF, replay=order)
    return st == "EXCLUDED"


# ----------------------------------------------------------------------------- engine 2


def _lin_substitution(r, A):
    """Solve A z = 0 (rows integer). Returns (free, expr, a) where expr[l] is a dict
    {free var index: Fraction} giving z_l as a combination of free vars, and a = the free
    label normalised to 1.  None if A forces z = 0."""
    if not A:
        free = list(range(r))
        expr = [{l: Fraction(1)} for l in range(r)]
        return free, expr
    M = fmpq_mat([[fmpq(x) for x in row] for row in A])
    R, rk = M.rref()
    L = R.tolist()
    piv = []
    for i in range(rk):
        piv.append(next(k for k in range(r) if L[i][k] != 0))
    free = [l for l in range(r) if l not in piv]
    expr = [None] * r
    for l in free:
        expr[l] = {l: Fraction(1)}
    for i, p in enumerate(piv):
        d = {}
        for l in free:
            x = L[i][l]
            if x != 0:
                d[l] = -Fraction(int(x.p), int(x.q))
        expr[p] = d
    return free, expr


def _to_mpoly(p, r, expr, a, ctx, varmap):
    """Substitute z_l -> expr[l] (a combination of free labels), z_a -> 1, into the integer
    poly p; returns an integer mpoly (scaled by a positive integer, which does not change the
    ideal or the zero set)."""
    gens = ctx.gens()
    zl = []
    for l in range(r):
        d = expr[l]
        den = 1
        for v in d.values():
            den = den * v.denominator // gcd(den, v.denominator)
        f = ctx.constant(0)
        for var, coef in d.items():
            ci = int(coef * den)
            f = f + (ci if var == a else ci * gens[varmap[var]])
        zl.append((f, den))
    D = 1
    for mon in p:
        dd = 1
        for l, e in enumerate(mon):
            dd *= zl[l][1] ** e
        D = D * dd // gcd(D, dd)
    out = ctx.constant(0)
    for mon, coef in p.items():
        term = ctx.constant(coef)
        dd = 1
        for l, e in enumerate(mon):
            if e:
                term = term * zl[l][0] ** e
                dd *= zl[l][1] ** e
        out = out + term * (D // dd)
    return out


def _primitive(f):
    if f.is_zero():
        return f
    cont = f.content()
    g = f / cont if cont != 1 else f
    if g.leading_coefficient() < 0:
        g = -g
    return g


def _is_unit_basis(G):
    return any(g.is_constant() and not g.is_zero() for g in G)


def _gb(polys, ctx, limits):
    vec = fmpz_mpoly_vec(polys, ctx)
    G, ok = vec.buchberger_naive(limits=limits)
    G = list(G)
    if _is_unit_basis(G):
        return [ctx.constant(1)], True  # {1} is a GB of the unit ideal
    if not ok:
        return G, False
    G = list(fmpz_mpoly_vec(G, ctx).autoreduction())
    return G, True


def _nf(f, G, ctx):
    if f.is_zero():
        return f
    return f.reduction_primitive_part(fmpz_mpoly_vec(G, ctx))


def engine2(s, limits=(400, 4000, 4096), eliminants=True, check_affine=True):
    """Groebner route. Returns dict: verdict in {'EXCLUDED','ZERODIM','POSDIM','UNRESOLVED'},
    reason, cert (for exclusions), eliminant info for ZERODIM, t."""
    t0 = time.perf_counter()
    c, r = parse(s)
    res = {"pattern": s, "r": r}
    ok, why, A = affine_tests(c, r)
    if check_affine and not ok:
        res.update(verdict="EXCLUDED", reason="affine:" + why, cert={"kind": "affine", "why": why})
        res["t"] = time.perf_counter() - t0
        return res
    sub = _lin_substitution(r, A)
    free, expr = sub
    a = free[0]  # normalised label (any free label; it is a positive z)
    others = [l for l in free if l != a]
    varmap = {l: k for k, l in enumerate(others)}
    names = tuple(f"z{l+1}" for l in others) + ("t",)
    ctx = fmpz_mpoly_ctx.get(names, "degrevlex")
    res["vars"] = list(names[:-1])
    res["a"] = a + 1
    conv = lambda p: _primitive(_to_mpoly(p, r, expr, a, ctx, varmap))
    gens = []
    for m in gram_minors(c, r):
        g = conv(m)
        if not g.is_zero():
            gens.append(g)
    facs = []
    for l in range(r):
        facs.append((f"z{l+1}", conv(p_var(l, r))))
    for l, t in itertools.combinations(range(r), 2):
        facs.append((f"z{l+1}-z{t+1}", conv(p_add(p_var(l, r), p_var(t, r), -1))))
    for T, f in triangle_factors(c, r):
        facs.append(("tri" + "".join(str(x + 1) for x in T), conv(f)))
    for Q, f in quad_dets(c, r):
        facs.append(("det" + "".join(str(x + 1) for x in Q), conv(f)))
    for nm, f in facs:
        if f.is_zero():
            res.update(verdict="EXCLUDED", reason=f"factor {nm} == 0 under linear relations",
                       cert={"kind": "factor-zero", "factor": nm})
            res["t"] = time.perf_counter() - t0
            return res
    if not gens:
        if len(names) == 1:  # every z fixed by the linear relations; all factors are nonzero consts
            res.update(verdict="ZERODIM", reason="z determined by linear relations", h=1,
                       eliminants={})
        else:
            res.update(verdict="POSDIM", reason="no equations")
        res["t"] = time.perf_counter() - t0
        return res
    G, ok = _gb(gens, ctx, limits)
    if _is_unit_basis(G):
        res.update(verdict="EXCLUDED", reason="GB(I)={1}", cert={"kind": "GB1", "stage": "I"})
        res["t"] = time.perf_counter() - t0
        return res
    if not ok:
        res.update(verdict="UNRESOLVED", reason="GB(I) limit")
        res["t"] = time.perf_counter() - t0
        return res
    P = ctx.constant(1)
    for nm, f in facs:
        if f.is_constant():
            continue
        P = _nf(P * f, G, ctx)
        if P.is_zero():
            res.update(verdict="EXCLUDED", reason=f"prod of factors in I (zero after {nm})",
                       cert={"kind": "NF0", "upto": nm})
            res["t"] = time.perf_counter() - t0
            return res
    tvar = ctx.gens()[-1]
    G2, ok2 = _gb(G + [tvar * P - 1], ctx, limits)
    if _is_unit_basis(G2):
        res.update(verdict="EXCLUDED", reason="GB(I+<tP-1>)={1}", cert={"kind": "GB1", "stage": "sat"})
        res["t"] = time.perf_counter() - t0
        return res
    if not ok2:
        res.update(verdict="UNRESOLVED", reason="GB(sat) limit")
        res["t"] = time.perf_counter() - t0
        return res
    nv = len(names)
    lms = [g.monoms()[0] for g in G2]
    pure = set()
    for m in lms:
        nz = [k for k in range(nv) if m[k]]
        if len(nz) == 1:
            pure.add(nz[0])
    if len(pure) < nv:
        res.update(verdict="POSDIM", reason="saturated ideal positive-dimensional")
        res["t"] = time.perf_counter() - t0
        return res
    res.update(verdict="ZERODIM")
    if eliminants:
        info = _quotient_info(G2, ctx, lms, nv)
        res.update(info)
    res["t"] = time.perf_counter() - t0
    return res


def _quotient_info(G, ctx, lms, nv):
    """Standard monomials, and minimal polynomials of each z variable and of a generic linear
    form in Q[z,t]/J (= eliminants of the saturated ideal)."""
    def divisible(m):
        return any(all(m[k] >= lm[k] for k in range(nv)) for lm in lms)
    std = []
    frontier = [tuple([0] * nv)]
    seen = set(frontier)
    while frontier:
        m = frontier.pop()
        if divisible(m):
            continue
        std.append(m)
        for k in range(nv):
            mm = tuple(m[i] + (1 if i == k else 0) for i in range(nv))
            if mm not in seen:
                seen.add(mm)
                frontier.append(mm)
    std.sort()
    idx = {m: i for i, m in enumerate(std)}
    h = len(std)
    gens = ctx.gens()
    vecG = fmpz_mpoly_vec(G, ctx)

    def mono(m):
        f = ctx.constant(1)
        for k, e in enumerate(m):
            if e:
                f = f * gens[k] ** e
        return f

    def mult_matrix(g):
        cols = []
        for m in std:
            f = g * mono(m)
            # exact normal form over Q: reduction_primitive_part loses the scale, so recover it
            # by reducing the scaled poly and solving: use divmod-like approach via scaling.
            nf = f.reduction_primitive_part(vecG)
            cols.append((f, nf))
        return cols

    # exact NF with scale: compute nf of f and of (f + X) where X is a fresh standard monomial
    # combination is awkward; instead do our own reduction over Q.
    out = {"h": h}
    Gq = []
    for g in G:
        Gq.append(g)

    def exact_nf(f):
        # pure-Python multivariate division over Q in degrevlex (G is a GB) -> dict std->Fraction
        terms = {tuple(m): Fraction(int(cf)) for m, cf in zip(f.monoms(), f.coeffs())}
        rem = {}
        Gd = [([tuple(m) for m in g.monoms()], [int(x) for x in g.coeffs()]) for g in G]
        key = lambda m: (sum(m), tuple(-x for x in reversed(m)))
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
                        terms[mm] = terms.get(mm, 0) - s_ * c2
                        if terms[mm] == 0:
                            del terms[mm]
                    break
            else:
                rem[lt] = rem.get(lt, 0) + cf
        return rem

    def mat_of(g):
        M = [[fmpq(0)] * h for _ in range(h)]
        for j, m in enumerate(std):
            rem = exact_nf(g * mono(m))
            for mm, cf in rem.items():
                M[idx[mm]][j] = fmpq(cf.numerator, cf.denominator)
        return fmpq_mat(M)

    elim = {}
    for k in range(nv - 1):
        Mk = mat_of(gens[k])
        mp = Mk.minpoly()
        elim[str(ctx.names()[k])] = str(mp)
    out["std_monomials"] = h
    out["eliminants"] = elim
    return out


# ----------------------------------------------------------------------------- reconstruction (F)


def f_univariate(res):
    """Wang Sec. 6 for an F result of engine1 (exact over Q). Returns dict with f(u) factors and
    g_j, or None if the M-choice of App. A.2 fails."""
    S, W = res["_S"], res["_W"]
    a = res["a"] - 1
    r, n = S.r, S.n
    Wr, piv = QF.rref_rows(W, n)
    # basis of V/W from monomials divisible by z_a: order by degree after normalisation, then lex
    cand = [m for m in S.m4 if m[a] > 0]
    cand.sort(key=lambda m: (4 - m[a], tuple(-x for x in m)))
    basis = []
    cur = [list(x) for x in Wr]
    rk = len(cur)
    for m in cand:
        v = [0] * n; v[S.i4[m]] = 1
        if QF.rank(cur + [v], n) > rk:
            cur.append(v); rk += 1; basis.append(m)
    h = len(basis)
    assert rk == n
    # coordinates: solve [W; B] -> express e_m
    Bmat = fmpq_mat([[fmpq(x) for x in row] for row in (cur)])  # n x n invertible
    Binv = Bmat.transpose().inv()

    def coords(m4):
        v = fmpq_mat([[fmpq(1 if k == S.i4[m4] else 0)] for k in range(n)])
        sol = Binv * v
        return [sol[len(Wr) + i, 0] for i in range(h)]

    def col_of(mon):  # mon exponent in r vars of degree <=4 (z_a may appear)
        d = sum(mon)
        m4 = tuple(x + (4 - d if k == a else 0) for k, x in enumerate(mon))
        return coords(m4)

    e = col_of(tuple([0] * r))
    bks = [tuple(x - (1 if k == a else 0) * 0 for k, x in enumerate(m)) for m in basis]
    # b_k = B_k with z_a = 1 : exponent with a-component zeroed
    bks = [tuple(0 if k == a else x for k, x in enumerate(m)) for m in basis]
    Ms = {}
    for j in range(r):
        if j == a:
            continue
        cols = [col_of(tuple(x + (1 if k == j else 0) for k, x in enumerate(b))) for b in bks]
        Ms[j] = fmpq_mat([[cols[kk][i] for kk in range(h)] for i in range(h)])
    out = {"h": h, "commute": all((Ms[i] * Ms[j] - Ms[j] * Ms[i]).is_zero() if hasattr(fmpq_mat, 'is_zero') else (Ms[i] * Ms[j] == Ms[j] * Ms[i]) for i in Ms for j in Ms)}
    ev = fmpq_mat([[x] for x in e])
    for cc in (0, 1, 2):
        M = fmpq_mat(h, h)
        for j in sorted(Ms):
            M = M * cc + Ms[j]
        f = M.charpoly()
        fz = flint.fmpq_poly(f.coeffs()) if not isinstance(f, flint.fmpq_poly) else f
        if flint.fmpq_poly.gcd(fz, fz.derivative()).degree() > 0:
            continue
        cols = [ev]
        for _ in range(h - 1):
            cols.append(M * cols[-1])
        T = fmpq_mat([[cols[k][i, 0] for k in range(h)] for i in range(h)])
        if T.det() == 0:
            continue
        Tinv = T.inv()
        gs = {}
        okg = True
        for j in Ms:
            coef = Tinv * (Ms[j] * ev)
            g = flint.fmpq_poly([coef[i, 0] for i in range(h)])
            # check g(M) == M_j
            acc = fmpq_mat(h, h)
            P_ = fmpq_mat([[fmpq(1 if i == k else 0) for k in range(h)] for i in range(h)])
            for i in range(h):
                acc = acc + P_ * coef[i, 0]
                P_ = P_ * M
            if acc != Ms[j]:
                okg = False
                break
            gs[j] = g
        if not okg:
            continue
        out.update(c=cc, f=fz, g=gs, factors=fz.factor())
        return out
    out["failed"] = True
    return out


# ----------------------------------------------------------------------------- core reconstruction


class NF:
    """Arithmetic in K = Q[u]/(q), q irreducible fmpq_poly."""

    def __init__(self, q):
        self.q = q

    def red(self, a):
        return a % self.q

    def mul(self, a, b):
        return (a * b) % self.q

    def inv(self, a):
        a = a % self.q
        if a == 0:
            raise ZeroDivisionError
        g, s, t = flint.fmpq_poly.xgcd(a, self.q)
        return (s / g[0]) % self.q if g.degree() == 0 else None


def _sturm_count(q):
    """Exact number of distinct real roots of q (Sturm sequence over Q)."""
    seq = [q, q.derivative()]
    while seq[-1].degree() > 0:
        rmd = -(seq[-2] % seq[-1])
        if rmd == 0:
            break
        seq.append(rmd)

    def changes(signs):
        signs = [x for x in signs if x != 0]
        return sum(1 for a, b in zip(signs, signs[1:]) if a * b < 0)
    lc = lambda p: p.coeffs()[-1] if p.degree() >= 0 else 0
    at_pinf = [(1 if lc(p) > 0 else -1) for p in seq]
    at_minf = [(1 if lc(p) > 0 else -1) * (1 if p.degree() % 2 == 0 else -1) for p in seq]
    return changes(at_minf) - changes(at_pinf)


def _real_roots(q, prec=64):
    """Real roots of the squarefree q as rigorous acb->arb balls: the number of real roots is
    fixed exactly by a Sturm count; precision is raised until exactly that many isolating
    complex-root balls meet the real axis (non-real roots are then excluded rigorously)."""
    nreal = _sturm_count(q)
    zq = flint.fmpz_poly([int(x) for x in q.numer().coeffs()])
    old = flint.ctx.prec
    try:
        while prec <= 65536:
            flint.ctx.prec = prec
            rts = [r for r, m in zq.complex_roots() if r.imag.contains(0)]
            if len(rts) == nreal:
                return sorted([r.real for r in rts], key=lambda x: float(x.mid()))
            prec *= 2
    finally:
        flint.ctx.prec = old
    raise RuntimeError("real root isolation failed")


def _sign_at(poly_in_u, k, q):
    """Sign of poly(u) at the k-th real root of q.  poly is nonzero mod q (q irreducible), so it
    cannot vanish at a root; precision is doubled until the arb ball excludes 0 (rigorous)."""
    prec = 64
    while prec <= 65536:
        old = flint.ctx.prec
        try:
            x = _real_roots(q, prec)[k]
            flint.ctx.prec = prec
            v = flint.arb(0)
            for cf in reversed(poly_in_u.coeffs()):
                v = v * x + flint.arb(cf)
            if v > 0:
                return 1
            if v < 0:
                return -1
        finally:
            flint.ctx.prec = old
        prec *= 2
    return 0


def reconstruct_F(res):
    """Wang Sec. 6.1 for an engine1 F-result.  Returns list of dicts, one per retained
    irreducible factor q: {'q', 'status', 'points' (in K_q, frame p1=(0,0),p2=(1,0),p3=(0,1),
    metric G), 'embeddings' (real roots where G is positive definite)}."""
    U = f_univariate(res)
    if U.get("failed"):
        return {"failed": "M-choice (App. A.2) failed; case UNRESOLVED"}
    c, r = parse(res["pattern"])
    a = res["a"] - 1
    out = []
    for q, mult in U["factors"][1]:
        q = flint.fmpq_poly(q)
        rec = {"q": str(q)}
        roots = _real_roots(q)
        if not roots:
            rec["status"] = "no real root"
            out.append(rec)
            continue
        K = NF(q)
        z = {}
        for j in range(r):
            z[j] = flint.fmpq_poly([1]) if j == a else K.red(U["g"][j])
        zs = [z[j] for j in range(r)]
        if any(x == 0 for x in zs) or len({str(x) for x in zs}) < r:
            rec["status"] = "zero or equal label values in K_q"
            out.append(rec)
            continue
        D = lambda i, j: flint.fmpq_poly([0]) if i == j else zs[c[i][j]]
        g11, g22 = D(0, 1), D(0, 2)
        g12 = (D(0, 1) + D(0, 2) - D(1, 2)) / 2
        detG = K.red(g11 * g22 - g12 * g12)
        if detG == 0:
            rec["status"] = "det G = 0"
            out.append(rec)
            continue
        pos = [roots[k] for k in range(len(roots)) if _sign_at(K.red(g11), k, q) > 0 and _sign_at(detG, k, q) > 0]
        if not pos:
            rec["status"] = "G not positive definite at any real root"
            out.append(rec)
            continue
        pts = [(flint.fmpq_poly([0]), flint.fmpq_poly([0])), (flint.fmpq_poly([1]), flint.fmpq_poly([0])),
               (flint.fmpq_poly([0]), flint.fmpq_poly([1]))]
        # solve 2G [x,y]^T = rhs
        idet = K.inv(4 * detG)
        for i in range(3, 6):
            b1 = D(0, i) + D(0, 1) - D(1, i)
            b2 = D(0, i) + D(0, 2) - D(2, i)
            # [2g11, 2g12; 2g12, 2g22]^{-1} = 1/(4det) [2g22, -2g12; -2g12, 2g11]
            x = K.mul(idet, 2 * g22 * b1 - 2 * g12 * b2)
            y = K.mul(idet, -2 * g12 * b1 + 2 * g11 * b2)
            pts.append((x, y))
        nrm = lambda v: K.red(v[0] * v[0] * g11 + 2 * v[0] * v[1] * g12 + v[1] * v[1] * g22)
        ok = True
        for i, j in PAIRS:
            dv = (pts[i][0] - pts[j][0], pts[i][1] - pts[j][1])
            if K.red(nrm(dv) - D(i, j)) != 0:
                ok = False
        if not ok:
            rec["status"] = "distance verification failed"
            out.append(rec)
            continue
        gp = True
        for i, j, k in TRIPLES:
            d = K.red((pts[j][0] - pts[i][0]) * (pts[k][1] - pts[i][1]) - (pts[j][1] - pts[i][1]) * (pts[k][0] - pts[i][0]))
            if d == 0:
                gp = False
        for Q in QUADS:
            M = [[flint.fmpq_poly([1]), pts[i][0], pts[i][1], nrm(pts[i])] for i in Q]
            if K.red(_pdet(M, K)) == 0:
                gp = False
        rec["status"] = "REALIZED" if gp else "not in general position"
        rec["points"] = [[str(p[0]), str(p[1])] for p in pts]
        rec["G"] = [str(g11), str(g12), str(g22)]
        rec["embeddings"] = [str(x) for x in pos]
        out.append(rec)
    return out


def _pdet(M, K):
    n = len(M)
    if n == 1:
        return M[0][0]
    s = flint.fmpq_poly([0])
    for j in range(n):
        minor = [row[:j] + row[j + 1:] for row in M[1:]]
        t = K.mul(M[0][j], _pdet(minor, K))
        s = s + t if j % 2 == 0 else s - t
    return K.red(s)
