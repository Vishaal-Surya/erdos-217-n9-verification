"""V-a independent implementation (from Wang arXiv:2611.03688 TeX, Lemmas 2.3, 4.1, 4.2, Procedure 4.4).

Written without reading src/ or tests/. Exact arithmetic over Z/Q via python-flint fmpz_mat
(rank / rref / nullspace are exact over Q for integer matrices).

Conventions:
  pattern string s: labels on pairs (0,1),(0,2),...,(4,5) (0-based points), digits 1..r.
  label variable z_l, l = 1..r  <->  exponent index l-1.
"""
import itertools
from functools import lru_cache
import flint

N = 6
PAIRS = [(i, j) for i in range(N) for j in range(i + 1, N)]
PIDX = {}
for _k, (_i, _j) in enumerate(PAIRS):
    PIDX[(_i, _j)] = _k
    PIDX[(_j, _i)] = _k


def lab(s, i, j):
    return int(s[PIDX[(i, j)]])


# ---------------------------------------------------------------- polynomials in pair variables
# A polynomial in the 15 pair variables y_p is a dict {sorted tuple of pair indices: int coeff}.
def pmul(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = tuple(sorted(ka + kb))
            out[k] = out.get(k, 0) + va * vb
    return {k: v for k, v in out.items() if v}


def padd(a, b, sb=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + sb * v
    return {k: v for k, v in out.items() if v}


def Dvar(i, j):
    return {} if i == j else {(PIDX[(i, j)],): 1}


def perm_sign(p):
    s = 1
    p = list(p)
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


def det_poly(M):
    n = len(M)
    out = {}
    for sig in itertools.permutations(range(n)):
        term = {(): perm_sign(sig)}
        for i in range(n):
            term = pmul(term, M[i][sig[i]])
            if not term:
                break
        out = padd(out, term)
    return out


# Gram-type matrix L_ij = D_0i + D_0j - D_ij, i,j = 1..5 (0-based points)
def _Lentry(i, j):
    return padd(padd(Dvar(0, i), Dvar(0, j)), Dvar(i, j), -1)


LMAT = [[_Lentry(i, j) for j in range(1, N)] for i in range(1, N)]
MINORS = []  # all 100 3x3 minors (duplicates harmless)
for R in itertools.combinations(range(5), 3):
    for C in itertools.combinations(range(5), 3):
        MINORS.append(det_poly([[LMAT[a][b] for b in C] for a in R]))

DET4 = {}
for Q in itertools.combinations(range(N), 4):
    DET4[Q] = det_poly([[Dvar(a, b) for b in Q] for a in Q])

TRI = {}
for (i, j, k) in itertools.combinations(range(N), 3):
    Dij, Dik, Djk = Dvar(i, j), Dvar(i, k), Dvar(j, k)
    s = padd(padd(Dij, Dik), Djk, -1)
    TRI[(i, j, k)] = padd(pmul({(): 4}, pmul(Dij, Dik)), pmul(s, s), -1)


# ---------------------------------------------------------------- monomials in z
@lru_cache(None)
def monos(r, d):
    out = []
    def rec(prefix, left, pos):
        if pos == r - 1:
            out.append(tuple(prefix + [left]))
            return
        for e in range(left, -1, -1):
            rec(prefix + [e], left - e, pos + 1)
    rec([], d, 0)
    return out


@lru_cache(None)
def mono_index(r, d):
    return {m: i for i, m in enumerate(monos(r, d))}


def subst(poly, s, r):
    """pair-variable polynomial -> {z exponent tuple: coeff}."""
    out = {}
    for k, v in poly.items():
        e = [0] * r
        for p in k:
            e[int(s[p]) - 1] += 1
        e = tuple(e)
        out[e] = out.get(e, 0) + v
    return {k: v for k, v in out.items() if v}


def zmul(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = tuple(x + y for x, y in zip(ka, kb))
            out[k] = out.get(k, 0) + va * vb
    return {k: v for k, v in out.items() if v}


def vec(poly, r, d):
    idx = mono_index(r, d)
    v = [0] * len(idx)
    for k, c in poly.items():
        v[idx[k]] += c
    return v


# ---------------------------------------------------------------- linear algebra helpers
def mat(rows, ncols):
    if not rows:
        return flint.fmpz_mat(0, ncols)
    return flint.fmpz_mat(rows)


def basis_rows(rows, ncols):
    """independent integer rows spanning the same Q-space (rref, primitive)."""
    if not rows:
        return []
    M, den, rk = flint.fmpz_mat(rows).rref()
    out = []
    for i in range(rk):
        row = [int(M[i, j]) for j in range(ncols)]
        out.append(row)
    return out


def rank(rows, ncols):
    if not rows:
        return 0
    return flint.fmpz_mat(rows).rank()


# ---------------------------------------------------------------- affine tests (Lemma 2.3)
def affine(s):
    """returns (passed, reason, list of linear forms (rows of A) as int lists over z_1..z_r)."""
    r = max(int(ch) for ch in s)
    rows = set()
    for a, b, c, d in itertools.permutations(range(N), 4):
        if lab(s, a, b) == lab(s, c, d) and lab(s, a, c) == lab(s, b, d):
            v = [0] * N
            v[a] += 1; v[d] += 1; v[b] -= 1; v[c] -= 1
            rows.add(tuple(v))
    for a, b, c in itertools.combinations(range(N), 3):
        if lab(s, a, b) == lab(s, a, c) == lab(s, b, c):
            for d in range(N):
                if d in (a, b, c):
                    continue
                if lab(s, a, d) == lab(s, b, d) == lab(s, c, d):
                    v = [0] * N
                    v[a] += 1; v[b] += 1; v[c] += 1; v[d] -= 3
                    rows.add(tuple(v))
    rows = [list(v) for v in rows]
    rk = rank(rows, N)
    if rk:
        for I in itertools.combinations(range(N), 3):
            sub = [[v[j] for j in range(N) if j not in I] for v in rows]
            if rank(sub, N - 3) < rk:
                return False, "support%s" % (I,), None
    # linear distance equations
    Rb = basis_rows(rows, N) if rows else []
    A = []
    for rho in Rb:
        for j in range(1, N):
            form = [0] * r
            for i in range(N):
                if rho[i] == 0:
                    continue
                if i != j:
                    form[lab(s, i, j) - 1] += rho[i]
                if i != 0:
                    form[lab(s, i, 0) - 1] -= rho[i]
            if any(form):
                A.append(form)
    A = basis_rows(A, r) if A else []
    ra = len(A)
    if ra:
        for l in range(r):
            e = [0] * r; e[l] = 1
            if rank(A + [e], r) == ra:
                return False, "zero z%d" % (l + 1), None
        for l in range(r):
            for m in range(l + 1, r):
                e = [0] * r; e[l] = 1; e[m] = -1
                if rank(A + [e], r) == ra:
                    return False, "equal z%d z%d" % (l + 1, m + 1), None
    return True, "ok", A


# ---------------------------------------------------------------- quartic space (section 4)
def W0(s, A):
    r = max(int(ch) for ch in s)
    gens = []
    seen = set()
    for mn in MINORS:
        c = subst(mn, s, r)
        if not c:
            continue
        for l in range(r):
            e = [0] * r; e[l] = 1
            q = zmul(c, {tuple(e): 1})
            v = tuple(vec(q, r, 4))
            if v not in seen:
                seen.add(v); gens.append(list(v))
    for form in A:
        lf = {}
        for l in range(r):
            if form[l]:
                e = [0] * r; e[l] = 1
                lf[tuple(e)] = form[l]
        for m in monos(r, 3):
            v = tuple(vec(zmul(lf, {m: 1}), r, 4))
            if v not in seen:
                seen.add(v); gens.append(list(v))
    return basis_rows(gens, len(monos(r, 4)))


def factors(s):
    r = max(int(ch) for ch in s)
    out = []
    for l in range(r):
        e = [0] * r; e[l] = 1
        out.append(("z%d" % (l + 1), {tuple(e): 1}, 1))
    for l in range(r):
        for t in range(l + 1, r):
            e1 = [0] * r; e1[l] = 1
            e2 = [0] * r; e2[t] = 1
            out.append(("z%d-z%d" % (l + 1, t + 1), {tuple(e1): 1, tuple(e2): -1}, 1))
    for (i, j, k), poly in TRI.items():
        out.append(("tri%d%d%d" % (i + 1, j + 1, k + 1), subst(poly, s, r), 2))
    return out


def divide(W, f, e, r):
    """one application of Lemma 4.1: returns new basis (list of rows)."""
    n4 = len(monos(r, 4))
    if len(W) == n4:
        return W
    Q = monos(r, 4 - e)
    FM = [vec(zmul(f, {m: 1}), r, 4) for m in Q]
    S = FM + W
    # left kernel of S  == nullspace of S^T
    ST = flint.fmpz_mat(S).transpose()
    X, nul = ST.nullspace()
    if nul == 0:
        return W
    nq = len(Q)
    new = []
    me = monos(r, e)
    for c in range(nul):
        coeffs = [int(X[i, c]) for i in range(nq)]
        if not any(coeffs):
            continue
        q = {Q[i]: coeffs[i] for i in range(nq) if coeffs[i]}
        for m in me:
            new.append(vec(zmul(q, {m: 1}), r, 4))
    if not new:
        return W
    return basis_rows(W + new, n4)


def in_span(W, v, n4):
    if not any(v):
        return True
    return rank(W + [v], n4) == len(W)


def targets(s, r):
    """yields (kind, key, vector) exclusion targets."""
    for m in monos(r, 4):
        v = [0] * len(monos(r, 4)); v[mono_index(r, 4)[m]] = 1
        yield ("monomial", m, v)
    for Q, poly in DET4.items():
        yield ("det4", Q, vec(subst(poly, s, r), r, 4))


def has_monomial(W, n4):
    """W in rref: contains a monomial iff some row has exactly one nonzero entry."""
    for row in W:
        if sum(1 for x in row if x) == 1:
            return True
    return False


def excluded(W, s, r):
    n4 = len(monos(r, 4))
    if has_monomial(W, n4):
        return "monomial"
    for Q, poly in DET4.items():
        if in_span(W, vec(subst(poly, s, r), r, 4), n4):
            return "det4%s" % (Q,)
    return None


def closure(s, W):
    r = max(int(ch) for ch in s)
    n4 = len(monos(r, 4))
    facs = factors(s)
    while True:
        changed = False
        for name, f, e in facs:
            if len(W) == n4:
                return W
            W2 = divide(W, f, e, r)
            if len(W2) != len(W):
                changed = True
                W = W2
        if not changed:
            return W


def finite_labels(W, r):
    """labels a (1-based) with V = W + z_a Q[z]_3."""
    ms = monos(r, 4)
    out = []
    for a in range(r):
        cols = [i for i, m in enumerate(ms) if m[a] == 0]
        sub = [[row[c] for c in cols] for row in W]
        if rank(sub, len(cols)) == len(cols):
            out.append(a + 1)
    return out


def classify(s):
    """full independent pipeline: affine -> W0 -> closure -> exclusion -> F/U."""
    ok, why, A = affine(s)
    r = max(int(ch) for ch in s)
    if not ok:
        return ("AFFINE", why, None)
    W = W0(s, A)
    W = closure(s, W)
    ex = excluded(W, s, r)
    if ex:
        return ("EXCL", ex, len(W))
    fl = finite_labels(W, r)
    h = len(monos(r, 4)) - len(W)
    if r in (5, 6) and fl:
        return ("F", fl[0], h)
    return ("U", r, h)


def replay(s, cert):
    """replay a builder certificate over Q. returns (ok, dims_match, msg)."""
    r = max(int(ch) for ch in s)
    n4 = len(monos(r, 4))
    ok, why, A = affine(s)
    if not ok:
        return False, False, "pattern fails affine (%s)" % why
    W = W0(s, A)
    facs = {name: (f, e) for name, f, e in factors(s)}
    dims_ok = True
    for name, before, after in cert["trace"]:
        if name not in facs:
            return False, False, "unknown factor %s" % name
        if len(W) != before:
            dims_ok = False
        f, e = facs[name]
        W = divide(W, f, e, r)
        if len(W) != after:
            dims_ok = False
    kind, key = cert["target"]
    if kind == "monomial":
        key = tuple(key)
        if len(key) != r or sum(key) != 4:
            return False, dims_ok, "bad monomial %s" % (key,)
        v = [0] * n4; v[mono_index(r, 4)[key]] = 1
    elif kind == "det4":
        Q = tuple(key)
        if Q not in DET4:
            return False, dims_ok, "bad det4 %s" % (Q,)
        v = vec(subst(DET4[Q], s, r), r, 4)
        if not any(v):
            # det identically zero would be a vacuous target; flag it
            return False, dims_ok, "det4 identically zero"
    else:
        return False, dims_ok, "unknown target"
    if in_span(W, v, n4):
        return True, dims_ok, "ok"
    return False, dims_ok, "target not in W"
