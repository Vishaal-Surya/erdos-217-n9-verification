"""Builder D: exact arithmetic in K_q = Q[u]/(q) and rigorous real-root sign evaluation.

Elements of K_q are flint.fmpq_poly reduced modulo the irreducible q (unique representative of
degree < deg q), so x == y in K_q iff the reduced polys are identical.  key(x) gives a hashable
canonical form.  Sign questions at a real root of q are decided by arb interval evaluation with
precision doubled until the ball excludes 0 (a nonzero element of the field cannot vanish at a
root of q, because q is irreducible); the number of real roots is fixed by an exact Sturm count.
"""
from __future__ import annotations

import flint
from flint import fmpq_poly, fmpq


class KField:
    def __init__(self, q):
        q = fmpq_poly(q)
        if q.degree() < 1:
            raise ValueError("q must have degree >= 1")
        # make q monic for a canonical description
        q = q / q.coeffs()[-1]
        self.q = q
        self.deg = q.degree()
        self.one = fmpq_poly([1])
        self.zero = fmpq_poly([])
        self._roots = None
        self._roots_prec = None

    # -- arithmetic
    def red(self, a):
        return a % self.q

    def mul(self, a, b):
        return (a * b) % self.q

    def inv(self, a):
        a = a % self.q
        if a == 0:
            raise ZeroDivisionError("inverse of 0 in K_q")
        if self.deg == 1 or a.degree() == 0:
            return fmpq_poly([1 / a.coeffs()[0]]) if a.degree() == 0 else self._xinv(a)
        return self._xinv(a)

    def _xinv(self, a):
        g, s, t = fmpq_poly.xgcd(a, self.q)
        if g.degree() != 0:
            raise ArithmeticError("q not irreducible (non-unit gcd)")
        return (s / g.coeffs()[0]) % self.q

    def div(self, a, b):
        return self.mul(a, self.inv(b))

    def const(self, x):
        return fmpq_poly([x]) if x != 0 else fmpq_poly([])

    @staticmethod
    def key(a):
        return tuple(a.coeffs())

    # -- real roots
    def sturm_count(self):
        q = self.q
        seq = [q, q.derivative()]
        while seq[-1].degree() > 0:
            rmd = -(seq[-2] % seq[-1])
            if rmd == 0:
                break
            seq.append(rmd)

        def changes(signs):
            s = [x for x in signs if x != 0]
            return sum(1 for a, b in zip(s, s[1:]) if a * b < 0)
        lc = lambda p: p.coeffs()[-1]
        at_p = [(1 if lc(p) > 0 else -1) for p in seq]
        at_m = [(1 if lc(p) > 0 else -1) * (1 if p.degree() % 2 == 0 else -1) for p in seq]
        return changes(at_m) - changes(at_p)

    def real_roots(self, prec=64):
        """Sorted list of arb balls, one per real root (isolating), at >= prec bits."""
        if self._roots is not None and self._roots_prec >= prec:
            return self._roots
        nreal = self.sturm_count()
        zq = flint.fmpz_poly([int(x) for x in self.q.numer().coeffs()])
        old = flint.ctx.prec
        try:
            p = max(prec, 64)
            while p <= 1 << 17:
                flint.ctx.prec = p
                rts = [r for r, m in zq.complex_roots() if r.imag.contains(0)]
                if len(rts) == nreal:
                    rr = sorted([r.real for r in rts], key=lambda x: float(x.mid()))
                    # isolation check: balls pairwise disjoint
                    ok = all(not rr[i].overlaps(rr[i + 1]) for i in range(len(rr) - 1))
                    if ok:
                        self._roots, self._roots_prec = rr, p
                        return rr
                p *= 2
        finally:
            flint.ctx.prec = old
        raise RuntimeError("real root isolation failed")

    def nreal(self):
        return self.sturm_count()

    def sign_at(self, a, k):
        """Exact sign of the field element a at the k-th real root (0 only if a == 0)."""
        a = a % self.q
        if a == 0:
            return 0
        if a.degree() == 0:
            c = a.coeffs()[0]
            return 1 if c > 0 else -1
        prec = 64
        while prec <= 1 << 17:
            rts = self.real_roots(prec)
            old = flint.ctx.prec
            try:
                flint.ctx.prec = self._roots_prec
                x = rts[k]
                v = flint.arb(0)
                for cf in reversed(a.coeffs()):
                    v = v * x + flint.arb(cf)
                if v > 0:
                    return 1
                if v < 0:
                    return -1
            finally:
                flint.ctx.prec = old
            prec = self._roots_prec * 2
        raise RuntimeError("sign evaluation did not separate from 0")

    def root_str(self, k, digits=30):
        rts = self.real_roots(128)
        old = flint.ctx.prec
        try:
            flint.ctx.prec = self._roots_prec
            return rts[k].str(digits, radius=True)
        finally:
            flint.ctx.prec = old

    def eval_float(self, a, k):
        rts = self.real_roots(128)
        old = flint.ctx.prec
        try:
            flint.ctx.prec = self._roots_prec
            v = flint.arb(0)
            for cf in reversed((a % self.q).coeffs()):
                v = v * rts[k] + flint.arb(cf)
            return float(v.mid())
        finally:
            flint.ctx.prec = old


def poly_from_list(lst):
    """['1/2','3',...] (ascending coefficients) -> fmpq_poly."""
    return fmpq_poly([fmpq(*map(int, s.split('/'))) if '/' in s else fmpq(int(s)) for s in lst])


def poly_to_list(p):
    return [str(c) for c in p.coeffs()]


# --------------------------------------------------------------------------- matrices over K


def kernel_over_K(K, M):
    """Right nullspace of M (list of lists of K elements).  Returns list of basis vectors in
    reduced form (identity on the free coordinates)."""
    m = len(M)
    n = len(M[0]) if m else 0
    A = [[K.red(x) for x in row] for row in M]
    piv_cols = []
    row = 0
    for col in range(n):
        pr = None
        for i in range(row, m):
            if A[i][col] != 0:
                pr = i
                break
        if pr is None:
            continue
        A[row], A[pr] = A[pr], A[row]
        iv = K.inv(A[row][col])
        A[row] = [K.mul(iv, x) for x in A[row]]
        for i in range(m):
            if i != row and A[i][col] != 0:
                f = A[i][col]
                A[i] = [K.red(x - f * y) for x, y in zip(A[i], A[row])]
        piv_cols.append(col)
        row += 1
        if row == m:
            break
    free = [c for c in range(n) if c not in piv_cols]
    basis = []
    for fc in free:
        v = [K.zero] * n
        v[fc] = K.one
        for i, pc in enumerate(piv_cols):
            v[pc] = K.red(-A[i][fc])
        basis.append(v)
    return basis, free
