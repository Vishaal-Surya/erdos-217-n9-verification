"""Exact crescent-configuration checker (Erdos #217).

A set of n points in R^2 is a crescent configuration iff
  * the points are distinct,
  * no 3 are collinear, no 4 are concyclic,
  * the C(n,2) squared distances fall into exactly n-1 equality classes whose
    sizes are exactly {1, 2, ..., n-1}.

All arithmetic is EXACT.  Every coordinate is an element of one number field
K = Q[t]/(m(t)), m irreducible over Q (checked with FLINT factorisation).
Elements are stored as their unique reduced representative: a polynomial in t
of degree < deg m with rational (fmpq) coefficients.  Because m is irreducible,
K is a field and this representative is canonical, so

    x == y  in K   <=>   the reduced coefficient vectors are identical.

Only ring operations (+, -, *) are ever used (no division, no square roots,
no floating point), so every zero test / equality test is decided exactly.
Rationals are the case m(t) = t; Q(sqrt d) is m(t) = t^2 - d.

The zero tests are field-intrinsic, hence independent of which real embedding
of K is meant (a Galois conjugation maps 0 to 0 and preserves equalities).

Input JSON formats: see src/README_checker.md.  CLI:
    python src/checker.py file.json [--quiet]
"""
from __future__ import annotations

import itertools
import json
import sys
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Iterable, List, Sequence, Tuple

import flint
import sympy as sp

# ---------------------------------------------------------------------------
# Exact number field arithmetic
# ---------------------------------------------------------------------------


class NumberField:
    """K = Q[t]/(m(t)) with m monic-or-not, irreducible over Q."""

    def __init__(self, minpoly_coeffs_low_first: Sequence, name: str = "t",
                 label: str | None = None, check_irreducible: bool = True):
        m = flint.fmpq_poly([flint.fmpq(Fraction(c).numerator, Fraction(c).denominator)
                             for c in minpoly_coeffs_low_first])
        if m.degree() < 1:
            raise ValueError("minimal polynomial must have degree >= 1")
        if check_irreducible:
            _, facs = m.factor()
            if len(facs) != 1 or facs[0][1] != 1:
                raise ValueError(f"polynomial {m} is not irreducible over Q: {facs}")
        self.m = m
        self.deg = m.degree()
        self.name = name
        self.label = label or f"Q[{name}]/({m})"

    @classmethod
    def rationals(cls) -> "NumberField":
        return cls([0, 1], name="t", label="Q")

    @classmethod
    def sqrt(cls, d: int) -> "NumberField":
        d = int(d)
        if d == 1 or d == 0:
            raise ValueError("d must be a squarefree integer != 0, 1")
        if abs(d) != 1 and any(e > 1 for e in sp.factorint(abs(d)).values()):
            raise ValueError(f"d={d} is not squarefree")
        return cls([-d, 0, 1], name=f"sqrt({d})", label=f"Q(sqrt({d}))")

    def __call__(self, value) -> "NFElem":
        """Embed an int / Fraction / list of coefficients (low degree first)."""
        if isinstance(value, NFElem):
            assert value.K is self
            return value
        if isinstance(value, (list, tuple)):
            p = flint.fmpq_poly([_to_fmpq(c) for c in value])
        else:
            p = flint.fmpq_poly([_to_fmpq(value)])
        return NFElem(self, p % self.m)

    def gen(self) -> "NFElem":
        return NFElem(self, flint.fmpq_poly([0, 1]) % self.m)

    def __eq__(self, other):
        return isinstance(other, NumberField) and self.m == other.m

    def __hash__(self):
        return hash(str(self.m))


def _to_fmpq(c) -> flint.fmpq:
    if isinstance(c, flint.fmpq):
        return c
    if isinstance(c, flint.fmpz):
        return flint.fmpq(c, 1)
    if isinstance(c, sp.Rational):
        return flint.fmpq(int(c.p), int(c.q))
    f = Fraction(c) if not isinstance(c, str) else Fraction(c)
    return flint.fmpq(f.numerator, f.denominator)


class NFElem:
    __slots__ = ("K", "p")

    def __init__(self, K: NumberField, p: flint.fmpq_poly):
        self.K = K
        self.p = p  # always reduced mod K.m

    def _coerce(self, o):
        if isinstance(o, NFElem):
            if o.K is not self.K and o.K != self.K:
                raise ValueError("elements of different number fields")
            return o.p
        return flint.fmpq_poly([_to_fmpq(o)])

    def __add__(self, o):
        return NFElem(self.K, self.p + self._coerce(o))

    __radd__ = __add__

    def __sub__(self, o):
        return NFElem(self.K, self.p - self._coerce(o))

    def __rsub__(self, o):
        return NFElem(self.K, self._coerce(o) - self.p)

    def __neg__(self):
        return NFElem(self.K, -self.p)

    def __mul__(self, o):
        return NFElem(self.K, (self.p * self._coerce(o)) % self.K.m)

    __rmul__ = __mul__

    def is_zero(self) -> bool:
        return self.p.degree() < 0  # zero polynomial has degree -1

    def __eq__(self, o):
        return (self - o).is_zero()

    def key(self) -> Tuple[Tuple[int, int], ...]:
        """Canonical hashable key (reduced coefficients, low degree first)."""
        cs = self.p.coeffs()
        return tuple((int(c.p), int(c.q)) for c in cs)

    def __hash__(self):
        return hash(self.key())

    def __str__(self):
        if self.K.deg == 1:
            cs = self.p.coeffs()
            return str(Fraction(int(cs[0].p), int(cs[0].q))) if cs else "0"
        s = str(self.p)
        return s.replace("x", self.K.name) if s else "0"

    __repr__ = __str__

    def to_sympy(self, gen_expr=None):
        t = sp.Symbol(self.K.name) if gen_expr is None else gen_expr
        return sum(sp.Rational(int(c.p), int(c.q)) * t ** i for i, c in enumerate(self.p.coeffs()))


# ---------------------------------------------------------------------------
# Geometry predicates (ring operations only)
# ---------------------------------------------------------------------------

Point = Tuple[NFElem, NFElem]


def sqdist(P: Point, Q: Point) -> NFElem:
    dx = P[0] - Q[0]
    dy = P[1] - Q[1]
    return dx * dx + dy * dy


def orient(A: Point, B: Point, C: Point) -> NFElem:
    """2x2 determinant |B-A, C-A|; zero iff A,B,C collinear."""
    return (B[0] - A[0]) * (C[1] - A[1]) - (B[1] - A[1]) * (C[0] - A[0])


def incircle(A: Point, B: Point, C: Point, D: Point) -> NFElem:
    """The 4x4 lifted determinant |x y x^2+y^2 1| (rows A,B,C,D), computed by
    translating D to the origin (column operations, determinant unchanged up to
    the sign (-1)) and expanding the 3x3 determinant.  For A,B,C not collinear it
    is zero iff D lies on the circle through A,B,C."""
    rows = []
    for P in (A, B, C):
        dx = P[0] - D[0]
        dy = P[1] - D[1]
        rows.append((dx, dy, dx * dx + dy * dy))
    (a, b, c), (d, e, f), (g, h, i) = rows
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


# ---------------------------------------------------------------------------
# Checker
# ---------------------------------------------------------------------------


@dataclass
class CheckResult:
    ok: bool
    n: int
    reasons: List[str] = field(default_factory=list)
    classes: List[dict] = field(default_factory=list)  # [{sqdist, count, pairs}]
    multiplicities: List[int] = field(default_factory=list)

    def summary(self) -> str:
        lines = [f"n = {self.n}: {'CRESCENT (ok)' if self.ok else 'NOT crescent'}"]
        for r in self.reasons:
            lines.append("  - " + r)
        lines.append(f"  class sizes (sorted): {sorted(self.multiplicities)}")
        for c in self.classes:
            lines.append(f"    d^2 = {c['sqdist']:<30}  x{c['count']}  pairs {c['pairs']}")
        return "\n".join(lines)

    def to_json(self):
        return {"ok": self.ok, "n": self.n, "reasons": self.reasons,
                "multiplicities": self.multiplicities, "classes": self.classes}


def distance_classes(points: Sequence[Point]):
    """Exact grouping of the C(n,2) squared distances. Returns dict key->(value, pairs)."""
    groups = {}
    for i, j in itertools.combinations(range(len(points)), 2):
        d = sqdist(points[i], points[j])
        k = d.key()
        if k not in groups:
            groups[k] = [d, []]
        groups[k][1].append((i, j))
    return groups


def check(points: Sequence[Point], max_reasons: int = 50) -> CheckResult:
    n = len(points)
    res = CheckResult(ok=True, n=n)
    if n < 2:
        res.ok = False
        res.reasons.append("need at least 2 points")
        return res
    # 1. distinct
    keys = {}
    for i, P in enumerate(points):
        k = (P[0].key(), P[1].key())
        if k in keys:
            res.ok = False
            res.reasons.append(f"duplicate point: {i} == {keys[k]}")
        else:
            keys[k] = i
    # 2. no 3 collinear
    coll = set()
    for i, j, k in itertools.combinations(range(n), 3):
        if orient(points[i], points[j], points[k]).is_zero():
            coll.add((i, j, k))
            res.ok = False
            if len(res.reasons) < max_reasons:
                res.reasons.append(f"collinear triple: {(i, j, k)}")
    # 3. no 4 concyclic (quadruples containing a collinear triple are already
    #    flagged; for them the lifted determinant is not a concyclicity test)
    for q in itertools.combinations(range(n), 4):
        if any(t in coll for t in itertools.combinations(q, 3)):
            continue
        if incircle(*(points[i] for i in q)).is_zero():
            res.ok = False
            if len(res.reasons) < max_reasons:
                res.reasons.append(f"concyclic quadruple: {q}")
    # 4. distance multiplicities
    groups = distance_classes(points)
    mult = sorted(len(v[1]) for v in groups.values())
    res.multiplicities = mult
    res.classes = [{"sqdist": str(v[0]), "count": len(v[1]), "pairs": v[1]}
                   for v in sorted(groups.values(), key=lambda v: len(v[1]))]
    if mult != list(range(1, n)):
        res.ok = False
        res.reasons.append(f"distance class sizes {mult} != {list(range(1, n))}")
    return res


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------


def _poly_in(expr, sym, K: NumberField) -> NFElem:
    expr = sp.expand(expr)
    try:
        P = sp.Poly(expr, sym, domain=sp.QQ)
    except (sp.PolynomialError, sp.polys.polyerrors.GeneratorsNeeded, sp.CoercionFailed) as e:
        raise ValueError(f"cannot read {expr} as a polynomial in {sym} over Q: {e}")
    coeffs = list(reversed(P.all_coeffs()))
    return K(coeffs)


def parse_rational(s) -> NFElem:
    K = NumberField.rationals()
    return K(Fraction(str(s)))


def load_points(spec: dict):
    """Return (K, points) from a JSON-like dict (format in README_checker.md)."""
    fld = spec.get("field", "Q")
    raw = spec["points"]
    if fld == "Q" or fld == {"type": "Q"}:
        K = NumberField.rationals()
        pts = [(K(Fraction(str(x))), K(Fraction(str(y)))) for x, y in raw]
        return K, pts
    if isinstance(fld, dict) and fld.get("type") == "sqrt":
        d = int(fld["d"])
        K = NumberField.sqrt(d)
        s = sp.Symbol("s")
        loc = {"s": s}

        def conv(z):
            e = sp.sympify(str(z), locals=loc, rational=True)
            e = e.xreplace({sp.sqrt(d): s, sp.sqrt(sp.Integer(d)): s})
            e = sp.expand(e)
            # rewrite any remaining sqrt(k*d^?) like sqrt(12) -> 2*sqrt(3) was handled by sympify
            return _poly_in(e, s, K)

        return K, [(conv(x), conv(y)) for x, y in raw]
    if isinstance(fld, dict) and fld.get("type") == "minpoly":
        # field generated by a root 'a' of an irreducible polynomial
        var = fld.get("var", "a")
        a = sp.Symbol(var)
        m = sp.Poly(sp.sympify(fld["minpoly"], locals={var: a}), a, domain=sp.QQ)
        K = NumberField(list(reversed(m.all_coeffs())), name=var, label=f"Q({var}), {m.as_expr()}=0")

        def conv(z):
            return _poly_in(sp.sympify(str(z), locals={var: a}, rational=True), a, K)

        return K, [(conv(x), conv(y)) for x, y in raw]
    if isinstance(fld, dict) and fld.get("type") == "algebraic":
        # arbitrary closed-form algebraic numbers (radicals, CRootOf, ...);
        # a common primitive element is found with sympy's exact algorithms
        exprs = [[sp.sympify(str(x), rational=True), sp.sympify(str(y), rational=True)] for x, y in raw]
        flat = [e for P in exprs for e in P]
        return _algebraic_field(flat, len(raw))
    raise ValueError(f"unknown field spec {fld!r}")


def _algebraic_field(flat, npts):
    from sympy.polys.numberfields.subfield import primitive_element, to_number_field
    x = sp.Symbol("x")
    nonrat = [e for e in flat if not e.is_Rational]
    if not nonrat:
        K = NumberField.rationals()
        vals = [K(Fraction(int(e.p), int(e.q))) for e in flat]
    else:
        # unique generators, then a primitive element theta of Q(all of them)
        uniq = list(dict.fromkeys(nonrat))
        mp, cs = primitive_element(uniq, x)
        theta = sum(c * e for c, e in zip(cs, uniq))
        mpoly = sp.Poly(mp, x, domain=sp.QQ)
        K = NumberField(list(reversed(mpoly.all_coeffs())), name="theta",
                        label=f"Q(theta), theta = {theta}, minpoly {mp}")
        vals = []
        for e in flat:
            if e.is_Rational:
                vals.append(K(Fraction(int(e.p), int(e.q))))
                continue
            an = to_number_field(e, theta)  # exact: rep is a polynomial in theta
            # an.rep: coefficients high degree first, as a polynomial in an.root (=theta)
            coeffs = [sp.Rational(c) for c in an.rep.to_list()]
            # sanity: the generator sympy used must have our minimal polynomial
            if sp.Poly(an.minpoly.as_expr(), an.minpoly.gen, domain=sp.QQ).monic().all_coeffs() != mpoly.monic().all_coeffs():
                raise ValueError("internal: primitive element / minpoly mismatch")
            vals.append(K(list(reversed(coeffs))))
        K.theta_expr = theta
    pts = [(vals[2 * i], vals[2 * i + 1]) for i in range(npts)]
    return K, pts


def check_spec(spec: dict) -> CheckResult:
    _, pts = load_points(spec)
    return check(pts)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    quiet = "--quiet" in argv
    as_json = "--json" in argv
    files = [a for a in argv if not a.startswith("--")]
    if not files:
        print(__doc__)
        return 2
    rc = 0
    for fn in files:
        with open(fn, encoding="utf-8-sig") as fh:
            data = json.load(fh)
        specs = data if isinstance(data, list) else [data]
        for spec in specs:
            K, pts = load_points(spec)
            r = check(pts)
            name = spec.get("name", fn)
            if as_json:
                print(json.dumps({"name": name, "field": K.label, **r.to_json()}))
            else:
                print(f"[{name}] field {K.label}")
                print(r.summary() if not quiet else f"  ok={r.ok} sizes={r.multiplicities} {r.reasons[:3]}")
            if not r.ok:
                rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main())
