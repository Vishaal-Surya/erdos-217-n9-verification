"""Independent exact crescent-configuration checker (VERIFIER implementation).

Written from the definition in PROMPT.md section 1 only; does not read src/ or tests/.

Definition: n points in R^2 form a crescent configuration iff
  (a) no 3 collinear, (b) no 4 concyclic,
  (c) exactly n-1 distinct distances, and
  (d) the class sizes are exactly {1, 2, ..., n-1}.
Duplicate points are rejected (a zero distance is not a distance between distinct points).

Two exact arithmetic modes, sharing one generic algorithm:
  * Q3Field: elements of Q(sqrt3) as pairs (a, b) of Fractions meaning a + b*sqrt3,
    with hand-written +,-,*,/ and an exact zero test (a+b*sqrt3 = 0 iff a = b = 0,
    since sqrt3 is irrational) and exact sign test.
  * QField: plain Fractions.

Concyclicity is tested two independent ways (neither is the 4x4 lifted determinant):
  1. circumcentre O of three points (Cramer's rule, exact division), then
     |O - D|^2 == |O - A|^2;
  2. the cross-ratio test: A,B,C,D (as complex numbers) lie on a common circle or line
     iff Im[(A-C)(B-D) * conj((A-D)(B-C))] == 0  (cross ratio real), division-free.
  The two must agree whenever A,B,C are not collinear; disagreement raises.
"""
from fractions import Fraction as Fr
from itertools import combinations
from collections import Counter


# ---------------------------------------------------------------- fields
class QField:
    name = "Q"

    @staticmethod
    def c(x):
        return Fr(x)

    zero = Fr(0)
    add = staticmethod(lambda x, y: x + y)
    sub = staticmethod(lambda x, y: x - y)
    mul = staticmethod(lambda x, y: x * y)

    @staticmethod
    def div(x, y):
        if y == 0:
            raise ZeroDivisionError
        return x / y

    is_zero = staticmethod(lambda x: x == 0)

    @staticmethod
    def sign(x):
        return (x > 0) - (x < 0)

    key = staticmethod(lambda x: x)


class Q3Field:
    """a + b*sqrt(3), stored as tuple (Fraction a, Fraction b)."""
    name = "Q(sqrt3)"

    @staticmethod
    def c(a, b=0):
        return (Fr(a), Fr(b))

    zero = (Fr(0), Fr(0))

    @staticmethod
    def add(x, y):
        return (x[0] + y[0], x[1] + y[1])

    @staticmethod
    def sub(x, y):
        return (x[0] - y[0], x[1] - y[1])

    @staticmethod
    def mul(x, y):
        a, b = x
        c, d = y
        return (a * c + 3 * b * d, a * d + b * c)

    @staticmethod
    def div(x, y):
        c, d = y
        nrm = c * c - 3 * d * d  # nonzero unless c=d=0 (sqrt3 irrational)
        if nrm == 0:
            assert c == 0 and d == 0
            raise ZeroDivisionError
        conj = (c / nrm, -d / nrm)
        return Q3Field.mul(x, conj)

    @staticmethod
    def is_zero(x):
        return x[0] == 0 and x[1] == 0

    @staticmethod
    def sign(x):
        a, b = x
        sa, sb = (a > 0) - (a < 0), (b > 0) - (b < 0)
        if sa == 0:
            return sb
        if sb == 0 or sa == sb:
            return sa
        # opposite signs: compare a^2 with 3 b^2
        return sa if a * a > 3 * b * b else -sa

    key = staticmethod(lambda x: x)  # normalized pair is canonical -> hashable


# ---------------------------------------------------------------- geometry
def _cross(F, A, B, C):
    ux, uy = F.sub(B[0], A[0]), F.sub(B[1], A[1])
    vx, vy = F.sub(C[0], A[0]), F.sub(C[1], A[1])
    return F.sub(F.mul(ux, vy), F.mul(uy, vx))


def _d2(F, A, B):
    dx, dy = F.sub(A[0], B[0]), F.sub(A[1], B[1])
    return F.add(F.mul(dx, dx), F.mul(dy, dy))


def collinear(F, A, B, C):
    return F.is_zero(_cross(F, A, B, C))


def circumcentre(F, A, B, C):
    # 2(B-A).O = |B|^2-|A|^2 ; 2(C-A).O = |C|^2-|A|^2
    two = F.c(2)
    a11, a12 = F.mul(two, F.sub(B[0], A[0])), F.mul(two, F.sub(B[1], A[1]))
    a21, a22 = F.mul(two, F.sub(C[0], A[0])), F.mul(two, F.sub(C[1], A[1]))
    nA = F.add(F.mul(A[0], A[0]), F.mul(A[1], A[1]))
    r1 = F.sub(F.add(F.mul(B[0], B[0]), F.mul(B[1], B[1])), nA)
    r2 = F.sub(F.add(F.mul(C[0], C[0]), F.mul(C[1], C[1])), nA)
    det = F.sub(F.mul(a11, a22), F.mul(a12, a21))
    ox = F.div(F.sub(F.mul(r1, a22), F.mul(a12, r2)), det)
    oy = F.div(F.sub(F.mul(a11, r2), F.mul(r1, a21)), det)
    return (ox, oy)


def concyclic_circumcentre(F, A, B, C, D):
    O = circumcentre(F, A, B, C)
    return F.is_zero(F.sub(_d2(F, O, D), _d2(F, O, A)))


def _cmul(F, z, w):
    return (F.sub(F.mul(z[0], w[0]), F.mul(z[1], w[1])),
            F.add(F.mul(z[0], w[1]), F.mul(z[1], w[0])))


def _csub(F, z, w):
    return (F.sub(z[0], w[0]), F.sub(z[1], w[1]))


def concyclic_or_collinear_crossratio(F, A, B, C, D):
    """Cross ratio (A,B;C,D) real  <=>  the 4 points lie on a circle or a line."""
    num = _cmul(F, _csub(F, A, C), _csub(F, B, D))
    den = _cmul(F, _csub(F, A, D), _csub(F, B, C))
    # Im(num * conj(den)) = num_im*den_re - num_re*den_im
    im = F.sub(F.mul(num[1], den[0]), F.mul(num[0], den[1]))
    return F.is_zero(im)


# ---------------------------------------------------------------- checker
def check(points, F=Q3Field, want_reasons=True):
    """Return (ok: bool, info: dict). points: list of (x, y) in field F."""
    n = len(points)
    info = {"n": n, "field": F.name}
    if n < 2:
        info["reason"] = "n<2"
        return False, info
    keys = [(F.key(p[0]), F.key(p[1])) for p in points]
    if len(set(keys)) != n:
        info["reason"] = "duplicate point"
        return False, info
    for A, B, C in combinations(range(n), 3):
        if collinear(F, points[A], points[B], points[C]):
            info["reason"] = f"collinear {(A, B, C)}"
            return False, info
    for quad in combinations(range(n), 4):
        P = [points[i] for i in quad]
        t1 = concyclic_circumcentre(F, *P)
        t2 = concyclic_or_collinear_crossratio(F, *P)
        if t1 != t2:  # no 3 collinear here, so both tests must agree
            raise AssertionError(f"concyclic tests disagree on {quad}")
        if t1:
            info["reason"] = f"concyclic {quad}"
            return False, info
    cnt = Counter(F.key(_d2(F, points[i], points[j]))
                  for i, j in combinations(range(n), 2))
    sizes = sorted(cnt.values())
    info["class_sizes"] = sizes
    info["distinct"] = len(cnt)
    if len(cnt) != n - 1:
        info["reason"] = f"{len(cnt)} distinct distances, need {n-1}"
        return False, info
    if sizes != list(range(1, n)):
        info["reason"] = f"class sizes {sizes} != 1..{n-1}"
        return False, info
    info["reason"] = "crescent"
    return True, info


# ---------------------------------------------------------------- helpers
def q3pt(x, y):
    """x, y each given as (a, b) meaning a + b sqrt3, or a plain rational."""
    def conv(v):
        if isinstance(v, tuple):
            return Q3Field.c(*v)
        return Q3Field.c(v, 0)
    return (conv(x), conv(y))


def lattice_to_q3(i, j):
    """Triangular lattice point i*(1,0) + j*(1/2, sqrt3/2)."""
    return (Q3Field.c(Fr(i) + Fr(j, 2), 0), Q3Field.c(0, Fr(j, 2)))
