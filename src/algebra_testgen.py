"""Exact generators of REALIZABLE 6-point patterns (general position, <=6 distinct distances)
for the soundness tests of src/algebra.py.  All arithmetic exact: coordinates in Q(sqrt3)
represented as Q3(a, b) = a + b*sqrt(3) with Fraction a, b.
"""
from __future__ import annotations

import itertools
import random
from fractions import Fraction as Fr

PAIRS = [(i, j) for i in range(6) for j in range(i + 1, 6)]


class Q3:
    __slots__ = ("a", "b")

    def __init__(self, a, b=0):
        self.a = Fr(a)
        self.b = Fr(b)

    def __add__(self, o):
        o = o if isinstance(o, Q3) else Q3(o)
        return Q3(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __sub__(self, o):
        o = o if isinstance(o, Q3) else Q3(o)
        return Q3(self.a - o.a, self.b - o.b)

    def __rsub__(self, o):
        return Q3(o) - self

    def __neg__(self):
        return Q3(-self.a, -self.b)

    def __mul__(self, o):
        o = o if isinstance(o, Q3) else Q3(o)
        return Q3(self.a * o.a + 3 * self.b * o.b, self.a * o.b + self.b * o.a)

    __rmul__ = __mul__

    def __eq__(self, o):
        o = o if isinstance(o, Q3) else Q3(o)
        return self.a == o.a and self.b == o.b

    def __hash__(self):
        return hash((self.a, self.b))

    def iszero(self):
        return self.a == 0 and self.b == 0

    def __repr__(self):
        return f"({self.a}+{self.b}r3)"


def sqd(p, q):
    dx, dy = p[0] - q[0], p[1] - q[1]
    return dx * dx + dy * dy


def det(M):
    n = len(M)
    if n == 1:
        return M[0][0]
    s = Q3(0)
    for j in range(n):
        minor = [row[:j] + row[j + 1:] for row in M[1:]]
        t = M[0][j] * det(minor)
        s = s + t if j % 2 == 0 else s - t
    return s


def general_position(pts):
    for i, j, k in itertools.combinations(range(len(pts)), 3):
        a, b, c = pts[i], pts[j], pts[k]
        cr = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
        if cr.iszero():
            return False
    for Q in itertools.combinations(range(len(pts)), 4):
        rows = []
        for i in Q:
            p = pts[i]
            rows.append([Q3(1), p[0], p[1], p[0] * p[0] + p[1] * p[1]])
        if det(rows).iszero():
            return False
    return True


def pattern_of(pts):
    """Label string (first-occurrence renamed) or None if not 6 distinct points."""
    ds = [sqd(pts[i], pts[j]) for i, j in PAIRS]
    if any(d.iszero() for d in ds):
        return None
    m, out = {}, []
    for d in ds:
        if d not in m:
            m[d] = len(m) + 1
        out.append(str(m[d]))
    return "".join(out)


def realizable(pts):
    """Return pattern string if pts (6 points in Q(sqrt3)^2) is in general position with
    <= 6 distinct distances; else None."""
    s = pattern_of(pts)
    if s is None or max(int(ch) for ch in s) > 6:
        return None
    if not general_position(pts):
        return None
    return s


def tri_point(a, b):  # triangular lattice (2a+b, b*sqrt3)/2
    return (Q3(Fr(2 * a + b, 2)), Q3(0, Fr(b, 2)))


def sq_point(a, b):
    return (Q3(a), Q3(b))


PALASTI8 = [(Q3(0), Q3(1)), (Q3(0, 1), Q3(0)), (Q3(0, 2), Q3(0)),
            (Q3(0, Fr(5, 2)), Q3(Fr(5, 2))), (Q3(0, Fr(3, 2)), Q3(Fr(9, 2))),
            (Q3(0, Fr(1, 2)), Q3(Fr(7, 2))), (Q3(0, Fr(3, 2)), Q3(Fr(7, 2))), (Q3(0, 1), Q3(2))]


def rot120(p):  # rotation by 120 degrees: cos=-1/2, sin=sqrt3/2
    c, s = Q3(Fr(-1, 2)), Q3(0, Fr(1, 2))
    return (c * p[0] - s * p[1], s * p[0] + c * p[1])


def rot90(p):
    return (-p[1], p[0])


def gen_lattice(kind, radius, nsamples, rng):
    """Random 6-subsets of a lattice patch; returns list of (pattern, points)."""
    if kind == "tri":
        cells = [(a, b) for a in range(-radius, radius + 1) for b in range(-radius, radius + 1)
                 if abs(a + b) <= radius]
        mk = tri_point
    else:
        cells = [(a, b) for a in range(radius + 1) for b in range(radius + 1)]
        mk = sq_point
    out = []
    tries = 0
    while len(out) < nsamples and tries < nsamples * 400:
        tries += 1
        sub = rng.sample(cells, 6)
        pts = [mk(a, b) for a, b in sub]
        s = realizable(pts)
        if s:
            out.append((s, pts))
    return out


def gen_palasti():
    out = []
    for sub in itertools.combinations(range(8), 6):
        pts = [PALASTI8[i] for i in sub]
        s = realizable(pts)
        if s:
            out.append((s, pts))
    return out


def gen_two_triangles(nsamples, rng):
    """Two concentric equilateral triangles {p, Rp, R^2p, q, Rq, R^2q} (5 distances),
    random rational p, q.  Also variants with a common centre shifted."""
    out = []
    tries = 0
    while len(out) < nsamples and tries < nsamples * 50:
        tries += 1
        p = (Q3(Fr(rng.randint(-9, 9), rng.randint(1, 5))), Q3(Fr(rng.randint(-9, 9), rng.randint(1, 5))))
        q = (Q3(Fr(rng.randint(-9, 9), rng.randint(1, 5))), Q3(Fr(rng.randint(-9, 9), rng.randint(1, 5))))
        pts = [p, rot120(p), rot120(rot120(p)), q, rot120(q), rot120(rot120(q))]
        perm = list(range(6)); rng.shuffle(perm)
        pts = [pts[i] for i in perm]
        s = realizable(pts)
        if s:
            out.append((s, pts))
    return out


def gen_square_sym(nsamples, rng):
    """Orbits under rotation by 90 degrees: {p, Rp, R^2p, R^3p, q, R^2 q} etc. (rational)."""
    out = []
    tries = 0
    while len(out) < nsamples and tries < nsamples * 100:
        tries += 1
        p = (Q3(Fr(rng.randint(-9, 9), rng.randint(1, 4))), Q3(Fr(rng.randint(-9, 9), rng.randint(1, 4))))
        q = (Q3(Fr(rng.randint(-9, 9), rng.randint(1, 4))), Q3(Fr(rng.randint(-9, 9), rng.randint(1, 4))))
        orb_p = [p, rot90(p), rot90(rot90(p)), rot90(rot90(rot90(p)))]
        orb_q = [q, rot90(q), rot90(rot90(q)), rot90(rot90(rot90(q)))]
        choice = rng.choice([(3, 3), (2, 4), (4, 2)])
        pts = rng.sample(orb_p, choice[0]) + rng.sample(orb_q, choice[1])
        s = realizable(pts)
        if s:
            out.append((s, pts))
    return out


def gen_mixed_lattice(nsamples, rng):
    """Subsets of a 'rational rotation' square lattice: Z^2 united with Z^2 rotated by the
    Pythagorean angle (3/5, 4/5) -- still rational, many repeated distances."""
    base = [(a, b) for a in range(-3, 4) for b in range(-3, 4)]
    pts_all = set()
    for a, b in base:
        pts_all.add((Fr(a), Fr(b)))
        pts_all.add((Fr(3 * a - 4 * b, 5), Fr(4 * a + 3 * b, 5)))
    pts_all = sorted(pts_all)
    out = []
    tries = 0
    while len(out) < nsamples and tries < nsamples * 300:
        tries += 1
        sub = rng.sample(pts_all, 6)
        pts = [(Q3(x), Q3(y)) for x, y in sub]
        s = realizable(pts)
        if s:
            out.append((s, pts))
    return out


def canonical(s):
    """Canonical representative under S6 x label renaming (brute force, 720 perms)."""
    c = [[0] * 6 for _ in range(6)]
    for (i, j), ch in zip(PAIRS, s):
        c[i][j] = c[j][i] = ch
    best = None
    for perm in itertools.permutations(range(6)):
        lab = [c[perm[i]][perm[j]] for i, j in PAIRS]
        m, out = {}, []
        for x in lab:
            if x not in m:
                m[x] = str(len(m) + 1)
            out.append(m[x])
        t = "".join(out)
        if best is None or t < best:
            best = t
    return best


def all_realizable(seed=0, n_tri=400, n_sq=300, n_two=100, n_sym=150, n_mix=150):
    rng = random.Random(seed)
    groups = {
        "palasti": gen_palasti(),
        "tri": gen_lattice("tri", 3, n_tri, rng),
        "square": gen_lattice("sq", 4, n_sq, rng),
        "two_triangles": gen_two_triangles(n_two, rng),
        "square_sym": gen_square_sym(n_sym, rng),
        "mixed_rational": gen_mixed_lattice(n_mix, rng),
    }
    return groups


def build_cache(path, seed=0):
    """Generate a larger realizable set and store {group: [[pattern, [[x,y] as 'a+b*r3' strings]]]}."""
    import json
    rng = random.Random(seed)
    groups = {
        "palasti": gen_palasti(),
        "tri_r3": gen_lattice("tri", 3, 1200, rng),
        "tri_r4": gen_lattice("tri", 4, 800, rng),
        "square_5": gen_lattice("sq", 4, 500, rng),
        "square_6": gen_lattice("sq", 5, 300, rng),
        "two_triangles": gen_two_triangles(60, rng),
        "square_sym": gen_square_sym(80, rng),
        "mixed_rational": gen_mixed_lattice(20, rng),
    }
    enc = lambda q: [str(q.a), str(q.b)]
    data = {k: [[s, [[enc(p[0]), enc(p[1])] for p in pts]] for s, pts in v] for k, v in groups.items()}
    with open(path, "w") as f:
        json.dump(data, f)
    return data
