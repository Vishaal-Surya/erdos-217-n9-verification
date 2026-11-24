"""Builder D: exact extension of a reconstructed core S (6 points in K_q^2, metric G) to an
n-point crescent configuration (n = 9 is the target; n = 7, 8 are used for known-answer tests).

Setting (see NOTES/extension.md for the proofs).  P is an n-point crescent configuration with
multiplicities 1..n-1; low values delta_1, delta_2 (3 low pairs), high values delta_3..delta_{n-1}
(n-3 values).  S subset P is a 6-point set all of whose pairs are high, r = |D(S)|, and
m = (n-3) - r is the number of high values absent from S.  We require m in {0, 1}.

  Lemma T (trilateration): x is determined by squared distances to a noncollinear triangle.
  Lemma H (high-linked): among the n-6 outside points at most one has >= 2 low pairs to S, so at
     least n-7 are high-linked (<= 1 low pair to S, hence >= 5 high pairs to S).  Each high-linked
     point is a trilateration point of a triangle of S with values in D(S)^3 or (m = 1 only) a
     circumcentre of a triangle of S.  m = 0: triangles (0,1,2),(3,4,5) suffice.
  Lemma L (last point): for |Q| = n-1, the last point is a circumcentre of a triangle of Q or a
     trilateration point of a triangle of Q with values in D(Q)^3; for the latter, a family of
     triangles covering every (n-4)-subset of Q suffices.
  Capacity tests (necessary): injective mu: D(Q) -> {1..n-1}, m_Q(d) <= mu(d), mu >= 3 on D(S),
     each selected high-linked point has >= 5 pairs to S with mu >= 3.

All tests are exact in K_q (zero tests of reduced polynomials), so one search covers every real
embedding of K_q (in particular every PD one).
"""
from __future__ import annotations

import itertools
from collections import Counter

from extend_kfield import KField


def fits(counts, caps):
    """Lemma 5.2: sorted counts a fit sorted capacities b (injectively, a_i <= assigned cap)."""
    a = sorted(counts)
    b = sorted(caps)
    s, t = len(a), len(b)
    if s > t:
        return False
    return all(a[i] <= b[t - s + i] for i in range(s))


def covering_triples(npts, s):
    """Greedy family of triples of range(npts) such that every s-subset contains one; verified."""
    allT = list(itertools.combinations(range(npts), 3))
    subsets = [frozenset(x) for x in itertools.combinations(range(npts), s)]
    fam = []
    uncovered = set(subsets)
    while uncovered:
        best = max(allT, key=lambda T: sum(1 for U in uncovered if set(T) <= U))
        fam.append(best)
        uncovered = {U for U in uncovered if not set(best) <= U}
    for U in subsets:  # exhaustive verification
        assert any(set(T) <= U for T in fam)
    return fam


class Geo:
    """Points (x, y) in K^2 with the metric G = [[g11, g12], [g12, g22]]."""

    def __init__(self, K: KField, G):
        self.K = K
        self.g11, self.g12, self.g22 = G

    def q(self, x, y):
        K = self.K
        return K.red(x * x * self.g11 + 2 * x * y * self.g12 + y * y * self.g22)

    def b(self, v, w):
        K = self.K
        return K.red(v[0] * w[0] * self.g11 + (v[0] * w[1] + v[1] * w[0]) * self.g12 + v[1] * w[1] * self.g22)

    def sqd(self, p, r):
        return self.q(p[0] - r[0], p[1] - r[1])

    def orient(self, p, r, s):
        return self.K.red((r[0] - p[0]) * (s[1] - p[1]) - (r[1] - p[1]) * (s[0] - p[0]))

    def incircle(self, p, r, s, t):
        """det of rows (1, x, y, |.|_G^2) for 4 points (nonzero iff not concyclic, given no
        3 collinear and G positive definite)."""
        K = self.K
        n0 = self.q(*p)
        rows = []
        for u in (r, s, t):
            rows.append((u[0] - p[0], u[1] - p[1], self.q(*u) - n0))
        (a, b_, c), (d, e, f), (g, h, i) = rows
        return K.red(a * K.red(e * i - f * h) - b_ * K.red(d * i - f * g) + c * K.red(d * h - e * g))

    def key(self, p):
        return (KField.key(p[0]), KField.key(p[1]))

    # -- Lemma T
    def tri_setup(self, a, b, c):
        """Precompute trilateration data for triangle (a, b, c); None if collinear."""
        K = self.K
        r1 = (b[0] - a[0], b[1] - a[1])
        r2 = (c[0] - a[0], c[1] - a[1])
        # rows 2 r_i^T G
        m11 = K.red(2 * (r1[0] * self.g11 + r1[1] * self.g12))
        m12 = K.red(2 * (r1[0] * self.g12 + r1[1] * self.g22))
        m21 = K.red(2 * (r2[0] * self.g11 + r2[1] * self.g12))
        m22 = K.red(2 * (r2[0] * self.g12 + r2[1] * self.g22))
        det = K.red(m11 * m22 - m12 * m21)
        if det == 0:
            return None
        idet = K.inv(det)
        i11, i12 = K.mul(m22, idet), K.mul(-m12, idet)
        i21, i22 = K.mul(-m21, idet), K.mul(m11, idet)
        na = self.q(*a)
        s1 = K.red(self.q(*b) - na)
        s2 = K.red(self.q(*c) - na)
        x0 = (K.red(i11 * s1 + i12 * s2), K.red(i21 * s1 + i22 * s2))  # circumcentre
        w1 = (i11, i21)
        w2 = (i12, i22)
        y0 = (x0[0] - a[0], x0[1] - a[1])
        return {"a": a, "x0": x0, "w1": w1, "w2": w2,
                "Q0": self.q(*y0), "B1": K.red(2 * self.b(y0, w1)), "B2": K.red(2 * self.b(y0, w2)),
                "Q1": self.q(*w1), "B12": K.red(2 * self.b(w1, w2)), "Q2": self.q(*w2)}

    def trilaterate(self, T, vals):
        """All x with |x-a|^2=al, |x-b|^2=be, |x-c|^2=ga for (al, be, ga) in vals^3, plus the
        circumcentre.  Returns list of points (exact, with the first distance verified)."""
        K = self.K
        out = [T["x0"]]
        Q0, B1, B2, Q1, B12, Q2 = T["Q0"], T["B1"], T["B2"], T["Q1"], T["B12"], T["Q2"]
        vq = {}
        for al in vals:
            for ga in vals:
                v = al - ga
                vq[(KField.key(al), KField.key(ga))] = (v, K.mul(v, Q2))
        for al in vals:
            ka = KField.key(al)
            for be in vals:
                u = al - be
                Au = K.red(Q0 + K.mul(u, B1) + K.mul(K.mul(u, u), Q1) - al)
                Bu = K.red(B2 + K.mul(u, B12))
                for ga in vals:
                    v, vQ2 = vq[(ka, KField.key(ga))]
                    if K.red(Au + v * (Bu + vQ2)) == 0:
                        x = (K.red(T["x0"][0] + u * T["w1"][0] + v * T["w2"][0]),
                             K.red(T["x0"][1] + u * T["w1"][1] + v * T["w2"][1]))
                        out.append(x)
        return out


class Extender:
    def __init__(self, K: KField, G, S, n=9, cover_cache={}, full=False):
        self.full = full  # cross-check mode: all triangles everywhere (no covering shortcuts)
        self.K, self.geo, self.S, self.n = K, Geo(K, G), list(S), n
        g = self.geo
        self.Skeys = {g.key(p) for p in S}
        vals = {}
        for i, j in itertools.combinations(range(6), 2):
            d = g.sqd(S[i], S[j])
            vals.setdefault(KField.key(d), d)
        self.DS = vals  # key -> value
        self.r = len(vals)
        self.m = (n - 3) - self.r
        if self.m not in (0, 1):
            raise ValueError(f"need |D(S)| in {{n-3, n-4}}; got r={self.r}, n={n}")
        self.highs = list(range(3, n))
        self.need = max(0, (n - 6) - 1)
        key = (n - 1, n - 4)
        if key not in cover_cache:
            cover_cache[key] = covering_triples(n - 1, n - 4)
        self.cover = cover_cache[key]
        self.stats = Counter()
        self.found = []
        self.near = []

    # -- general position, incremental
    def gp_add(self, Q, x):
        g = self.geo
        for i, j in itertools.combinations(range(len(Q)), 2):
            if g.orient(Q[i], Q[j], x) == 0:
                return False
        for i, j, k in itertools.combinations(range(len(Q)), 3):
            if g.incircle(Q[i], Q[j], Q[k], x) == 0:
                return False
        return True

    # -- capacity test on Q = S + outside (outside = selected high-linked points)
    def capacity(self, Q, dist_keys, n_out_sel):
        """dist_keys[(i,j)] = key of squared distance.  The first 6 points of Q are S, the next
        n_out_sel points are the selected high-linked points."""
        cnt = Counter(dist_keys.values())
        DS = self.DS
        core = [cnt[d] for d in DS]
        nonc = {d: c for d, c in cnt.items() if d not in DS}
        lows = [1, 2]
        xs = list(range(6, 6 + n_out_sel))
        xv = {x: [dist_keys[(i, x)] for i in range(6)] for x in xs}
        if self.m == 0:
            return (fits(core, self.highs) and fits(list(nonc.values()), lows)
                    and all(sum(v in DS for v in xv[x]) >= 5 for x in xs))
        for k in self.highs:
            if not fits(core, [h for h in self.highs if h != k]):
                continue
            if fits(list(nonc.values()), lows) and all(sum(v in DS for v in xv[x]) >= 5 for x in xs):
                return True
            for d0, c0 in nonc.items():
                if c0 > k:
                    continue
                rest = [c for d, c in nonc.items() if d != d0]
                if fits(rest, lows) and all(sum((v in DS or v == d0) for v in xv[x]) >= 5 for x in xs):
                    return True
        return False

    def high_linked_candidates(self):
        g = self.geo
        S = self.S
        vals = list(self.DS.values())
        tris = [(0, 1, 2), (3, 4, 5)] if (self.m == 0 and not self.full) else list(itertools.combinations(range(6), 3))
        cands = {}
        for t in tris:
            T = g.tri_setup(S[t[0]], S[t[1]], S[t[2]])
            assert T is not None  # S is in general position
            for x in g.trilaterate(T, vals):
                self.stats["hl_trials_hit"] += 1
                k = g.key(x)
                if k not in self.Skeys:
                    cands.setdefault(k, x)
        self.stats["hl_trials"] += len(tris) * (len(vals) ** 3 + 1)
        return list(cands.values())

    def dists(self, Q):
        g = self.geo
        return {(i, j): KField.key(g.sqd(Q[i], Q[j])) for i, j in itertools.combinations(range(len(Q)), 2)}

    def run(self):
        g, S, n = self.geo, self.S, self.n
        dS = self.dists(S)
        # stage 1: high-linked candidates, 7-point tests
        if self.need > 0:
            cands = self.high_linked_candidates()
            self.stats["hl_candidates"] = len(cands)
            good = []
            for x in cands:
                if not self.gp_add(S, x):
                    continue
                Q = S + [x]
                d = dict(dS)
                for i in range(6):
                    d[(i, 6)] = KField.key(g.sqd(S[i], x))
                if self.capacity(Q, d, 1):
                    good.append(x)
            self.stats["hl_pass7"] = len(good)
        else:
            good = []
        # stage 2: sets of `need` high-linked points
        Qs = []
        if self.need == 0:
            Qs.append(S)
        else:
            for combo in itertools.combinations(range(len(good)), self.need):
                Q = list(S)
                ok = True
                for idx in combo:
                    x = good[idx]
                    if not self.gp_add(Q, x):
                        ok = False
                        break
                    Q.append(x)
                if not ok:
                    continue
                if self.capacity(Q, self.dists(Q), self.need):
                    Qs.append(Q)
        self.stats["Q_pass"] += len(Qs)
        # stage 3: last point
        for Q in Qs:
            self.last_point(Q)
        return self.found

    def last_point(self, Q):
        g, n = self.geo, self.n
        assert len(Q) == n - 1
        dQ = self.dists(Q)
        vals = {}
        for (i, j), kk in dQ.items():
            if kk not in vals:
                vals[kk] = g.sqd(Q[i], Q[j])
        vlist = list(vals.values())
        Qkeys = {g.key(p) for p in Q}
        cands = {}
        cover = set(self.cover) if not self.full else set(itertools.combinations(range(n - 1), 3))
        for t in itertools.combinations(range(n - 1), 3):
            T = g.tri_setup(Q[t[0]], Q[t[1]], Q[t[2]])
            if T is None:
                raise AssertionError("collinear triple in a GP set")
            pts = g.trilaterate(T, vlist) if t in cover else [T["x0"]]
            self.stats["last_trials"] += (len(vlist) ** 3 + 1) if t in cover else 1
            for x in pts:
                k = g.key(x)
                if k not in Qkeys:
                    cands.setdefault(k, x)
        self.stats["last_candidates"] += len(cands)
        cntQ = Counter(dQ.values())
        target = list(range(1, n))
        for x in cands.values():
            cnt = Counter(cntQ)
            for i in range(n - 1):
                cnt[KField.key(g.sqd(Q[i], x))] += 1
            if sorted(cnt.values()) != target:
                continue
            self.stats["last_mult_ok"] += 1
            if not self.gp_add(Q, x):
                self.near.append(Q + [x])  # crescent multiplicities but NOT in general position
                continue
            P = Q + [x]
            core_ok = all(cnt[d] >= 3 for d in self.DS)
            self.found.append({"points": P, "core_is_high": core_ok})
            self.stats["FOUND"] += 1
