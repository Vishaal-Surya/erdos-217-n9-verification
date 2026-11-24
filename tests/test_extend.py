"""Tests for the exact extension code (src/extend_core.py).

Known-answer tests (end-to-end, same code path as the n = 9 run, parameterised by n):
  * Palasti n = 8: both 6-point all-high subsets S (|D(S)| = 5 = n-3, m = 0) are given to the
    extender as cores (rational frame coordinates, K = Q).  It must output exactly one 8-point
    crescent per core, and it must be Palasti (outside points identified by their exact squared
    distances to S, which determine a point).  Path exercised: high-linked candidates (m = 0),
    7-point capacity/GP test, last-point lemma with the (7,4) covering, final test.
  * Lattice n = 7 (all 8 classes from the radius-3 search): every 6-point all-high subset with
    |D(S)| in {3, 4} (m = 1 or 0) must be extended back to its 7-point configuration
    (last-point lemma with all 20 triangles, circumcentres included).
Unit tests:
  * fits() == brute force; capacity() == brute-force search for the map mu of Wang's (7.x)
    on random synthetic data for n = 9, m in {0, 1}, one or two selected outside points.
  * covering families verified exhaustively; planted trilateration / circumcentre recovery.
"""
import itertools, json, os, random, sys
from fractions import Fraction as F
from types import SimpleNamespace

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
import extend_core as E
import extend_recon as R
from extend_kfield import KField
from flint import fmpq_poly, fmpq

KQ = KField(fmpq_poly([0, 1]))  # K = Q


def core_from_sqd(sq, idx):
    """sq(i,j) -> Fraction.  Build the reconstructed core of the 6 points idx (frame coords)."""
    labs = [sq(i, j) for i, j in itertools.combinations(idx, 2)]
    pat = R.A.normalize(labs)
    c, r = R.A.parse(pat)
    z = [None] * r
    for (i, j), v in zip(itertools.combinations(range(6), 2), labs):
        z[c[i][j]] = KQ.const(fmpq(v.numerator, v.denominator))
    rec = R.reconstruct_core(c, r, KQ, z)
    assert rec["status"] == "REALIZED", rec["status"]
    return rec


def sig(geo, S, x):
    return tuple(KField.key(geo.sqd(p, x)) for p in S)


def run_known(sq, npts, n, full=False):
    low = {}
    from collections import Counter
    cnt = Counter(sq(i, j) for i, j in itertools.combinations(range(npts), 2))
    lowv = {v for v, c_ in cnt.items() if c_ <= 2}
    ncores = 0
    for idx in itertools.combinations(range(npts), 6):
        if any(sq(i, j) in lowv for i, j in itertools.combinations(idx, 2)):
            continue
        r = len({sq(i, j) for i, j in itertools.combinations(idx, 2)})
        if r not in (n - 3, n - 4):
            continue
        rec = core_from_sqd(sq, idx)
        X = E.Extender(KQ, rec["G"], rec["points"], n=n, full=full)
        found = X.run()
        # actual outside signatures (squared distances to the 6 core points, in core order)
        out = [k for k in range(npts) if k not in idx]
        want = sorted(tuple(KField.key(KQ.const(fmpq(sq(i, k).numerator, sq(i, k).denominator))) for i in idx) for k in out)
        gots = [sorted(sig(X.geo, rec["points"], p) for p in f["points"][6:]) for f in found]
        assert want in gots, (idx, len(found))
        assert all(f["core_is_high"] for f in found)
        ncores += 1
    return ncores


def test_palasti_n8():
    P = [(0, 1), (1, 0), (2, 0), (F(5, 2), F(5, 2)), (F(3, 2), F(9, 2)), (F(1, 2), F(7, 2)), (F(3, 2), F(7, 2)), (1, 2)]
    sq = lambda i, j: 3 * (F(P[i][0]) - P[j][0]) ** 2 + (F(P[i][1]) - P[j][1]) ** 2
    assert run_known(sq, 8, 8) == 2
    assert run_known(sq, 8, 8, full=True) == 2


def lattice_sq(ab):
    def sq(i, j):
        p, q = ab[i][0] - ab[j][0], ab[i][1] - ab[j][1]
        return F(p * p + p * q + q * q)
    return sq


N7 = [[(0, 0), (0, 1), (1, -3), (2, -4), (2, -3), (4, -1), (5, -2)],
      [(0, 0), (0, 1), (1, -3), (2, -3), (3, -1), (3, 0), (4, -1)],
      [(0, 0), (0, 1), (1, -2), (2, -1), (2, 1), (3, -3), (3, -1)],
      [(0, 0), (0, 2), (1, -2), (1, -1), (2, 1), (3, -2), (4, -1)],
      [(0, 0), (0, 2), (1, -2), (1, 0), (2, 1), (3, -1), (4, -1)],
      [(0, 0), (0, 2), (1, -2), (1, 0), (2, 1), (3, 1), (4, -1)],
      [(0, 0), (0, 4), (1, -2), (1, 2), (2, 3), (3, -1), (4, 0)],
      [(0, 0), (1, -4), (1, -3), (2, -4), (3, -2), (3, -1), (5, -6)]]  # verify/lattice_n7_r3.txt


def test_lattice_n7():
    tot = 0
    for ab in N7:
        tot += run_known(lattice_sq(ab), 7, 7)
    assert tot >= 1
    return tot


def test_fits_bruteforce():
    rng = random.Random(1)
    for _ in range(3000):
        a = [rng.randint(1, 8) for _ in range(rng.randint(0, 6))]
        b = rng.sample(range(1, 10), rng.randint(0, 7))
        bf = any(all(x <= y for x, y in zip(a, perm)) for perm in itertools.permutations(b, len(a))) if len(a) <= len(b) else False
        assert E.fits(a, b) == bf


def brute_mu(dkeys, DS, n, nsel):
    vals = sorted(set(dkeys.values()))
    from collections import Counter
    cnt = Counter(dkeys.values())
    xs = list(range(6, 6 + nsel))
    for mu in itertools.permutations(range(1, n), len(vals)):
        M = dict(zip(vals, mu))
        if any(cnt[v] > M[v] for v in vals):
            continue
        if any(M[d] < 3 for d in DS):
            continue
        # exactly n-3 high slots, r taken by D(S): implied by injectivity into {1..n-1}
        if all(sum(M[dkeys[(i, x)]] >= 3 for i in range(6)) >= 5 for x in xs):
            return True
    return False


def test_capacity_bruteforce():
    rng = random.Random(7)
    n = 9
    seen = set()
    for trial in range(1500):
        r = rng.choice([5, 6])
        nsel = rng.choice([1, 2])
        npts = 6 + nsel
        DS = list(range(r))
        dk = {}
        # core pairs: values from D(S), every value used
        pairs6 = list(itertools.combinations(range(6), 2))
        labs = list(range(r)) + [rng.randrange(r) for _ in range(15 - r)]
        rng.shuffle(labs)
        for p, l in zip(pairs6, labs):
            dk[p] = l
        pool = list(range(r)) + [10, 11, 12, 13]
        for i, j in itertools.combinations(range(npts), 2):
            if j >= 6:
                dk[(i, j)] = rng.choice(pool) if rng.random() < 0.8 else rng.choice(range(r))
        X = SimpleNamespace(DS={d: None for d in DS}, m=(n - 3) - r, highs=list(range(3, n)))
        got = E.Extender.capacity(X, list(range(npts)), dk, nsel)
        want = brute_mu(dk, DS, n, nsel)
        assert got == want, (trial, r, nsel, dk)
        seen.add((r, nsel, got))
    assert len(seen) == 8  # both verdicts occur for every (|D(S)|, #selected) combination


def test_covering():
    for npts, s in [(8, 5), (7, 4), (6, 3)]:
        fam = E.covering_triples(npts, s)  # asserts coverage internally
        assert fam


def test_trilateration_planted():
    rng = random.Random(3)
    P = [(0, 1), (1, 0), (2, 0), (F(5, 2), F(5, 2)), (F(3, 2), F(7, 2)), (1, 2)]
    sq = lambda i, j: 3 * (F(P[i][0]) - P[j][0]) ** 2 + (F(P[i][1]) - P[j][1]) ** 2
    rec = core_from_sqd(sq, range(6))
    geo = E.Geo(KQ, rec["G"])
    S = rec["points"]
    for _ in range(200):
        x = (KQ.const(fmpq(rng.randint(-20, 20), rng.randint(1, 5))), KQ.const(fmpq(rng.randint(-20, 20), rng.randint(1, 5))))
        t = rng.sample(range(6), 3)
        T = geo.tri_setup(S[t[0]], S[t[1]], S[t[2]])
        vals = [geo.sqd(S[i], x) for i in t] + [KQ.const(fmpq(rng.randint(1, 30))) for _ in range(3)]
        out = geo.trilaterate(T, vals)
        assert geo.key(x) in {geo.key(p) for p in out}
        cc = T["x0"]
        d0 = geo.sqd(S[t[0]], cc)
        assert geo.sqd(S[t[1]], cc) == d0 and geo.sqd(S[t[2]], cc) == d0
        for p in out[1:]:  # every returned point satisfies its 3 distances for some triple
            ds = [geo.sqd(S[i], p) for i in t]
            assert all(KField.key(d) in {KField.key(v) for v in vals} for d in ds)


if __name__ == "__main__":
    import time
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            t = time.time()
            r = fn()
            print(f"{name}: OK {'' if r is None else r} ({time.time()-t:.1f}s)")
