"""Soundness tests for src/algebra.py (Builder A).

Mandatory: NEITHER engine may exclude a realizable pattern (exact point sets in general
position with <= 6 distinct distances, data/algebra_realizable.json, built by
src/algebra_testgen.build_cache).  Plus: local tests (Lemma 2.1/2.2) implemented here
independently, plausibility on random label strings, E1-exclusion => E2-exclusion
consistency, certificate replay, invariance under point relabelling.

Run:  .venv/Scripts/python -m pytest tests/test_algebra.py -q
Set ALGEBRA_FULL=1 to test every cached realizable pattern (default: a 300-pattern sample).
"""
import itertools
import json
import os
import random
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
import algebra as A  # noqa: E402
import algebra_testgen as G  # noqa: E402

CACHE = os.path.join(ROOT, "data", "algebra_realizable.json")
PAIRS = [(i, j) for i in range(6) for j in range(i + 1, 6)]


def _realizable_patterns():
    if not os.path.exists(CACHE):
        G.build_cache(CACHE)
    with open(CACHE) as f:
        data = json.load(f)
    pats = sorted({s for v in data.values() for s, _ in v})
    if os.environ.get("ALGEBRA_FULL"):
        return pats
    rng = random.Random(3)
    return sorted(rng.sample(pats, min(300, len(pats))))


# --- independent local tests (Lemma 2.1 (i)-(iii), Lemma 2.2), written separately from algebra.local_ok


def local_tests(s):
    lab = {}
    for (i, j), ch in zip(PAIRS, s):
        lab[(i, j)] = lab[(j, i)] = ch
    for a in range(6):
        vals = [lab[(a, x)] for x in range(6) if x != a]
        if any(vals.count(v) >= 4 for v in vals):
            return False
    for a, b in itertools.combinations(range(6), 2):
        xs = [x for x in range(6) if x not in (a, b)]
        if sum(lab[(a, x)] == lab[(b, x)] for x in xs) >= 3:
            return False
        keys = [(lab[(a, x)], lab[(b, x)]) for x in xs]
        if any(keys.count(k) >= 3 for k in keys):
            return False
    for q in itertools.combinations(range(6), 4):
        a, b, c, d = q
        if lab[(a, b)] == lab[(c, d)] and lab[(a, c)] == lab[(b, d)] and lab[(a, d)] == lab[(b, c)]:
            return False
    return True


def random_pattern(rng, rmax=6):
    while True:
        s = A.normalize([rng.randint(1, rmax) for _ in range(15)])
        if max(int(ch) for ch in s) <= 6:
            return s


def permuted(s, perm):
    c = {}
    for (i, j), ch in zip(PAIRS, s):
        c[(i, j)] = c[(j, i)] = ch
    return A.normalize([c[(perm[i], perm[j])] for i, j in PAIRS])


REAL = _realizable_patterns()


def test_cache_nonempty_and_local():
    assert len(REAL) >= 100
    for s in REAL:
        assert local_tests(s), s
        assert A.local_ok(A.parse(s)[0]), s


@pytest.mark.parametrize("s", REAL)
def test_engine1_never_excludes_realizable(s):
    r = A.engine1(s, check_local=True)
    assert r["verdict"] in ("F", "U"), (s, r.get("stage"), r.get("reason"))
    # mod-p and Q agree on these
    assert not r.get("modp_discrepancy"), s


@pytest.mark.parametrize("s", REAL[:120])
def test_engine2_never_excludes_realizable(s):
    r = A.engine2(s, eliminants=False)
    assert r["verdict"] != "EXCLUDED", (s, r.get("reason"))


def test_local_tests_agree():
    rng = random.Random(11)
    for _ in range(3000):
        s = random_pattern(rng)
        assert local_tests(s) == A.local_ok(A.parse(s)[0]), s


def test_random_strings_mostly_excluded():
    rng = random.Random(5)
    n = ex = 0
    for _ in range(150):
        s = random_pattern(rng)
        r = A.engine1(s, check_local=True)
        n += 1
        ex += r["verdict"] == "EXCLUDED"
    assert ex / n > 0.9, (ex, n)


def test_relabel_invariance():
    rng = random.Random(2)
    pats = REAL[:15]
    while len(pats) < 40:
        s = random_pattern(rng)
        if local_tests(s):
            pats.append(s)
    for s in pats:
        v0 = A.engine1(s)
        for _ in range(3):
            perm = list(range(6)); rng.shuffle(perm)
            t = permuted(s, perm)
            v1 = A.engine1(t)
            assert (v0["verdict"], v0.get("stage") or "") == (v1["verdict"], v1.get("stage") or ""), (s, t)
            if v0["verdict"] != "EXCLUDED":
                assert v0["h"] == v1["h"]


def test_certificates_replay():
    rng = random.Random(9)
    k = 0
    while k < 25:
        s = random_pattern(rng)
        if not local_tests(s):
            continue
        r = A.engine1(s)
        if r.get("stage") == "quartic":
            assert A.verify_certificate(s, r["cert"]), s
            k += 1


def test_e1_exclusion_implies_e2_not_feasible():
    """E1 exclusions are consequences of the ideal + nonzero factors, so E2 must never report
    ZERODIM/POSDIM for an E1-excluded pattern (UNRESOLVED allowed)."""
    rng = random.Random(13)
    k = 0
    while k < 40:
        s = random_pattern(rng)
        if not local_tests(s):
            continue
        r1 = A.engine1(s)
        if r1["verdict"] != "EXCLUDED":
            continue
        k += 1
        r2 = A.engine2(s, eliminants=False)
        assert r2["verdict"] in ("EXCLUDED", "UNRESOLVED"), (s, r1["reason"], r2)


def test_palasti_subsets():
    for s, _ in G.gen_palasti():
        assert A.engine1(s)["verdict"] in ("F", "U")
        assert A.engine2(s, eliminants=False)["verdict"] != "EXCLUDED"
