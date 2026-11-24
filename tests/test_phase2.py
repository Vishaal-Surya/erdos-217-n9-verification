"""Builder E Phase 2 tests: pattern enumeration + affine filter.

Run:  .venv/Scripts/python -m pytest tests/test_phase2.py -q
"""
import gzip
import os
import random
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
import affine  # noqa: E402
import canon  # noqa: E402

DATA = os.path.join(ROOT, "data")


def rgs(n, maxb=6):
    """All restricted growth strings of length n with <= maxb blocks."""
    def rec(prefix, t):
        if len(prefix) == n:
            yield tuple(prefix)
            return
        for l in range(min(t + 1, maxb)):
            prefix.append(l)
            yield from rec(prefix, max(t, l + 1))
            prefix.pop()
    yield from rec([0], 1)


def brute(k):
    np_ = k * (k - 1) // 2
    out = set()
    for lab in rgs(np_):
        if canon.local_ok(lab, k):
            out.add(canon.canonical(lab, k))
    return out


def read_lines(name):
    path = os.path.join(DATA, name)
    op = gzip.open if name.endswith(".gz") else open
    with op(path, "rt") as fh:
        return [ln.strip() for ln in fh if ln.strip()]


@pytest.mark.parametrize("k,expected", [(3, 3), (4, 22), (5, 932)])
def test_small_counts_python_bruteforce_matches_java(k, expected):
    got = brute(k)
    assert len(got) == expected
    assert got == set(read_lines(f"patterns{k}_local.txt"))


def test_six_point_file_counts():
    assert len(read_lines("patterns6_local.txt.gz")) == 625472
    assert len(read_lines("patterns6_affine.txt.gz")) == 344835


def test_java_canonical_forms_agree_with_python_sample():
    pats = read_lines("patterns6_local.txt.gz")
    rnd = random.Random(217)
    sample = rnd.sample(pats, 1500)
    for s in sample:
        assert canon.canonical(s, 6) == s
        assert canon.local_ok(s, 6)
        # random relabel + permutation maps back to the same canonical form
        perm = list(range(6)); rnd.shuffle(perm)
        ren = list("abcdef"); rnd.shuffle(ren)
        c = canon.matrix(s, 6)
        lab = [ren[int(c[perm[i]][perm[j]]) - 1] for i in range(6) for j in range(i + 1, 6)]
        assert canon.canonical(lab, 6) == s


def test_affine_engines_agree_sample():
    pats = read_lines("patterns6_local.txt.gz")
    rnd = random.Random(1)
    for s in rnd.sample(pats, 3000):
        assert affine.affine_test(s, "flint")[0] == affine.affine_test(s, "fraction")[0]


def test_affine_survivors_subset():
    loc = set(read_lines("patterns6_local.txt.gz"))
    aff = read_lines("patterns6_affine.txt.gz")
    assert all(s in loc for s in aff)
    rnd = random.Random(2)
    for s in rnd.sample(aff, 2000):
        assert affine.affine_test(s, "fraction")[0]


# ---- hand examples for the local rules ----
def test_local_rules_hand_examples():
    # K4 all-equal labels: excluded by Lemma 2.2 (all three matchings equal)
    assert not canon.local_ok([0] * 6, 4)
    # (i): 4 equal pairs at a point (star K_{1,4} on 5 points); other pairs distinct
    lab = [0, 0, 0, 0, 1, 2, 3, 4, 5, 1]  # pairs 12 13 14 15 | 23 24 25 34 35 45
    assert not canon.local_ok(lab, 5)
    # (iii): three points x with c_1x = c_2x
    #   pairs 12 13 14 15 23 24 25 34 35 45 ; c13=c23=a, c14=c24=b, c15=c25=c
    lab = [0, 1, 2, 3, 1, 2, 3, 4, 5, 4]
    assert not canon.local_ok(lab, 5)


def test_R_and_A_hold_on_realizable_examples():
    """For every stored realizable example, the affine rows R annihilate the exact lattice
    coordinates (affine relations are basis-independent) and A z = 0 at the true squared
    distances. Direct check that R and A are built correctly (not just that the filter passes)."""
    import ast
    import glob
    paths = sorted(glob.glob(os.path.join(DATA, "soundness_realizable6_*.txt.gz")))
    if not paths:
        pytest.skip("run tests/soundness_run.py first")
    from soundness_run import qdist, PALASTI_FORM  # noqa
    nrows = 0
    aff = set(read_lines("patterns6_affine.txt.gz"))
    lines = []
    for path in paths:
        with gzip.open(path, "rt") as fh:
            lines += list(fh)
    if True:
        for ln in lines:
            cpat, src, pts = ln.rstrip("\n").split("\t")
            f = PALASTI_FORM if src == "palasti" else ast.literal_eval(src)
            pts = ast.literal_eval(pts)
            vals = [qdist(f, pts[i], pts[j]) for i in range(6) for j in range(i + 1, 6)]
            distinct = sorted(set(vals))
            s = "".join(str(distinct.index(v) + 1) for v in vals)
            c = affine.to_matrix(s)
            R = affine.build_R(c)
            for rho in R:
                assert sum(rho[i] * pts[i][0] for i in range(6)) == 0
                assert sum(rho[i] * pts[i][1] for i in range(6)) == 0
                assert sum(rho) == 0
            A = affine.build_A(c, R, len(distinct))
            for row in A:
                assert sum(a * z for a, z in zip(row, distinct)) == 0
            nrows += len(R)
            assert affine.affine_test(s, "fraction")[0]
            assert canon.canonical(vals, 6) == cpat and cpat in aff
    assert nrows > 0


def test_realizable_palasti_subsets_survive():
    import itertools
    from soundness_run import PALASTI_UV, PALASTI_FORM, qdist, general_position, pattern_of  # noqa
    aff = set(read_lines("patterns6_affine.txt.gz"))
    n = 0
    for sub in itertools.combinations(PALASTI_UV, 6):
        labs = pattern_of(PALASTI_FORM, list(sub))
        if len(set(labs)) <= 6 and general_position(PALASTI_FORM, list(sub)):
            assert canon.canonical(labs, 6) in aff
            n += 1
    assert n > 0


sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
