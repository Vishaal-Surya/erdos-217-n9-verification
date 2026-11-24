import itertools
import json
import os
import random
import subprocess
import sys
from fractions import Fraction

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

import checker as C  # noqa: E402
import lattice_search as L  # noqa: E402

PALASTI = {"field": {"type": "sqrt", "d": 3},
           "points": [["0", "1"], ["sqrt(3)", "0"], ["2*sqrt(3)", "0"], ["5*sqrt(3)/2", "5/2"],
                      ["3*sqrt(3)/2", "9/2"], ["sqrt(3)/2", "7/2"], ["3*sqrt(3)/2", "7/2"],
                      ["sqrt(3)", "2"]]}


def spec_with(points, base=PALASTI):
    return {"field": base["field"], "points": points}


def run(spec):
    return C.check_spec(spec)


# ---------------------------------------------------------------- number field

def test_field_arithmetic_exact():
    K = C.NumberField.sqrt(3)
    s = K.gen()
    assert s * s == K(3)
    assert (s + 1) * (s - 1) == K(2)
    assert not (s - K(Fraction(17320508, 10000000))).is_zero()
    with pytest.raises(ValueError):
        C.NumberField.sqrt(12)
    with pytest.raises(ValueError):
        C.NumberField([-4, 0, 1])  # x^2 - 4 reducible


def test_sqrt_parsing_normalises():
    K, pts = C.load_points({"field": {"type": "sqrt", "d": 3},
                            "points": [["sqrt(12)", "s/2 + 1/3"], ["2*s", "(sqrt(3)+2/3)/2"]]})
    assert pts[0][0] == pts[1][0]
    assert pts[0][1] == pts[1][1]
    with pytest.raises(Exception):
        C.load_points({"field": {"type": "sqrt", "d": 3}, "points": [["sqrt(2)", "0"]]})


# ---------------------------------------------------------------- ACCEPT

def test_accept_palasti_8():
    r = run(PALASTI)
    assert r.ok, r.summary()
    assert r.multiplicities == list(range(1, 8))


def test_accept_palasti_8_via_minpoly_and_algebraic_modes():
    # same configuration, field given as Q(a), a^2 = 3
    pts = [[x.replace("sqrt(3)", "a"), y] for x, y in PALASTI["points"]]
    r = run({"field": {"type": "minpoly", "minpoly": "a**2 - 3", "var": "a"}, "points": pts})
    assert r.ok
    r = run({"field": {"type": "algebraic"}, "points": PALASTI["points"]})
    assert r.ok


def test_accept_n3_isosceles():
    r = run({"field": "Q", "points": [["0", "0"], ["2", "0"], ["1", "3"]]})
    # d^2: 4, 10, 10
    assert r.ok and r.multiplicities == [1, 2]


def test_reject_n3_equilateral_and_scalene():
    assert not run({"field": {"type": "sqrt", "d": 3},
                    "points": [["0", "0"], ["2", "0"], ["1", "sqrt(3)"]]}).ok
    assert not run({"field": "Q", "points": [["0", "0"], ["3", "0"], ["1", "2"]]}).ok


def test_accept_n4_hand_example():
    # Rhombus of two equilateral triangles: d^2 = 4 x5, 12 x1 -> NOT crescent.
    r = run({"field": {"type": "sqrt", "d": 3},
             "points": [["0", "0"], ["2", "0"], ["1", "sqrt(3)"], ["3", "sqrt(3)"]]})
    assert not r.ok and r.multiplicities == [1, 5]
    # Hand-verified n=4 crescent: A(0,0) B(2,0) C(1,sqrt3) equilateral (d^2 = 4, x3);
    # D(1,-1): AD^2 = BD^2 = 2 (x2), CD^2 = 4 + 2 sqrt3 (x1).  No 3 collinear; D is at distance
    # 1 + 1/sqrt3 from the circumcentre (1, 1/sqrt3) of ABC, not the radius 2/sqrt3.
    r = run({"field": {"type": "sqrt", "d": 3},
             "points": [["0", "0"], ["2", "0"], ["1", "sqrt(3)"], ["1", "-1"]]})
    assert r.ok, r.summary()
    assert r.multiplicities == [1, 2, 3]
    assert r.classes[0]["sqdist"].replace(" ", "") in ("2*sqrt(3)+4", "4+2*sqrt(3)")


def test_accept_general_algebraic_field():
    # n=3 isosceles with apex height cbrt(2): d^2 = 4, 1+2^(2/3) (x2)
    r = run({"field": {"type": "algebraic"}, "points": [["0", "0"], ["2", "0"], ["1", "2**(1/3)"]]})
    assert r.ok
    r = run({"field": {"type": "minpoly", "minpoly": "a**3 - 2", "var": "a"},
             "points": [["0", "0"], ["2", "0"], ["1", "a"]]})
    assert r.ok
    # mixed radicals: field Q(sqrt2, cbrt3) of degree 6
    r = run({"field": {"type": "algebraic"},
             "points": [["0", "0"], ["2", "0"], ["1", "sqrt(2) + 3**(1/3)"]]})
    assert r.ok
    # equilateral written with nested radicals must be detected exactly
    r = run({"field": {"type": "algebraic"},
             "points": [["0", "0"], ["2", "0"], ["1", "sqrt(4 - 2*sqrt(3)) + 1"]]})
    assert not r.ok and r.multiplicities == [3]


# ---------------------------------------------------------------- REJECT

def test_reject_three_collinear():
    pts = [list(p) for p in PALASTI["points"]]
    # move point 7 (sqrt3, 2) to (sqrt3, 1/2): collinear? Use exact: put it on line through
    # points 1=(sqrt3,0) and 2=(2sqrt3,0): y=0 -> (3*sqrt3, 0)
    pts[7] = ["3*sqrt(3)", "0"]
    r = run(spec_with(pts))
    assert not r.ok
    assert any("collinear triple: (1, 2, 7)" in x for x in r.reasons)


def test_reject_four_concyclic():
    # square-ish: four points on the unit circle (rational points) + others
    pts = [["1", "0"], ["0", "1"], ["-1", "0"], ["3/5", "4/5"], ["7", "3"]]
    r = run({"field": "Q", "points": pts})
    assert not r.ok
    assert any("concyclic quadruple: (0, 1, 2, 3)" in x for x in r.reasons)


def test_reject_four_concyclic_in_sqrt_field():
    # Palasti with point 0 replaced by a point on the circumcircle of points 1,2,3
    # circle through (s,0),(2s,0),(5s/2,5/2): centre x = 3s/2; reflect (5s/2,5/2) in x=3s/2 -> (s/2, 5/2)
    pts = [list(p) for p in PALASTI["points"]]
    pts[0] = ["sqrt(3)/2", "5/2"]
    r = run(spec_with(pts))
    assert not r.ok
    assert any(x.startswith("concyclic quadruple: (0, 1, 2, 3)") for x in r.reasons)


def test_reject_wrong_multiplicities():
    pts = [list(p) for p in PALASTI["points"]]
    pts[0] = ["1/7", "1"]  # rational offset of (0,1)
    r = run(spec_with(pts))
    assert not r.ok
    assert any("distance class sizes" in x for x in r.reasons)


def test_reject_duplicate():
    pts = [list(p) for p in PALASTI["points"]] + [["sqrt(12)/2", "2"]]
    r = run(spec_with(pts))
    assert not r.ok
    assert any("duplicate point: 8 == 7" in x for x in r.reasons)


# ---------------------------------------------------------------- cross-checks

def brute_classes(pts):
    """Independent brute force over Fractions: pairwise equality comparisons."""
    pairs = list(itertools.combinations(range(len(pts)), 2))
    d = {(i, j): (pts[i][0] - pts[j][0]) ** 2 + (pts[i][1] - pts[j][1]) ** 2 for i, j in pairs}
    classes = []
    for pr in pairs:
        for cl in classes:
            if d[cl[0]] == d[pr]:
                cl.append(pr)
                break
        else:
            classes.append([pr])
    return sorted((str(d[cl[0]]), sorted(cl)) for cl in classes)


@pytest.mark.parametrize("seed", range(40))
def test_random_rational_class_table(seed):
    rng = random.Random(seed)
    n = rng.randint(3, 9)
    pts = [(Fraction(rng.randint(-4, 4), rng.choice([1, 2])), Fraction(rng.randint(-4, 4), rng.choice([1, 2])))
           for _ in range(n)]
    K = C.NumberField.rationals()
    epts = [(K(x), K(y)) for x, y in pts]
    groups = C.distance_classes(epts)
    mine = sorted((str(v[0]), sorted(v[1])) for v in groups.values())
    assert mine == brute_classes(pts)
    # general-position predicates vs Fraction brute force
    r = C.check(epts)
    for i, j, k in itertools.combinations(range(n), 3):
        (ax, ay), (bx, by), (cx, cy) = pts[i], pts[j], pts[k]
        col = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax) == 0
        assert col == (f"collinear triple: {(i, j, k)}" in r.reasons) or len(r.reasons) >= 50


# ---------------------------------------------------------------- lattice search

def lattice_pts_to_field(conf):
    K = C.NumberField.sqrt(3)
    s = K.gen()
    return [(K(Fraction(2 * a + b, 2)), s * Fraction(b, 2)) for a, b in conf]


def test_lattice_predicates_match_generic_checker():
    rng = random.Random(1)
    P = L.patch(3)
    for _ in range(3000):
        quad = rng.sample(P, 4)
        F = lattice_pts_to_field(quad)
        assert L.collinear(*quad[:3]) == C.orient(*F[:3]).is_zero()
        if not any(L.collinear(*t) for t in itertools.combinations(quad, 3)):
            assert L.concyclic(*quad) == C.incircle(*F).is_zero()
        assert L.norm(quad[0][0] - quad[1][0], quad[0][1] - quad[1][1]) == C.sqdist(F[0], F[1])
    # deliberately concyclic: hexagon vertices
    hexa = [(1, 0), (0, 1), (-1, 1), (-1, 0)]
    assert L.concyclic(*hexa) and C.incircle(*lattice_pts_to_field(hexa)).is_zero()


def test_lattice_search_finds_valid_n7():
    sols, _, _, _ = L.search(3, 7, find_all=False)
    assert sols
    assert C.check(lattice_pts_to_field(sols[0])).ok


def test_known_examples_file_all_accepted():
    with open(os.path.join(ROOT, "data", "known_examples.json"), encoding="utf-8") as fh:
        rec = json.load(fh)
    count = 0
    for k, res in rec["results"].items():
        for ex in res["examples"]:
            r = C.check_spec(ex["checker_input"])
            assert r.ok and r.multiplicities == list(range(1, int(k))), (k, ex)
            count += 1
    assert "7" in rec["results"] and rec["results"]["7"]["examples"]
    assert count > 0


def test_cli_runs(tmp_path):
    f = tmp_path / "p.json"
    f.write_text(json.dumps(PALASTI))
    out = subprocess.run([sys.executable, os.path.join(ROOT, "src", "checker.py"), str(f)],
                         capture_output=True, text=True)
    assert out.returncode == 0 and "CRESCENT (ok)" in out.stdout
    bad = dict(PALASTI, points=PALASTI["points"][:-1] + [["0", "1"]])
    f.write_text(json.dumps(bad))
    out = subprocess.run([sys.executable, os.path.join(ROOT, "src", "checker.py"), str(f)],
                         capture_output=True, text=True)
    assert out.returncode == 1
