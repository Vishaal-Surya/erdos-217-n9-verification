"""Verifier test cases for verify/checker_v.py. Prints PASS/FAIL per case."""
import sys
import os
from fractions import Fraction as Fr
from itertools import combinations, product

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from checker_v import (check, Q3Field, QField, q3pt, circumcentre, Q3Field as F3)

results = []


def case(name, ok_expected, pts, F=Q3Field, reason_prefix=None):
    ok, info = check(pts, F)
    good = (ok == ok_expected)
    if good and reason_prefix is not None:
        good = info["reason"].startswith(reason_prefix)
    results.append(good)
    print(f"{'PASS' if good else 'FAIL'}  [{F.name}] {name}: ok={ok} reason={info['reason']}"
          + (f" sizes={info.get('class_sizes')}" if 'class_sizes' in info else ""))
    return ok, info


h = Fr(1, 2)
PALASTI = [q3pt(0, 1), q3pt((0, 1), 0), q3pt((0, 2), 0), q3pt((0, Fr(5, 2)), Fr(5, 2)),
           q3pt((0, Fr(3, 2)), Fr(9, 2)), q3pt((0, h), Fr(7, 2)),
           q3pt((0, Fr(3, 2)), Fr(7, 2)), q3pt((0, 1), 2)]

# ---- ACCEPT
case("Palasti n=8 as printed", True, PALASTI)


def rot(pts, c, s, F=F3):
    return [(F.sub(F.mul(c, x), F.mul(s, y)), F.add(F.mul(s, x), F.mul(c, y))) for x, y in pts]


def shift(pts, dx, dy, F=F3):
    return [(F.add(x, dx), F.add(y, dy)) for x, y in pts]


case("Palasti rotated 30deg + shifted", True,
     shift(rot(PALASTI, F3.c(0, h), F3.c(h)), F3.c(7, -2), F3.c(Fr(1, 3), 5)))
case("Palasti rotated 60deg", True, rot(PALASTI, F3.c(h), F3.c(0, h)))
case("Palasti mirrored (x->-x)", True, [(F3.sub(F3.zero, x), y) for x, y in PALASTI])
for k in range(8):
    sub7 = PALASTI[:k] + PALASTI[k + 1:]
    ok, info = check(sub7)
    print(f"INFO  Palasti minus point {k}: ok={ok} sizes={info.get('class_sizes')}")

# ---- REJECT (Q3)
case("Palasti with duplicate point", False, PALASTI[:7] + [PALASTI[0]], reason_prefix="duplicate")
case("Palasti, last point -> (3sqrt3,0) [collinear w/ (sqrt3,0),(2sqrt3,0)]", False,
     PALASTI[:7] + [q3pt((0, 3), 0)], reason_prefix="collinear")
case("Palasti, last point perturbed to (sqrt3, 15/7)", False,
     PALASTI[:7] + [q3pt((0, 1), Fr(15, 7))])
case("Palasti first 7 + 8th point (n=8 truncated to 7 -> wrong mults)", False, PALASTI[:7])


# 4-concyclic variant: replace point k by the antipode of some other point on the
# circumcircle of three Palasti points; accept the first variant whose rejection reason
# is 'concyclic' (i.e. general-position-in-lines but 4 on a circle).
def antipode_variants():
    for (a, b, c) in combinations(range(8), 3):
        O = circumcentre(F3, PALASTI[a], PALASTI[b], PALASTI[c])
        for src in (a, b, c):
            D = (F3.sub(F3.mul(F3.c(2), O[0]), PALASTI[src][0]),
                 F3.sub(F3.mul(F3.c(2), O[1]), PALASTI[src][1]))
            for k in range(8):
                if k in (a, b, c):
                    continue
                pts = PALASTI[:k] + [D] + PALASTI[k + 1:]
                ok, info = check(pts)
                if info["reason"].startswith("concyclic"):
                    return (a, b, c, src, k), pts
    return None, None


tag, pts = antipode_variants()
if pts is None:
    print("FAIL  could not construct concyclic Palasti variant"); results.append(False)
else:
    case(f"Palasti with 4 concyclic (circle of {tag[:3]}, replaced pt {tag[4]})", False, pts,
         reason_prefix="concyclic")

# equilateral triangle (n=3, sizes [3]) in Q3
case("equilateral triangle", False, [q3pt(0, 0), q3pt(2, 0), q3pt(1, (0, 1))],
     reason_prefix="1 distinct")
# n=3 crescent (isosceles, non-equilateral) in Q3
case("isosceles triangle (0,0),(2,0),(1,sqrt3+1)", True, [q3pt(0, 0), q3pt(2, 0), q3pt(1, (1, 1))])

# ---- RATIONAL MODE
R = lambda x, y: (Fr(x), Fr(y))
case("square (4 concyclic)", False, [R(0, 0), R(1, 0), R(0, 1), R(1, 1)], F=QField,
     reason_prefix="concyclic")
case("isosceles triangle", True, [R(0, 0), R(2, 0), R(1, 3)], F=QField)
case("3 collinear", False, [R(0, 0), R(1, 1), R(3, 3), R(0, 5)], F=QField, reason_prefix="collinear")
case("duplicate", False, [R(0, 0), R(2, 0), R(0, 0)], F=QField, reason_prefix="duplicate")


# small integer-grid search: first n=4 crescent, and first n=4 general-position
# config with exactly 3 distances but sizes != [1,2,3]; cross-check Q vs Q3 modes.
grid = [R(x, y) for x in range(5) for y in range(5)]
found_ok, found_bad = None, None
agree = True
cnt = 0
for quad in combinations(grid, 4):
    okq, iq = check(list(quad), QField)
    pts3 = [(F3.c(x), F3.c(y)) for x, y in quad]
    ok3, i3 = check(pts3, Q3Field)
    cnt += 1
    if okq != ok3 or iq["reason"] != i3["reason"]:
        agree = False
    if okq and found_ok is None:
        found_ok = quad
    if (not okq and iq.get("distinct") == 3 and found_bad is None):
        found_bad = quad
print(f"{'PASS' if agree else 'FAIL'}  Q vs Q(sqrt3) mode agree on all {cnt} 4-subsets of 5x5 grid")
results.append(agree)
if found_ok:
    case(f"rational n=4 crescent {[(int(x), int(y)) for x, y in found_ok]}", True, list(found_ok), F=QField)
if found_bad:
    case(f"rational n=4, 3 distances, wrong sizes {[(int(x), int(y)) for x, y in found_bad]}",
         False, list(found_bad), F=QField, reason_prefix="class sizes")

print()
print(f"SUMMARY: {sum(results)}/{len(results)} PASS")
sys.exit(0 if all(results) else 1)
