"""Build known-answer targets for the coverage-search tests (Builder C).

From exact coordinates (data/known_examples.json triangular-lattice n=7,8 examples and Palasti n=8)
compute, in exact integer/rational arithmetic:
  * the class-size string (pair (i,j) in lex order -> multiplicity of its squared distance),
    checking it is a crescent pattern (sizes 1..n-1, general position);
  * the low-pair arrangement name;
  * the k-point (k=n-3) all-high core patterns (canonical, via src/canon.py), to build minimal U files.
Writes data/coverage_targets_n{7,8}.txt and data/coverage_coresU_n{7,8}*.txt.
"""
import itertools
import json
import os
import sys
from fractions import Fraction as Fr

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "src", "coverage"))
import canon  # noqa: E402
import arrangements  # noqa: E402

DATA = os.path.join(ROOT, "data")


def tri(ab):  # triangular lattice point (a,b) -> exact (x, y/sqrt3-free) handled via quadratic form
    return tuple(ab)


def d2_tri(p, q):
    a, b = p[0] - q[0], p[1] - q[1]
    return a * a + a * b + b * b


def gp_ok(pts, d2):
    n = len(pts)
    # collinearity / concyclicity via Cayley-Menger-type tests on squared distances (exact):
    # 3 points collinear iff 16*area^2 = 0 ; 4 concyclic iff Ptolemy det = 0 (det of [d2] 4x4 = 0 for concyclic)
    for i, j, k in itertools.combinations(range(n), 3):
        a, b, c = d2(i, j), d2(i, k), d2(j, k)
        if 2 * (a * b + a * c + b * c) - (a * a + b * b + c * c) == 0:
            return False
    for q in itertools.combinations(range(n), 4):
        M = [[Fr(d2(x, y)) if x != y else Fr(0) for y in q] for x in q]
        if det(M) == 0:
            return False
    return True


def det(M):
    M = [r[:] for r in M]
    n = len(M)
    s = Fr(1)
    for i in range(n):
        piv = next((r for r in range(i, n) if M[r][i] != 0), None)
        if piv is None:
            return Fr(0)
        if piv != i:
            M[i], M[piv] = M[piv], M[i]
            s = -s
        s *= M[i][i]
        for r in range(i + 1, n):
            f = M[r][i] / M[i][i]
            for c in range(i, n):
                M[r][c] -= f * M[i][c]
    return s


def class_string(n, d2):
    vals = [d2(i, j) for i, j in itertools.combinations(range(n), 2)]
    mult = {v: vals.count(v) for v in set(vals)}
    sizes = sorted(mult.values())
    assert sizes == list(range(1, n)), sizes
    return "".join(str(mult[v]) for v in vals)


def arrangement_name(s, n):
    pairs = list(itertools.combinations(range(n), 2))
    e1 = [p for p, ch in zip(pairs, s) if ch == "1"]
    e2 = [p for p, ch in zip(pairs, s) if ch == "2"]
    c = arrangements.canon_of((e1[0], e2[0], e2[1]))
    for name, rep in arrangements.OURS.items():
        r0 = tuple(tuple(x - 1 - (9 - n) for x in e) for e in rep)
        if arrangements.canon_of(r0) == c:
            return name
    raise AssertionError


def cores(s, n):
    k = n - 3
    pairs = list(itertools.combinations(range(n), 2))
    lab = dict(zip(pairs, s))
    out = set()
    for T in itertools.combinations(range(n), k):
        sub = [lab[(a, b)] for a, b in itertools.combinations(T, 2)]
        if any(x in "12" for x in sub):
            continue
        out.add(canon.canonical(sub, k))
    return out


def main():
    d = json.load(open(os.path.join(DATA, "known_examples.json")))
    ex = {7: [], 8: []}
    for n in (7, 8):
        for t, e in enumerate(d["results"][str(n)]["examples"]):
            pts = [tuple(p) for p in e["lattice_ab"]]
            dd = lambda i, j, pts=pts: d2_tri(pts[i], pts[j])
            assert gp_ok(pts, dd)
            ex[n].append((f"tri{n}_{t}", class_string(n, dd)))
    # Palasti n=8: points (r*sqrt3, s) with rational r,s ; squared distance 3 dr^2 + ds^2
    pal = [(0, 1), (1, 0), (2, 0), (Fr(5, 2), Fr(5, 2)), (Fr(3, 2), Fr(9, 2)), (Fr(1, 2), Fr(7, 2)), (Fr(3, 2), Fr(7, 2)), (1, 2)]
    pal = [(Fr(a), Fr(b)) for a, b in pal]
    dp = lambda i, j: 3 * (pal[i][0] - pal[j][0]) ** 2 + (pal[i][1] - pal[j][1]) ** 2
    assert gp_ok(pal, dp)
    ex[8].append(("palasti8", class_string(8, dp)))
    for n in (7, 8):
        with open(os.path.join(DATA, f"coverage_targets_n{n}.txt"), "w") as f:
            f.write(f"# name class-size-string(lex pairs) arrangement canonical(min over {n}! perms, python)\n")
            allc = set()
            for name, s in ex[n]:
                cf = min("".join(s[pairs_idx] for pairs_idx in pm) for pm in canon._perm_maps(n))
                f.write(f"{name} {s} {arrangement_name(s, n)} {cf}\n")
                cs = cores(s, n)
                allc |= cs
                with open(os.path.join(DATA, f"coverage_coresU_n{n}_{name}.txt"), "w") as g:
                    g.write("\n".join(sorted(cs)) + "\n")
                print(n, name, s, arrangement_name(s, n), "cores:", len(cs))
            with open(os.path.join(DATA, f"coverage_coresU_n{n}_all.txt"), "w") as g:
                g.write("\n".join(sorted(allc)) + "\n")


if __name__ == "__main__":
    main()
