"""Derive the low-pair arrangements from scratch (Builder C, coverage search).

A low structure on n points = (e1, {e2, e2'}): three distinct unordered pairs, e1 carrying delta1
(multiplicity 1) and the two others carrying delta2 (multiplicity 2). Two structures are equivalent
iff some permutation of the n points maps one to the other (e1 -> e1, {e2,e2'} -> {e2,e2'}).

We enumerate ALL structures on n points, compute orbits by brute force (canonical form = minimum
over all n! permutations of the sorted encoding), and check that the representatives used by
CoverageSearch.java (ours) and by Wang's Table (theirs) each hit every orbit exactly once.
Also reports, for each representative, the largest k such that points 1..k span no low pair.

Usage: python arrangements.py [n=9]
"""
import itertools
import sys


def enc(e1, e2s):
    return (tuple(sorted(e1)), tuple(sorted(tuple(sorted(e)) for e in e2s)))


def canon(e1, e2s, n, perms):
    best = None
    for p in perms:
        a = tuple(sorted((p[e1[0]], p[e1[1]])))
        bs = tuple(sorted(tuple(sorted((p[x], p[y]))) for x, y in e2s))
        key = (a, bs)
        if best is None or key < best:
            best = key
    return best


def orbits(n):
    pairs = list(itertools.combinations(range(n), 2))
    perms = list(itertools.permutations(range(n)))
    # speed: canonical form only depends on the <= 6 involved vertices; relabel them to 0..m-1 first
    seen = {}
    for e1 in pairs:
        rest = [q for q in pairs if q != e1]
        for e2s in itertools.combinations(rest, 2):
            verts = sorted(set(e1) | set(e2s[0]) | set(e2s[1]))
            m = {v: i for i, v in enumerate(verts)}
            k = len(verts)
            sub_perms = PERMS.setdefault(k, list(itertools.permutations(range(k))))
            c = canon((m[e1[0]], m[e1[1]]), [(m[x], m[y]) for x, y in e2s], k, sub_perms)
            seen.setdefault(c, 0)
            seen[c] += 1
    return seen


PERMS = {}


def canon_of(rep):
    e1, e2a, e2b = rep
    verts = sorted(set(e1) | set(e2a) | set(e2b))
    m = {v: i for i, v in enumerate(verts)}
    k = len(verts)
    sub_perms = PERMS.setdefault(k, list(itertools.permutations(range(k))))
    return canon((m[e1[0]], m[e1[1]]), [(m[x], m[y]) for x, y in (e2a, e2b)], k, sub_perms)


def prefix_high(rep, n):
    k = 0
    for t in range(1, n + 1):
        if all(not (max(e) < t) for e in rep):  # 0-based: pair inside 0..t-1
            k = t
    return k


# 1-based, n = 9.  (delta1, delta2, delta2)
OURS = {
    "disjoint":           ((4, 7), (5, 8), (6, 9)),
    "triangle":           ((8, 9), (7, 8), (7, 9)),
    "star":               ((6, 9), (7, 9), (8, 9)),
    "chain_d1_middle":    ((7, 9), (6, 9), (7, 8)),
    "chain_d1_end":       ((6, 9), (7, 9), (7, 8)),
    "P3K2_d1_disjoint":   ((5, 8), (6, 9), (7, 9)),
    "P3K2_d1_meeting":    ((6, 9), (7, 9), (5, 8)),
}
WANG = {
    "disjoint": ((4, 7), (5, 8), (6, 9)),
    "triangle": ((8, 9), (1, 8), (1, 9)),
    "star": ((1, 9), (2, 9), (3, 9)),
    "chain_d1_middle": ((8, 9), (1, 8), (2, 9)),
    "chain_d1_end": ((1, 8), (8, 9), (2, 9)),
    "P3K2_d1_disjoint": ((3, 9), (1, 8), (2, 8)),
    "P3K2_d1_meeting": ((1, 8), (2, 8), (3, 9)),
}


def shift(rep, n):
    return tuple(tuple(x - (9 - n) for x in e) for e in rep)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    orb = orbits(n)
    total = sum(orb.values())
    print(f"n={n}: {total} labelled low structures, {len(orb)} orbits; orbit sizes {sorted(orb.values())}")
    for name, table in (("ours", OURS), ("wang", WANG)):
        if name == "wang" and n != 9:
            continue
        hit = {}
        for key, rep in table.items():
            r = shift(rep, n) if name == "ours" else rep
            r0 = tuple(tuple(x - 1 for x in e) for e in r)
            assert all(0 <= x < n for e in r0 for x in e), (key, r)
            c = canon_of(r0)
            assert c in orb
            hit.setdefault(c, []).append(key)
            print(f"  {name:4s} {key:18s} rep={r} orbit_size={orb[c]} all-high prefix 1..{prefix_high(r0, n)}")
        ok = len(hit) == len(orb) and all(len(v) == 1 for v in hit.values())
        print(f"  {name}: representatives hit every orbit exactly once: {ok}")
        assert ok
    # cross-match ours <-> wang by orbit
    if n == 9:
        for key in OURS:
            a = canon_of(tuple(tuple(x - 1 for x in e) for e in OURS[key]))
            b = canon_of(tuple(tuple(x - 1 for x in e) for e in WANG[key]))
            print(f"  same orbit ours/wang for {key}: {a == b}")
            assert a == b


if __name__ == "__main__":
    main()
