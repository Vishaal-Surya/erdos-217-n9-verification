"""Builder E: independent Python canonical form + local tests for k-point distance patterns.

Canonical form: minimum over ALL k! point permutations of the first-occurrence-renamed label
string over lex pairs (1,2),(1,3),...,(k-1,k); labels written as digits '1'..'6'.
Independent re-implementation of src/enum/PatternEnum.java (used for cross-checks).
"""
import itertools
from functools import lru_cache


@lru_cache(None)
def _perm_maps(k):
    idx = {}
    p = 0
    for i in range(k):
        for j in range(i + 1, k):
            idx[(i, j)] = idx[(j, i)] = p
            p += 1
    maps = []
    for pi in itertools.permutations(range(k)):
        maps.append(tuple(idx[(pi[i], pi[j])] for i in range(k) for j in range(i + 1, k)))
    return tuple(maps)


def k_from_len(n):
    k = 2
    while k * (k - 1) // 2 < n:
        k += 1
    assert k * (k - 1) // 2 == n
    return k


def relabel(seq):
    ren = {}
    out = []
    for x in seq:
        if x not in ren:
            ren[x] = len(ren)
        out.append(ren[x])
    return tuple(out)


def canonical(labels, k=None):
    """labels: sequence of hashable labels over lex pairs. Returns canonical digit string."""
    labels = list(labels)
    if k is None:
        k = k_from_len(len(labels))
    best = min(relabel([labels[m] for m in mp]) for mp in _perm_maps(k))
    return "".join(str(x + 1) for x in best)


def matrix(labels, k):
    c = [[None] * k for _ in range(k)]
    p = 0
    for i in range(k):
        for j in range(i + 1, k):
            c[i][j] = c[j][i] = labels[p]
            p += 1
    return c


def local_ok(labels, k):
    """Wang Lemma 2.1 (i),(ii),(iii) and Lemma 2.2 on a complete labelling."""
    c = matrix(labels, k)
    for a in range(k):
        vals = [c[a][x] for x in range(k) if x != a]
        if any(vals.count(v) > 3 for v in set(vals)):
            return False
    for a in range(k):
        for b in range(a + 1, k):
            xs = [x for x in range(k) if x not in (a, b)]
            keys = [(c[a][x], c[b][x]) for x in xs]
            if any(keys.count(q) > 2 for q in set(keys)):
                return False
            if sum(1 for u, v in keys if u == v) > 2:
                return False
    for a, b, cc, d in itertools.combinations(range(k), 4):
        if c[a][b] == c[cc][d] and c[a][cc] == c[b][d] and c[a][d] == c[b][cc]:
            return False
    return True
