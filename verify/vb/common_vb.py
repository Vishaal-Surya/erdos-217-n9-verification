"""V-b common utilities (independent of src/). Patterns are label lists over pairs in
lexicographic order (0,1),(0,2),...,(n-2,n-1) of 0-based points."""
import gzip, itertools
from flint import fmpz_mat


def pidx_table(n):
    t = {}
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            t[(i, j)] = k
            t[(j, i)] = k
            k += 1
    return t


PI = {n: pidx_table(n) for n in range(2, 10)}
PAIRS = {n: [(i, j) for i in range(n) for j in range(i + 1, n)] for n in range(2, 10)}


def rename(seq):
    """first-occurrence renaming of a label sequence (labels hashable)."""
    m = {}
    out = []
    for x in seq:
        if x not in m:
            m[x] = len(m)
        out.append(m[x])
    return tuple(out)


def relabel_points(lab, n, perm):
    """perm: new position of old point i is perm[i]. returns sequence over new lex pairs."""
    pi = PI[n]
    out = [None] * len(lab)
    for (i, j), k in zip(PAIRS[n], range(len(lab))):
        out[pi[(perm[i], perm[j])]] = lab[k]
    return out


PERMS6 = list(itertools.permutations(range(6)))


def orbit6(lab):
    """all renamed orderings of a 6-point pattern (as tuples)."""
    return {rename(relabel_points(lab, 6, p)) for p in PERMS6}


def canon6(lab):
    return min(orbit6(lab))


def read_pattern_file(path):
    op = gzip.open if path.endswith('.gz') else open
    out = []
    with op(path, 'rt') as f:
        for line in f:
            t = line.split()
            if t:
                out.append(tuple(int(c) for c in t[0]))
    return out


def sub_pattern(lab, n, pts):
    """labels of the sub-pattern on sorted point list pts (lex order)."""
    pi = PI[n]
    return [lab[pi[(pts[a], pts[b])]] for a in range(len(pts)) for b in range(a + 1, len(pts))]


# ---------------- local tests (Wang Lemma 2.1, 2.2) on a complete pattern ----------------
def local_ok_full(lab, n):
    pi = PI[n]
    c = lambda a, b: lab[pi[(a, b)]]
    for a in range(n):
        cnt = {}
        for x in range(n):
            if x != a:
                l = c(a, x)
                cnt[l] = cnt.get(l, 0) + 1
                if cnt[l] > 3:
                    return False, 'L21i'
    for a in range(n):
        for b in range(a + 1, n):
            cnt = {}
            eq = 0
            for x in range(n):
                if x == a or x == b:
                    continue
                k = (c(a, x), c(b, x))
                cnt[k] = cnt.get(k, 0) + 1
                if cnt[k] > 2:
                    return False, 'L21ii'
                if k[0] == k[1]:
                    eq += 1
                    if eq > 2:
                        return False, 'L21iii'
    for q in itertools.combinations(range(n), 4):
        a, b, cc, d = q
        if c(a, b) == c(cc, d) and c(a, cc) == c(b, d) and c(a, d) == c(b, cc):
            return False, 'L22'
    return True, ''


# ---------------- affine tests (Wang Lemma 2.3, parallelogram rows only) ----------------
def affine_reject(lab, n, centroid=False):
    """returns None if passes, else reason string. Exact integer ranks via flint."""
    pi = PI[n]
    c = lambda a, b: lab[pi[(a, b)]]
    rows = set()
    for a in range(n):
        for b in range(n):
            if b == a:
                continue
            for cc in range(n):
                if cc in (a, b):
                    continue
                if c(a, b) is None:
                    continue
                for d in range(n):
                    if d in (a, b, cc):
                        continue
                    if c(a, b) == c(cc, d) and c(a, cc) == c(b, d):
                        r = [0] * n
                        r[a] += 1; r[d] += 1; r[b] -= 1; r[cc] -= 1
                        rows.add(tuple(r))
    if centroid:
        for a, b, cc in itertools.combinations(range(n), 3):
            if c(a, b) == c(a, cc) == c(b, cc):
                for d in range(n):
                    if d not in (a, b, cc) and c(a, d) == c(b, d) == c(cc, d):
                        r = [0] * n
                        r[a] += 1; r[b] += 1; r[cc] += 1; r[d] -= 3
                        rows.add(tuple(r))
    if not rows:
        return None
    rows = sorted(rows)
    R = fmpz_mat([list(r) for r in rows])
    rk = R.rank()
    for I in itertools.combinations(range(n), 3):
        cols = [j for j in range(n) if j not in I]
        Rs = fmpz_mat([[r[j] for j in cols] for r in rows])
        if Rs.rank() < rk:
            return 'T1:%s' % (I,)
    labels = sorted(set(lab))
    li = {l: k for k, l in enumerate(labels)}
    L = len(labels)
    Arows = set()
    for r in rows:
        for j in range(1, n):
            v = [0] * L
            for i in range(n):
                if r[i] == 0:
                    continue
                if i != j:
                    v[li[c(i, j)]] += r[i]
                if i != 0:
                    v[li[c(i, 0)]] -= r[i]
            if any(v):
                Arows.add(tuple(v))
    if not Arows:
        return None
    Arows = sorted(Arows)
    A = fmpz_mat([list(r) for r in Arows])
    ra = A.rank()
    for l in range(L):
        e = [0] * L; e[l] = 1
        if fmpz_mat([list(r) for r in Arows] + [e]).rank() == ra:
            return 'T2:z%d=0' % labels[l]
    for l in range(L):
        for m in range(l + 1, L):
            e = [0] * L; e[l] = 1; e[m] = -1
            if fmpz_mat([list(r) for r in Arows] + [e]).rank() == ra:
                return 'T2:z%d=z%d' % (labels[l], labels[m])
    return None


# ---------------- canonical form of complete patterns (labels are canonical already) -------
def canon_complete(lab, n):
    """labels must be intrinsic (e.g. class sizes). min over permutations, with an invariant
    refinement: points are ordered by (sorted labels at point); only perms within cells."""
    pi = PI[n]
    inv = []
    for a in range(n):
        inv.append(tuple(sorted(lab[pi[(a, x)]] for x in range(n) if x != a)))
    order = sorted(range(n), key=lambda a: inv[a])
    cells = []
    for a in order:
        if cells and inv[cells[-1][0]] == inv[a]:
            cells[-1].append(a)
        else:
            cells.append([a])
    best = None
    for choice in itertools.product(*[list(itertools.permutations(cl)) for cl in cells]):
        seq = [x for cl in choice for x in cl]  # seq[k] = old point at new position k
        s = tuple(lab[pi[(seq[i], seq[j])]] for i in range(n) for j in range(i + 1, n))
        if best is None or s < best:
            best = s
    return (tuple(inv[a] for a in order), best)
