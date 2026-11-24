"""V-b spot check of Wang Cor 4.5 for the computed lists: every exactly-realizable 6-point GP set
with <= 6 distances must have its pattern (my canon6) in F u U (data/E1_F, data/E1_U).
Sources (exact, Q(sqrt3) via verify/checker_v.Q3Field):
  A. two concentric C3 orbits {a, Ra, R^2a, b, Rb, R^2b} (R = rotation by 120 deg), a,b on a grid
  B. all 6-subsets of my 8 radius-3 lattice n=7 crescents and of Palasti n=8 (verify/*.txt)
  C. two C3 orbits about DIFFERENT centres is not symmetric -> skipped; instead C3 orbit + mirror
     images: {a, Ra, R^2a} u {sigma a, ...} handled by A with b = mirror(a) included in grid.
"""
import sys, os, itertools, random, collections, json
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from common_vb import canon6, read_pattern_file
import checker_v as C
from fractions import Fraction as Fr

F3 = C.Q3Field


def rot(p):  # 120 degrees: (x,y) -> (-x/2 - s y/2, s x/2 - y/2), s = sqrt3
    x, y = p
    h = F3.c(Fr(-1, 2))
    s2 = F3.c(0, Fr(1, 2))
    return (F3.add(F3.mul(h, x), F3.mul(F3.c(0, Fr(-1, 2)), y)), F3.add(F3.mul(s2, x), F3.mul(h, y)))


def pattern_of(pts):
    """returns (labels over lex pairs, nlabels) or None if not GP / duplicate."""
    n = len(pts)
    ds = []
    for i, j in itertools.combinations(range(n), 2):
        d = C._d2(F3, pts[i], pts[j])
        if F3.is_zero(d):
            return None
        ds.append(d)
    for a, b, c in itertools.combinations(range(n), 3):
        if C.collinear(F3, pts[a], pts[b], pts[c]):
            return None
    for q in itertools.combinations(range(n), 4):
        if C.concyclic_circumcentre(F3, *[pts[i] for i in q]):
            return None
    reps = []
    lab = []
    for d in ds:
        for k, r in enumerate(reps):
            if F3.is_zero(F3.sub(d, r)):
                lab.append(k)
                break
        else:
            reps.append(d)
            lab.append(len(reps) - 1)
    return lab, len(reps)


def main():
    FU = {}
    for name, path in (('F', 'data/E1_F.txt.gz'), ('U', 'data/E1_U.txt.gz')):
        for p in read_pattern_file(path):
            FU[canon6(list(p))] = name
    print('F u U canonical', len(FU), collections.Counter(FU.values()), flush=True)
    seen = {}
    miss = []
    src_count = collections.Counter()

    def test(pts, src):
        r = pattern_of(pts)
        if r is None:
            return
        lab, nl = r
        if nl > 6:
            return
        cf = canon6(lab)
        src_count[src] += 1
        if cf not in seen:
            seen[cf] = src
            if cf not in FU:
                miss.append((src, ''.join(map(str, cf)), [[str(c) for c in p] for p in pts]))
    # B. lattice crescents
    for fn in ('verify/lattice_n7_r3.txt', 'verify/lattice_n8_r3.txt'):
        for line in open(fn):
            if 'lattice(i,j) =' in line:
                L = eval(line.split('=', 1)[1].split('sizes')[0])
                P = [C.lattice_to_q3(i, j) for i, j in L]
                for S in itertools.combinations(P, 6):
                    test(list(S), 'lattice:' + fn[-12:])
    print('after lattice:', len(seen), 'distinct patterns, misses', len(miss), flush=True)
    # A. two concentric C3 orbits on a grid of Q(sqrt3) points
    grid = []
    for x in range(-3, 4):
        for y in range(-3, 4):
            for t in (0, 1):
                # point (x + t*sqrt3/2 ... ) keep small: a = (x/2, y*sqrt3/2) or (x/2, y/2)
                grid.append((F3.c(Fr(x, 2)), F3.c(0, Fr(y, 2)) if t else F3.c(Fr(y, 2))))
    random.seed(217)
    for a, b in itertools.combinations(grid, 2):
        pts = [a, rot(a), rot(rot(a)), b, rot(b), rot(rot(b))]
        test(pts, 'C3')
    print('after C3:', len(seen), 'distinct patterns, misses', len(miss), dict(src_count), flush=True)
    # D. random 6-subsets of the radius-3 triangular patch with <= 6 distances (sampled)
    patch = [(i, j) for i in range(-3, 4) for j in range(-3, 4) if abs(i) <= 3 and abs(j) <= 3 and abs(i + j) <= 3]
    P = [C.lattice_to_q3(i, j) for i, j in patch]
    tried = 0
    while tried < 300000:
        S = random.sample(range(len(P)), 6)
        tried += 1
        # cheap prefilter on integer squared distances (lattice norm i^2+ij+j^2)
        ds = set()
        for u, v in itertools.combinations(S, 2):
            di, dj = patch[u][0] - patch[v][0], patch[u][1] - patch[v][1]
            ds.add(di * di + di * dj + dj * dj)
        if len(ds) <= 6:
            test([P[k] for k in S], 'tri_random')
    print('after random triangular:', len(seen), 'distinct patterns, misses', len(miss), dict(src_count), flush=True)
    json.dump(dict(distinct=len(seen), by_source=collections.Counter(seen.values()), raw=dict(src_count),
                   inF=sum(1 for c in seen if FU.get(c) == 'F'), inU=sum(1 for c in seen if FU.get(c) == 'U'),
                   misses=miss), open('verify/vb/cor45_vb.json', 'w'), indent=1)
    print('DONE distinct', len(seen), 'misses', len(miss))


if __name__ == '__main__':
    main()
