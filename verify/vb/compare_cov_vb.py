"""compare V-b coverage survivors with builder's data/coverage_E1_found.txt and level7."""
import sys, os, collections, itertools
sys.path.insert(0, os.path.dirname(__file__))
from common_vb import canon_complete, rename, relabel_points, affine_reject, PI

mine = collections.defaultdict(set)
for l in open('verify/vb/cov_E1_found.txt'):
    n, s = l.split()
    mine[n].add(canon_complete([int(c) for c in s], 9))
theirs = collections.defaultdict(set)
for l in open('data/coverage_E1_found.txt'):
    n, s = l.split()
    theirs[n].add(canon_complete([int(c) for c in s], 9))
allm = set().union(*mine.values())
allt = set().union(*theirs.values())
print('mine total', len(allm), 'builder total', len(allt), 'equal sets:', allm == allt)
for n in sorted(set(mine) | set(theirs)):
    print(n, len(mine[n]), len(theirs[n]), mine[n] == theirs[n])
# level 7
perms7 = list(itertools.permutations(range(7)))
def c7(lab):
    return min(rename(relabel_points(lab, 7, p)) for p in perms7)
m7 = {c7([int(c) for c in l.split()[0]]) for l in open('verify/vb/cov_E1_level7.txt')}
t7 = {c7([int(c) for c in l.split()[0]]) for l in open('data/coverage_E1_level7.txt')}
print('level7 mine', len(m7), 'builder', len(t7), 'equal:', m7 == t7)
# affine on builder's 487 with my implementation, reasons
reasons = collections.Counter()
for cf in allt:
    r = affine_reject(list(cf[1]), 9)
    reasons[(r or 'PASS').split(':')[0]] += 1
print('my affine on builder 487:', dict(reasons))
