"""Builder C: independent (Python, src/canon.py) re-check of the Java survivors against U.
For each complete 9-point survivor: class sizes 1..8, low pairs 1/2, local tests (canon.local_ok),
every all-high 6-subset canonical form in U. For each 7-point level pattern: every 6-subset in U.
Usage: python verify_survivors.py data/E1_U.txt.gz data/coverage_E1_found.txt data/coverage_E1_level7.txt
"""
import gzip, itertools, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import canon

uf, found, lev7 = sys.argv[1:4]
U = {canon.canonical(l.split()[0], 6) for l in gzip.open(uf, "rt") if l.strip()}
bad = 0
rows = [l.split() for l in open(found) if l.strip()]
for arr, s in rows:
    n = 9
    pairs = list(itertools.combinations(range(n), 2))
    assert sorted(s.count(str(i)) for i in range(1, 9)) == list(range(1, 9))
    if not canon.local_ok(s, n):
        bad += 1; continue
    lab = dict(zip(pairs, s))
    for T in itertools.combinations(range(n), 6):
        sub = [lab[p] for p in itertools.combinations(T, 2)]
        if "1" in sub or "2" in sub:
            continue
        if canon.canonical(sub, 6) not in U:
            bad += 1; break
print(f"survivors checked={len(rows)} inconsistent={bad}")
p7 = [l.strip() for l in open(lev7) if l.strip()]
b7 = 0
for s in p7:
    lab = dict(zip(itertools.combinations(range(7), 2), s))
    if not canon.local_ok(s, 7):
        b7 += 1; continue
    for T in itertools.combinations(range(7), 6):
        if canon.canonical([lab[p] for p in itertools.combinations(T, 2)], 6) not in U:
            b7 += 1; break
print(f"level7 patterns={len(p7)} inconsistent={b7}; canonical fixed points: {all(canon.canonical(s,7)==s for s in p7)}")
