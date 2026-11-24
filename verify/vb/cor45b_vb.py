"""Cor 4.5 spot check, part B: all GP 6-subsets with <= 6 distances of the exact 9-point sets
(class sizes 1..8, not GP) found by extend_vb.py (number fields of degree up to 12, non-lattice).
Their patterns must lie in F u U."""
import sys, os, json, gzip, itertools, collections
sys.path.insert(0, os.path.dirname(__file__))
from extend_vb import KQ, ex_d2, ex_check_config
from common_vb import canon6, read_pattern_file

FU = {}
for name, path in (('F', 'data/E1_F.txt.gz'), ('U', 'data/E1_U.txt.gz')):
    for p in read_pattern_file(path):
        FU[canon6(list(p))] = name
recon = {}
with gzip.open('data/F_recon_final.jsonl.gz', 'rt') as f:
    for line in f:
        d = json.loads(line)
        for fi, fac in enumerate(d['factors']):
            if fac['status'] == 'REALIZED':
                recon['%s#%d' % (d['pattern'], fi)] = (fac['q'], fac['G'])
seen = {}
miss = []
nsub = 0
for l in open('verify/vb/ext_full.jsonl'):
    d = json.loads(l)
    for it in d.get('found') or []:
        if it.get('sizes') != list(range(1, 9)):
            continue
        q, Gs = recon[d['key']]
        K = KQ(q)
        G = [K.el(c) for c in Gs]
        pts = []
        for s in it['points']:
            k = s.index('|')
            pts.append((K.el(s[:k]), K.el(s[k + 1:])))
        col = set(map(tuple, it['collinear']))
        cyc = set(map(tuple, it['concyclic']))
        dist = {}
        for a, b in itertools.combinations(range(9), 2):
            dist[(a, b)] = ex_d2(K, G, pts[a], pts[b])
        for S in itertools.combinations(range(9), 6):
            if any(t in col for t in itertools.combinations(S, 3)) or any(t in cyc for t in itertools.combinations(S, 4)):
                continue
            reps, lab = [], []
            for a, b in itertools.combinations(S, 2):
                v = dist[(a, b)]
                for k2, r in enumerate(reps):
                    if K.iszero(r - v):
                        lab.append(k2)
                        break
                else:
                    reps.append(v)
                    lab.append(len(reps) - 1)
            if len(reps) > 6:
                continue
            nsub += 1
            cf = canon6(lab)
            if cf not in seen:
                seen[cf] = d['key']
                if cf not in FU:
                    miss.append((d['key'], S, ''.join(map(str, cf))))
res = dict(subsets_tested=nsub, distinct=len(seen), inF=sum(FU.get(c) == 'F' for c in seen),
           inU=sum(FU.get(c) == 'U' for c in seen), misses=miss)
print(res)
json.dump(res, open('verify/vb/cor45b_vb.json', 'w'), indent=1)
