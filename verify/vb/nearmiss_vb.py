"""Exact re-check (own K_q code) of the builder's near-miss 9-point sets (F_extend near_non_gp)."""
import sys, os, json, gzip, collections
sys.path.insert(0, os.path.dirname(__file__))
from extend_vb import KQ, ex_check_config, real_roots, ev, to_f

out = []
with gzip.open('data/F_extend.jsonl.gz', 'rt') as f:
    for line in f:
        d = json.loads(line)
        if not d.get('near_non_gp'):
            continue
        K = KQ(d['q'])
        G = [K.el(c) for c in d['G']]
        for cfg in d['near_non_gp']:
            pts = [(K.el(p[0]), K.el(p[1])) for p in cfg]
            ex = ex_check_config(K, G, pts)
            sizes = ex['sizes']
            # numeric normalised delta_k (class of size k), at each PD real root
            norm = []
            for th in real_roots(K.q):
                g11, g12, g22 = (ev(x, th) for x in G)
                if not (g11 > 0 and g11 * g22 - g12 * g12 > 0):
                    continue
                vals = {len(cl): to_f(ev(v, th)) for v, cl in zip(ex['values'], ex['classes'])}
                s = vals[8]
                norm.append([round(vals[k] / s, 4) for k in range(1, 9)])
            rec = dict(key=d['key'], q=d['q'], sizes=sizes, n_collinear=len(ex['collinear']),
                       concyclic=ex['concyclic'], delta_over_delta8=norm)
            out.append(rec)
            print(json.dumps(rec))
json.dump(out, open('verify/vb/nearmiss_vb.json', 'w'), indent=1)
