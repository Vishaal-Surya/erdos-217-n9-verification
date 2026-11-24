"""Orchestrator's third exact check of near-miss config 2 (core 111232433454125, q = u^2 - 17/5 u + 1).
Frame coords in K = Q(sqrt21), u = (17 - 3 sqrt21)/10 (the PD root ~0.32523). Metric G from reconstruction.
Counts: squared-distance class sizes, collinear triples (affine-invariant det), concyclic quadruples
(det rows (1, x, y, x^T G x) -- the G-metric lift, equivalent to Euclidean concyclicity under x_E = Lx, G = L^T L).
Also prints exact Euclidean coordinates via Cholesky."""
import itertools, json, sys
import sympy as sp
s21 = sp.sqrt(21)
u = (17 - 3*s21)/10
assert sp.expand(u**2 - sp.Rational(17,5)*u + 1) == 0
rec = None
for l in open('verify/vb/ext_full.jsonl'):
    r = json.loads(l)
    if r['key'].startswith('111232433454125'):
        for f in r['found']:
            if f['sizes'] == list(range(1,9)) and not f['collinear']:
                rec = f; break
    if rec: break
def el(cs): return sum(sp.Rational(c)*u**i for i, c in enumerate(cs))
pts = []
for p in rec['points']:
    i = p.index('|'); pts.append(sp.Matrix([el(p[:i]), el(p[i+1:])]))
G = sp.Matrix([[sp.Rational(2,3)-u/3, sp.Rational(1,2)-u/2], [sp.Rational(1,2)-u/2, sp.Rational(2,3)-u/3]])
Z = lambda e: sp.radsimp(sp.expand(e)) == 0
def d2(a, b): v = a-b; return sp.radsimp(sp.expand((v.T*G*v)[0]))
D = {}
for i, j in itertools.combinations(range(9), 2):
    D[(i,j)] = d2(pts[i], pts[j])
classes = {}
for k, v in D.items():
    for key in classes:
        if Z(key - v): classes[key].append(k); break
    else: classes[v] = [k]
sizes = sorted(len(v) for v in classes.values())
col = [t for t in itertools.combinations(range(9), 3)
       if Z((pts[t[1]]-pts[t[0]]).row_join(pts[t[2]]-pts[t[0]]).det())]
def lift(p): return [1, p[0], p[1], (p.T*G*p)[0]]
cyc = [q for q in itertools.combinations(range(9), 4) if Z(sp.Matrix([lift(pts[i]) for i in q]).det())]
print('class sizes', sizes); print('collinear triples', len(col), col); print('concyclic quadruples', len(cyc), cyc)
s = sp.sqrt(G[0,0]); t = G[0,1]/s; w = sp.sqrt(sp.radsimp(G[1,1]-t**2))
L = sp.Matrix([[s, t], [0, w]])
print('Euclidean coords (x_E = L x), numeric:')
for p in pts:
    e = L*p; print('  (%.12f, %.12f)' % (float(e[0]), float(e[1])))
print('squared distances by class (exact, size: value):')
for v, ks in sorted(classes.items(), key=lambda kv: len(kv[1])): print(' ', len(ks), sp.nsimplify(v), '~', float(v))
