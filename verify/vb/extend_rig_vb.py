"""V-b RIGOROUS extension check (Task 3, follow-up): same case split as extend_vb.py, but every
decision is certified: a value is declared nonzero / two values unequal only if arb balls prove
it (start 128 bits, escalate to 256, 512); everything else is decided by exact arithmetic in
K_q = Q[u]/(q) (affine frame with metric G). No float64 anywhere in decisions.

Case split (proved in NOTES/verification.md, V-b item 3):
 m=0: every outside point = trilateration point of a core triangle with values in D(S)^3;
      all 3 outside points from that set; exactly 3 pairs carry values outside D(S).
 m=1: two outside points y,z are trilateration points (D(S)^3) or circumcentres of core
      triangles; the third w is a trilateration point of a core triangle with values in D(Q)^3
      or a circumcentre of a triangle of Q = S+{y,z}.
Usage: python extend_rig_vb.py in.jsonl.gz out_prefix
"""
import sys, os, json, gzip, itertools, time, collections
sys.path.insert(0, os.path.dirname(__file__))
from flint import fmpq_poly, arb, ctx
from extend_vb import KQ, ex_d2, ex_trilat, ex_circum, ex_check_config, ev

PRECS = (128, 256, 512)
LIMIT = None
CNT = collections.Counter()


def key(p):
    return str(p)


class Emb:
    """one real PD embedding; arb data at several precisions, computed lazily."""

    def __init__(self, K, G, cpts, z, root_index, nroots_check):
        self.K, self.G, self.cpts, self.z = K, G, cpts, z
        self.i = root_index
        self.cache = {}

    def theta(self, prec):
        if ('th', prec) not in self.cache:
            old = ctx.prec
            ctx.prec = prec
            rts = real_roots_prec(self.K.q)
            ctx.prec = old
            if prec == PRECS[0]:
                self.cache[('th', prec)] = rts[self.i]
            else:
                base = self.theta(PRECS[0])
                m = [r for r in rts if r.overlaps(base)]
                assert len(m) == 1
                self.cache[('th', prec)] = m[0]
        return self.cache[('th', prec)]

    def val(self, el, prec):
        old = ctx.prec
        ctx.prec = prec
        r = ev(el, self.theta(prec))
        ctx.prec = old
        return r

    def pt(self, P, prec):
        return (self.val(P[0], prec), self.val(P[1], prec))

    def Gb(self, prec):
        if ('G', prec) not in self.cache:
            self.cache[('G', prec)] = [self.val(g, prec) for g in self.G]
        return self.cache[('G', prec)]


def real_roots_prec(q):
    rts = [r for r, m in q.complex_roots()]
    out = []
    for r in rts:
        if r.imag.contains(0):
            if sum(1 for s in rts if s.real.overlaps(r.real) and s.imag.overlaps(-r.imag)) == 1:
                out.append(r.real)
    return out


def d2b(G, p, r):
    vx, vy = p[0] - r[0], p[1] - r[1]
    return G[0] * vx * vx + 2 * G[1] * vx * vy + G[2] * vy * vy


def gdot(G, p, r):
    return G[0] * p[0] * r[0] + G[1] * (p[0] * r[1] + p[1] * r[0]) + G[2] * p[1] * r[1]


def tri_coeffs(G, a, b, c):
    """x = x0 + u w1 + v w2 ; |x-a|^2 = Q0 + B1 u + B2 v + Q1 u^2 + B12 u v + Q2 v^2."""
    r1 = (G[0] * (b[0] - a[0]) + G[1] * (b[1] - a[1]), G[1] * (b[0] - a[0]) + G[2] * (b[1] - a[1]))
    r2 = (G[0] * (c[0] - a[0]) + G[1] * (c[1] - a[1]), G[1] * (c[0] - a[0]) + G[2] * (c[1] - a[1]))
    det = 4 * (r1[0] * r2[1] - r1[1] * r2[0])
    if det.contains(0):
        return None

    def solve(h1, h2):
        return ((h1 * 2 * r2[1] - h2 * 2 * r1[1]) / det, (2 * r1[0] * h2 - 2 * r2[0] * h1) / det)
    na = gdot(G, a, a)
    x0 = solve(gdot(G, b, b) - na, gdot(G, c, c) - na)
    w1 = solve(arb(1), arb(0))
    w2 = solve(arb(0), arb(1))
    d = (x0[0] - a[0], x0[1] - a[1])
    return dict(x0=x0, w1=w1, w2=w2, Q0=gdot(G, d, d), B1=2 * gdot(G, d, w1), B2=2 * gdot(G, d, w2),
                Q1=gdot(G, w1, w1), B12=2 * gdot(G, w1, w2), Q2=gdot(G, w2, w2))


def trilat_scan(E, tri_pts_exact, vals_exact, vals_keys, explained, on_exact_point):
    """scan all (al,be,ga) in vals^3 for the triangle; certified rejection by balls (escalating),
    otherwise exact. tri_pts_exact: 3 exact points. Calls on_exact_point(P) for every exact
    solution. explained: set of (tri_id, i, j, k) already known to be solutions (skips exact)."""
    prec = PRECS[0]
    ctx.prec = prec
    G = E.Gb(prec)
    a, b, c = (E.pt(p, prec) for p in tri_pts_exact)
    co = tri_coeffs(G, a, b, c)
    assert co is not None
    V = [E.val(v, prec) for v in vals_exact]
    n = len(V)
    Q0, B1, B2, Q1, B12, Q2 = co['Q0'], co['B1'], co['B2'], co['Q1'], co['B12'], co['Q2']
    for i in range(n):
        al = V[i]
        base = Q0 - al
        Up = [(al - V[j], B1 * (al - V[j]) + Q1 * (al - V[j]) * (al - V[j])) for j in range(n)]
        Vp = [(al - V[k], B2 * (al - V[k]) + Q2 * (al - V[k]) * (al - V[k])) for k in range(n)]
        for j in range(n):
            u, uu = Up[j]
            bu = base + uu
            B12u = B12 * u
            for k in range(n):
                v, vv = Vp[k]
                f = bu + vv + B12u * v
                CNT['trilat_trials'] += 1
                if not f.contains(0):
                    CNT['trilat_reject_ball128'] += 1
                    continue
                if (i, j, k) in explained:
                    CNT['trilat_explained_exact'] += 1
                    continue
                # escalate precision for this single triple
                done = False
                for p2 in PRECS[1:]:
                    ctx.prec = p2
                    G2 = E.Gb(p2)
                    a2, b2, c2 = (E.pt(p, p2) for p in tri_pts_exact)
                    co2 = tri_coeffs(G2, a2, b2, c2)
                    A, Bv, Cv = (E.val(vals_exact[t], p2) for t in (i, j, k))
                    u2, v2 = A - Bv, A - Cv
                    f2 = co2['Q0'] - A + co2['B1'] * u2 + co2['B2'] * v2 + co2['Q1'] * u2 * u2 + co2['B12'] * u2 * v2 + co2['Q2'] * v2 * v2
                    if not f2.contains(0):
                        CNT['trilat_reject_ball%d' % p2] += 1
                        done = True
                        break
                ctx.prec = prec
                if done:
                    continue
                CNT['trilat_exact_tests'] += 1
                P = ex_trilat(E.K, E.G, *tri_pts_exact, vals_exact[i], vals_exact[j], vals_exact[k])
                if P is None:
                    CNT['trilat_reject_exact'] += 1
                    continue
                on_exact_point(P)


def components(balls):
    """number of connected components of the overlap graph (lower bound on #distinct values)."""
    n = len(balls)
    parent = list(range(n))

    def f(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i in range(n):
        for j in range(i + 1, n):
            if balls[i].overlaps(balls[j]):
                parent[f(i)] = f(j)
    return len({f(i) for i in range(n)})


def final_check(E, pts9, stats_prefix):
    """certified rejection if >8 ball-components among 36 distance balls (each component is a
    set of values pairwise-distinct from other components). else exact full check."""
    ctx.prec = PRECS[0]
    G = E.Gb(PRECS[0])
    B = [E.pt(p, PRECS[0]) for p in pts9]
    balls = [d2b(G, B[i], B[j]) for i, j in itertools.combinations(range(9), 2)]
    if components(balls) > 8:
        CNT[stats_prefix + '_reject_ball'] += 1
        return None
    CNT[stats_prefix + '_exact'] += 1
    ex = ex_check_config(E.K, E.G, pts9)
    if ex.get('dup'):
        return dict(sizes=None)
    return ex


def search_embedding(E, m, found):
    K = E.K
    core = E.cpts
    z = E.z
    zkeys = [key(v) for v in z]
    zidx = {k: i for i, k in enumerate(zkeys)}
    corekeys = {(key(p[0]), key(p[1])) for p in core}
    cands = {}      # key -> dict(P, cls: list of z-index or -1, ndvals: keys of non-D values)

    def add_cand(P):
        kk = (key(P[0]), key(P[1]))
        if kk in corekeys or kk in cands:
            return kk
        ds = [ex_d2(K, E.G, P, c) for c in core]
        dk = [key(d) for d in ds]
        cls = [zidx.get(k, -1) for k in dk]
        cands[kk] = dict(P=P, cls=cls, dk=dk, ds=ds)
        CNT['cand_distinct'] += 1
        return kk

    tris = list(itertools.combinations(range(6), 3))
    for t in tris:
        explained = set()
        # triples explained by already-known candidates (exact: P satisfies the 3 equations)
        for c in cands.values():
            ci = [c['cls'][x] for x in t]
            if min(ci) >= 0:
                explained.add(tuple(ci))
        trilat_scan(E, [core[x] for x in t], z, zkeys, explained, add_cand)
        if m == 1:
            CNT['circum_core'] += 1
            add_cand(ex_circum(K, E.G, *[core[x] for x in t]))
    C = list(cands.values())
    nd = [sum(1 for x in c['cls'] if x < 0) for c in C]
    prec = PRECS[0]
    ctx.prec = prec
    Gb = E.Gb(prec)
    Zb = [E.val(v, prec) for v in z]
    Cb = [E.pt(c['P'], prec) for c in C]

    def pair_dist_class(i, j):
        """index in z of d(Ci,Cj) or -1; certified (balls, else exact). also exact value lazily."""
        ctx.prec = prec
        db = d2b(Gb, Cb[i], Cb[j])
        ov = [k for k, zb in enumerate(Zb) if db.overlaps(zb)]
        if not ov:
            CNT['pairdist_ball'] += 1
            return -1
        CNT['pairdist_exact'] += 1
        dk = key(ex_d2(K, E.G, C[i]['P'], C[j]['P']))
        return zidx.get(dk, -1)

    if m == 0:
        good = [i for i in range(len(C)) if nd[i] <= 3]
        pe = {}
        for i, j in itertools.combinations(good, 2):
            if nd[i] + nd[j] > 3:
                continue
            e = 1 if pair_dist_class(i, j) < 0 else 0
            if nd[i] + nd[j] + e <= 3:
                pe[(i, j)] = e
        CNT['pairs_ok'] += len(pe)
        for (i, j), e in list(pe.items()):
            for k in good:
                if k <= j or (i, k) not in pe or (j, k) not in pe:
                    continue
                if nd[i] + nd[j] + nd[k] + e + pe[(i, k)] + pe[(j, k)] != 3:
                    continue
                CNT['triples_budget'] += 1
                pts = list(core) + [C[i]['P'], C[j]['P'], C[k]['P']]
                ex = final_check(E, pts, 'final')
                if ex is not None and ex.get('sizes') == list(range(1, 9)):
                    found.append(ex)
    else:
        corecnt = collections.Counter()
        pi = [(a, b) for a in range(6) for b in range(a + 1, 6)]
        for (a, b) in pi:
            corecnt[zidx[key(ex_d2(K, E.G, core[a], core[b]))]] += 1
        for i, j in itertools.combinations(range(len(C)), 2):
            # exact multiset of values of Q = S + {Ci, Cj}
            cnt = collections.Counter({('z', k): v for k, v in corecnt.items()})
            vals = {}
            for c in (C[i], C[j]):
                for x, dk, d in zip(c['cls'], c['dk'], c['ds']):
                    kk = ('z', x) if x >= 0 else ('n', dk)
                    cnt[kk] += 1
                    vals[kk] = d
            dij = ex_d2(K, E.G, C[i]['P'], C[j]['P'])
            CNT['pairdist_exact'] += 1
            dk = key(dij)
            kk = ('z', zidx[dk]) if dk in zidx else ('n', dk)
            cnt[kk] += 1
            vals[kk] = dij
            corek = [k for k in cnt if k[0] == 'z']
            other = [k for k in cnt if k[0] == 'n']
            if len(other) > 3 or not capacity_ok([cnt[k] for k in corek], [cnt[k] for k in other]):
                continue
            CNT['pairs_ok'] += 1
            DQ = [z[k[1]] for k in corek] + [vals[k] for k in other]
            DQk = [key(v) for v in DQ]
            Q = list(core) + [C[i]['P'], C[j]['P']]
            ws = {}
            wcls = []
            DQidx = {k2: n2 for n2, k2 in enumerate(DQk)}

            def add_w(P):
                kk2 = (key(P[0]), key(P[1]))
                if kk2 not in ws:
                    ws[kk2] = P
                    wcls.append([DQidx.get(key(ex_d2(K, E.G, P, c)), -1) for c in core])
            for t in tris:
                expl = set()
                for cl in wcls:
                    ci = tuple(cl[x] for x in t)
                    if min(ci) >= 0:
                        expl.add(ci)
                trilat_scan(E, [core[x] for x in t], DQ, DQk, expl, add_w)
            for t in itertools.combinations(range(8), 3):
                ctx.prec = prec
                qa, qb, qc = (E.pt(Q[x], prec) for x in t)
                if tri_coeffs(Gb, qa, qb, qc) is None:
                    # possibly collinear: decide exactly
                    dt = K.mul(Q[t[1]][0] - Q[t[0]][0], Q[t[2]][1] - Q[t[0]][1]) - K.mul(Q[t[1]][1] - Q[t[0]][1], Q[t[2]][0] - Q[t[0]][0])
                    if K.iszero(dt):
                        CNT['circum_collinear_exact'] += 1
                        continue
                CNT['circum_Q'] += 1
                add_w(ex_circum(K, E.G, *[Q[x] for x in t]))
            for k2, c in enumerate(C):
                if k2 != i and k2 != j:
                    add_w(c['P'])
            qkeys = {(key(p[0]), key(p[1])) for p in Q}
            for kk2, P in ws.items():
                if kk2 in qkeys:
                    continue
                CNT['last_cands'] += 1
                ex = final_check(E, Q + [P], 'final')
                if ex is not None and ex.get('sizes') == list(range(1, 9)):
                    found.append(ex)


def capacity_ok(counts_core, counts_other):
    items = [(c, True) for c in counts_core] + [(c, False) for c in counts_other]
    if len(items) > 8:
        return False
    match_of = {}

    def try_assign(i, seen):
        c, core = items[i]
        for s in range(1, 9):
            if s < c or (core and s < 3) or s in seen:
                continue
            seen.add(s)
            if s not in match_of or try_assign(match_of[s], seen):
                match_of[s] = i
                return True
        return False
    return all(try_assign(i, set()) for i in range(len(items)))


def main():
    global LIMIT
    inp, out = sys.argv[1], sys.argv[2]
    for a in sys.argv:
        if a.startswith('--limit='):
            LIMIT = int(a.split('=')[1])
    done = set()
    if os.path.exists(out + '.jsonl'):
        for l in open(out + '.jsonl'):
            done.add(json.loads(l)['key'])
    fo = open(out + '.jsonl', 'a')
    t0 = time.time()
    with gzip.open(inp, 'rt') as f:
        for ln, line in enumerate(f):
            if LIMIT is not None and ln >= LIMIT:
                break
            d = json.loads(line)
            pat = [int(c) for c in d['pattern']]
            for fi, fac in enumerate(d['factors']):
                kkey = '%s#%d' % (d['pattern'], fi)
                if kkey in done or fac['status'] != 'REALIZED':
                    continue
                K = KQ(fac['q'])
                G = [K.el(c) for c in fac['G']]
                cpts = [(K.el(p[0]), K.el(p[1])) for p in fac['points']]
                z = [K.el(c) for c in fac['z']]
                ctx.prec = PRECS[0]
                rts = real_roots_prec(K.q)
                rec = dict(key=kkey, m=6 - d['r'], emb=[], found=[])
                for ri, th in enumerate(rts):
                    g11, g12, g22 = (ev(x, th) for x in G)
                    if not (g11 > 0 and g11 * g22 - g12 * g12 > 0):
                        assert g11 < 0 or g11 * g22 - g12 * g12 < 0, 'PD undecided'
                        continue
                    E = Emb(K, G, cpts, z, ri, len(rts))
                    CNT['embeddings'] += 1
                    found = []
                    search_embedding(E, 6 - d['r'], found)
                    rec['emb'].append(ri)
                    for ex in found:
                        sz = ex.get('sizes')
                        it = dict(root=ri, sizes=sz, ncol=len(ex.get('collinear') or []),
                                  ncyc=len(ex.get('concyclic') or []))
                        if sz == list(range(1, 9)):
                            CNT['exact_sizes_ok'] += 1
                            if not ex['collinear'] and not ex['concyclic']:
                                CNT['CRESCENT'] += 1
                                it['CRESCENT'] = True
                        rec['found'].append(it)
                fo.write(json.dumps(rec) + '\n')
                fo.flush()
            if ln % 100 == 0:
                print(ln, '%.0fs' % (time.time() - t0), dict(CNT), flush=True)
    print('DONE', '%.0fs' % (time.time() - t0), dict(CNT), flush=True)


if __name__ == '__main__':
    main()
