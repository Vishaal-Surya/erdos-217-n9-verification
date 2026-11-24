"""V-b independent extension check (Task 3).

Input: data/F_recon_final.jsonl.gz (builder's reconstructed F cores; only q, z, G, points are used,
and they are re-verified exactly here). Design (differs from NOTES/extension.md):
 * per real PD embedding (roots of q computed here with arb), not one abstract-field search;
 * candidates for outside points come from my own case split (no 'high-linked' Lemma H, no covering
   family): m=0: EVERY outside point is a trilateration point of SOME triangle of S with values in
   D(S)^3 (it has >=3 high pairs to S, H=D(S)); all three outside points are taken from this set,
   with the exact budget 'exactly 3 pairs outside D(S)'.
   m=1: >=2 outside points have <=1 low pair to S, hence >=3 partners with D(S) values
   (trilateration) or >=3 partners at the missing value d0 (circumcentre); the third point w is a
   trilateration point (all 56 triangles of Q, values D(Q)^3) or a circumcentre of Q=S+{y,z}.
 * screening in float64 Euclidean coordinates with tolerance tracking; every survivor of the
   class-size test is re-checked EXACTLY in K_q = Q[u]/(q) (affine frame with metric G).
Usage: python extend_vb.py in.jsonl.gz out_prefix [--limit N] [--stride k --offset j]
"""
import sys, os, json, gzip, itertools, time, collections
import numpy as np
from fractions import Fraction
from flint import fmpq_poly, fmpq, arb, ctx

ctx.prec = 200


def _fq(c):
    f = Fraction(c)
    return fmpq(f.numerator, f.denominator)
TOL = 1e-8
STATS = collections.Counter()
MARG = {'acc_max': 0.0, 'rej_min': 1e9}


# ------------------------------- exact K_q --------------------------------------
class KQ:
    def __init__(self, qcoeffs):
        self.q = fmpq_poly([_fq(c) for c in qcoeffs])

    def el(self, coeffs):
        return fmpq_poly([_fq(c) for c in coeffs]) % self.q if coeffs else fmpq_poly([])

    def mul(self, a, b):
        return (a * b) % self.q

    def inv(self, a):
        g, s, t = a.xgcd(self.q)
        assert g.degree() == 0 and not g.is_zero(), 'not invertible'
        return (s / g.coeffs()[0]) % self.q

    def iszero(self, a):
        return (a % self.q).is_zero()


def ex_d2(K, G, p, r):
    vx, vy = p[0] - r[0], p[1] - r[1]
    return (K.mul(G[0], K.mul(vx, vx)) + 2 * K.mul(G[1], K.mul(vx, vy)) + K.mul(G[2], K.mul(vy, vy))) % K.q


def ex_norm(K, G, p):
    return ex_d2(K, G, p, (fmpq_poly([]), fmpq_poly([])))


def ex_Gv(K, G, v):
    return (K.mul(G[0], v[0]) + K.mul(G[1], v[1]), K.mul(G[1], v[0]) + K.mul(G[2], v[1]))


def ex_trilat(K, G, a, b, c, al, be, ga):
    """x with |x-a|^2=al,|x-b|^2=be,|x-c|^2=ga (if consistent) else None. affine frame."""
    r1 = ex_Gv(K, G, (b[0] - a[0], b[1] - a[1]))
    r2 = ex_Gv(K, G, (c[0] - a[0], c[1] - a[1]))
    na = ex_norm(K, G, a)
    h1 = (ex_norm(K, G, b) - na + al - be) % K.q
    h2 = (ex_norm(K, G, c) - na + al - ga) % K.q
    # 2 r1.x = h1, 2 r2.x = h2
    det = (K.mul(r1[0], r2[1]) - K.mul(r1[1], r2[0])) % K.q
    di = K.inv(2 * det)
    x = K.mul(K.mul(h1, r2[1]) - K.mul(h2, r1[1]), di)
    y = K.mul(K.mul(r1[0], h2) - K.mul(r2[0], h1), di)
    P = (x, y)
    if not K.iszero(ex_d2(K, G, P, a) - al):
        return None
    return P


def ex_circum(K, G, a, b, c):
    r1 = ex_Gv(K, G, (b[0] - a[0], b[1] - a[1]))
    r2 = ex_Gv(K, G, (c[0] - a[0], c[1] - a[1]))
    na = ex_norm(K, G, a)
    h1 = (ex_norm(K, G, b) - na) % K.q
    h2 = (ex_norm(K, G, c) - na) % K.q
    det = (K.mul(r1[0], r2[1]) - K.mul(r1[1], r2[0])) % K.q
    di = K.inv(2 * det)
    x = K.mul(K.mul(h1, r2[1]) - K.mul(h2, r1[1]), di)
    y = K.mul(K.mul(r1[0], h2) - K.mul(r2[0], h1), di)
    return (x, y)


def ex_check_config(K, G, pts):
    """exact: class sizes, collinear triples, concyclic quadruples."""
    n = len(pts)
    d = {}
    for i, j in itertools.combinations(range(n), 2):
        d[(i, j)] = ex_d2(K, G, pts[i], pts[j])
    classes = []
    for k, v in d.items():
        if K.iszero(v):
            return dict(dup=True)
        for cl in classes:
            if K.iszero(cl[0] - v):
                cl[1].append(k)
                break
        else:
            classes.append([v, [k]])
    sizes = sorted(len(c[1]) for c in classes)
    col = []
    for a, b, c in itertools.combinations(range(n), 3):
        dt = K.mul(pts[b][0] - pts[a][0], pts[c][1] - pts[a][1]) - K.mul(pts[b][1] - pts[a][1], pts[c][0] - pts[a][0])
        if K.iszero(dt):
            col.append((a, b, c))
    cyc = []
    nrm = [ex_norm(K, G, p) for p in pts]
    for quad in itertools.combinations(range(n), 4):
        M = [[fmpq_poly([1]), pts[i][0], pts[i][1], nrm[i]] for i in quad]
        if K.iszero(det4(K, M)):
            cyc.append(quad)
    return dict(sizes=sizes, collinear=col, concyclic=cyc,
                values=[c[0] for c in classes], classes=[c[1] for c in classes])


def det4(K, M):
    # Laplace / Bareiss-free cofactor expansion (exact, small)
    def det3(m):
        return (K.mul(m[0][0], K.mul(m[1][1], m[2][2]) - K.mul(m[1][2], m[2][1]))
                - K.mul(m[0][1], K.mul(m[1][0], m[2][2]) - K.mul(m[1][2], m[2][0]))
                + K.mul(m[0][2], K.mul(m[1][0], m[2][1]) - K.mul(m[1][1], m[2][0]))) % K.q
    s = fmpq_poly([])
    for j in range(4):
        minor = [[M[i][k] for k in range(4) if k != j] for i in range(1, 4)]
        term = K.mul(M[0][j], det3(minor))
        s = s + term if j % 2 == 0 else s - term
    return s % K.q


# ------------------------------- numerics ---------------------------------------
def ev(poly, th):
    c = poly.coeffs()
    r = arb(0)
    for a in reversed(c):
        r = r * th + arb(a.p) / arb(a.q)
    return r


def to_f(a):
    return float(a.mid())


def real_roots(q):
    rts = [r for r, m in q.complex_roots()]
    out = []
    for r in rts:
        if r.imag.contains(0):
            # real iff its conjugate ball overlaps no OTHER root ball
            conj_overlap = sum(1 for s in rts if s.real.overlaps(r.real) and s.imag.overlaps(-r.imag))
            if conj_overlap == 1:
                out.append(r.real)
    return out


def trilat_vec(a, b, c, A, B, Cc):
    """vectorised Euclidean trilateration; A,B,C arrays of squared distances."""
    M = 2 * np.array([b - a, c - a])
    if abs(np.linalg.det(M)) <= 1e-12 * (np.abs(M).max() ** 2 + 1e-300):
        STATS['degenerate_triangle_skipped'] += 1
        return np.zeros((len(A), 2)), np.full(len(A), np.inf)
    Minv = np.linalg.inv(M)
    h1 = b @ b - a @ a + A - B
    h2 = c @ c - a @ a + A - Cc
    X = (Minv @ np.vstack([h1, h2])).T
    res = ((X - a) ** 2).sum(1) - A
    return X, res


def circum(a, b, c):
    M = 2 * np.array([b - a, c - a])
    h = np.array([b @ b - a @ a, c @ c - a @ a])
    if abs(np.linalg.det(M)) <= 1e-12 * (np.abs(M).max() ** 2 + 1e-300):
        STATS['degenerate_triangle_skipped'] += 1
        return np.array([np.nan, np.nan])
    return np.linalg.solve(M, h)


def cluster_vals(vals, scale):
    """group values numerically. returns representative list and index per value."""
    order = np.argsort(vals)
    reps = []
    idx = np.empty(len(vals), int)
    for k in order:
        v = vals[k]
        if reps and abs(v - reps[-1][0]) <= TOL * scale:
            idx[k] = len(reps) - 1
            note_acc(abs(v - reps[-1][0]) / scale)
        else:
            if reps:
                note_rej(abs(v - reps[-1][0]) / scale)
            reps.append([v])
            idx[k] = len(reps) - 1
    return [r[0] for r in reps], idx


def note_acc(x):
    if x > MARG['acc_max']:
        MARG['acc_max'] = x


def note_rej(x):
    if x < MARG['rej_min']:
        MARG['rej_min'] = x


def match(v, Dv, scale):
    """index of value in Dv equal to v (numerically) or -1."""
    diff = np.abs(np.asarray(Dv) - v)
    k = int(np.argmin(diff))
    if diff[k] <= TOL * scale:
        note_acc(diff[k] / scale)
        return k
    note_rej(diff[k] / scale)
    return -1


def capacity_ok(counts_core, counts_other):
    """injective sizes in 1..8, core values need >=3."""
    items = [(c, True) for c in counts_core] + [(c, False) for c in counts_other]
    if len(items) > 8:
        return False
    # bipartite matching, tiny
    sizes = list(range(1, 9))
    match_of = {}

    def try_assign(i, seen):
        c, core = items[i]
        for s in sizes:
            if s < c or (core and s < 3) or s in seen:
                continue
            seen.add(s)
            if s not in match_of or try_assign(match_of[s], seen):
                match_of[s] = i
                return True
        return False
    for i in range(len(items)):
        if not try_assign(i, set()):
            return False
    return True


def dedupe_points(X, prov, core, scale):
    keep = []
    keepP = []
    s = np.sqrt(scale)
    for x, p in zip(X, prov):
        if not np.all(np.isfinite(x)):
            continue
        if np.min(np.sqrt(((core - x) ** 2).sum(1))) <= 1e-7 * s:
            continue
        if keep:
            dd = np.sqrt(((np.array(keep) - x) ** 2).sum(1))
            if np.min(dd) <= 1e-7 * s:
                continue
        keep.append(x)
        keepP.append(p)
    return np.array(keep).reshape(-1, 2), keepP


def class_sizes_num(P, scale):
    vals = np.array([((P[i] - P[j]) ** 2).sum() for i, j in itertools.combinations(range(len(P)), 2)])
    reps, idx = cluster_vals(vals, scale)
    return sorted(collections.Counter(idx).values())


# --------------------------- per-embedding search -------------------------------
def search_embedding(core, Dv, m, scale):
    """core: 6x2 Euclidean; Dv: core values (list). returns list of (points 9x2, provenance)."""
    tri = list(itertools.combinations(range(6), 3))
    nv = len(Dv)
    Dv_a = np.array(Dv)
    grid = np.array(list(itertools.product(range(nv), repeat=3)))
    Xs, provs = [], []
    for t in tri:
        a, b, c = core[list(t)]
        X, res = trilat_vec(a, b, c, Dv_a[grid[:, 0]], Dv_a[grid[:, 1]], Dv_a[grid[:, 2]])
        ok = np.abs(res) <= TOL * scale
        STATS['trilat_trials'] += len(res)
        rr = np.abs(res[~ok]) / scale
        rr = rr[np.isfinite(rr)]
        if len(rr):
            note_rej(rr.min())
        if ok.any():
            note_acc((np.abs(res[ok]) / scale).max())
        for k in np.nonzero(ok)[0]:
            Xs.append(X[k]); provs.append(('T', t, tuple(grid[k])))
        if m == 1:
            Xs.append(circum(a, b, c)); provs.append(('C', t))
    X, prov = dedupe_points(np.array(Xs), provs, core, scale)
    STATS['cand_distinct'] += len(X)
    # classify core-partner values
    info = []
    for x in X:
        cl = [match(((x - p) ** 2).sum(), Dv, scale) for p in core]
        info.append(cl)
    nd = np.array([sum(1 for c in cl if c < 0) for cl in info]) if len(X) else np.zeros(0, int)
    found = []
    if m == 0:
        good = [i for i in range(len(X)) if nd[i] <= 3]
        STATS['cand_budget_ok'] += len(good)
        pairs = []
        for i, j in itertools.combinations(good, 2):
            e = match(((X[i] - X[j]) ** 2).sum(), Dv, scale) < 0
            if nd[i] + nd[j] + e <= 3:
                pairs.append((i, j, e))
        STATS['pairs_ok'] += len(pairs)
        pe = {(i, j): e for i, j, e in pairs}
        for i, j, e in pairs:
            for k in good:
                if k <= j:
                    continue
                if (i, k) not in pe or (j, k) not in pe:
                    continue
                if nd[i] + nd[j] + nd[k] + e + pe[(i, k)] + pe[(j, k)] != 3:
                    continue
                STATS['triples_budget'] += 1
                P = np.vstack([core, X[[i, j, k]]])
                if class_sizes_num(P, scale) == list(range(1, 9)):
                    found.append((P, [prov[i], prov[j], prov[k]], None))
    else:
        good = list(range(len(X)))
        STATS['cand_budget_ok'] += len(good)
        cc = collections.Counter(c for cl in info for c in cl if c >= 0)
        pairs = []
        for i, j in itertools.combinations(good, 2):
            Q = np.vstack([core, X[[i, j]]])
            vals = np.array([((Q[a] - Q[b]) ** 2).sum() for a, b in itertools.combinations(range(8), 2)])
            reps, idx = cluster_vals(vals, scale)
            cnt = collections.Counter(idx)
            core_ids = set(idx[k] for k, (a, b) in enumerate(itertools.combinations(range(8), 2)) if a < 6 and b < 6)
            if not capacity_ok([cnt[r] for r in core_ids], [cnt[r] for r in cnt if r not in core_ids]):
                continue
            if len(cnt) - len(core_ids) > 3:
                continue
            pairs.append((i, j, Q, reps))
        STATS['pairs_ok'] += len(pairs)
        for i, j, Q, reps in pairs:
            # w candidates: X (good points) + all 56 trilaterations over D(Q)^3 + 56 circumcentres
            Ws, wprov = [], []
            for k in range(len(X)):
                if k != i and k != j:
                    Ws.append(X[k]); wprov.append(prov[k])
            R = np.array(reps)
            g = np.array(list(itertools.product(range(len(R)), repeat=3)))
            for t in itertools.combinations(range(8), 3):
                a, b, c = Q[list(t)]
                if max(t) < 6:
                    Xw, res = trilat_vec(a, b, c, R[g[:, 0]], R[g[:, 1]], R[g[:, 2]])
                else:
                    # w has >=5 D(Q)-partners in Q, at most 2 outside S -> >=3 in S: core triangles suffice
                    Xw, res = np.zeros((0, 2)), np.zeros(0)
                ok = np.abs(res) <= TOL * scale
                STATS['last_trials'] += len(res)
                rr = np.abs(res[~ok]) / scale
                rr = rr[np.isfinite(rr)]
                if len(rr):
                    note_rej(rr.min())
                for kk in np.nonzero(ok)[0]:
                    Ws.append(Xw[kk]); wprov.append(('LT', t, tuple(R[g[kk]])))
                Ws.append(circum(a, b, c)); wprov.append(('LC', t))
            W, wp = dedupe_points(np.array(Ws), wprov, Q, scale)
            STATS['last_cands'] += len(W)
            for w, p in zip(W, wp):
                P = np.vstack([Q, w])
                if class_sizes_num(P, scale) == list(range(1, 9)):
                    found.append((P, [prov[i], prov[j], p], None))
    return found


def exact_rebuild(K, G, cpts, z, pattern, provs):
    """rebuild outside points exactly from provenance. z: exact core label values (list),
    Dv order = sorted labels as used numerically (labels 1..r -> z[l-1])."""
    pts = list(cpts)
    for p in provs:
        if p[0] == 'T':
            t, vals = p[1], p[2]
            x = ex_trilat(K, G, pts[t[0]], pts[t[1]], pts[t[2]], z[vals[0]], z[vals[1]], z[vals[2]])
        elif p[0] == 'C':
            t = p[1]
            x = ex_circum(K, G, pts[t[0]], pts[t[1]], pts[t[2]])
        elif p[0] == 'LC':
            t = p[1]
            x = ex_circum(K, G, pts[t[0]], pts[t[1]], pts[t[2]])
        else:  # LT: values given numerically -> find exact value among exact pair distances of Q
            return None
        if x is None:
            return None
        pts.append(x)
    return pts


def main():
    inp, out = sys.argv[1], sys.argv[2]
    limit = None
    stride, offset = 1, 0
    for a in sys.argv:
        if a.startswith('--limit='):
            limit = int(a.split('=')[1])
        if a.startswith('--stride='):
            stride = int(a.split('=')[1])
        if a.startswith('--offset='):
            offset = int(a.split('=')[1])
    done = set()
    if os.path.exists(out + '.jsonl'):
        for l in open(out + '.jsonl'):
            done.add(json.loads(l)['key'])
    fo = open(out + '.jsonl', 'a')
    t0 = time.time()
    with gzip.open(inp, 'rt') as f:
        for ln, line in enumerate(f):
            if ln % stride != offset:
                continue
            if limit is not None and ln >= limit:
                break
            d = json.loads(line)
            pat = [int(c) for c in d['pattern']]
            r = d['r']
            for fi, fac in enumerate(d['factors']):
                key = '%s#%d' % (d['pattern'], fi)
                if key in done:
                    continue
                rec = process_factor(d, pat, r, fac, key)
                fo.write(json.dumps(rec) + '\n')
                fo.flush()
            if ln % 200 == 0:
                print(ln, '%.0fs' % (time.time() - t0), dict(STATS), MARG, flush=True)
    print('DONE', '%.0fs' % (time.time() - t0), dict(STATS), MARG, flush=True)


def process_factor(d, pat, r, fac, key):
    K = KQ(fac['q'])
    rec = dict(key=key, status=fac['status'], deg=fac['deg'])
    rts = real_roots(K.q)
    rec['nreal_vb'] = len(rts)
    rec['nreal_builder'] = fac.get('nreal')
    if 'G' not in fac or fac['G'] is None:
        rec['pd_vb'] = None
        return rec
    G = [K.el(c) for c in fac['G']]
    pd = []
    for i, th in enumerate(rts):
        g11, g12, g22 = (ev(x, th) for x in G)
        det = g11 * g22 - g12 * g12
        if g11 > 0 and det > 0:
            pd.append(i)
        elif not (g11 < 0 or det < 0):
            rec.setdefault('pd_undecided', []).append(i)
    rec['pd_vb'] = len(pd)
    if fac['status'] != 'REALIZED':
        return rec
    rec['pd_builder'] = len(fac['pd_roots'])
    z = [K.el(c) for c in fac['z']]
    cpts = [(K.el(p[0]), K.el(p[1])) for p in fac['points']]
    # exact core verification: 15 distances = z[label], GP
    ok = True
    pi = [(i, j) for i in range(6) for j in range(i + 1, 6)]
    for (i, j), l in zip(pi, pat):
        if not K.iszero(ex_d2(K, G, cpts[i], cpts[j]) - z[l - 1]):
            ok = False
    chk = ex_check_config(K, G, cpts)
    rec['core_exact_ok'] = ok and not chk.get('dup') and not chk['collinear'] and not chk['concyclic'] \
        and len(chk['sizes']) == r
    m = 6 - r
    rec['m'] = m
    found_all = []
    for i in pd:
        th = rts[i]
        Gn = np.array([[to_f(ev(G[0], th)), to_f(ev(G[1], th))], [to_f(ev(G[1], th)), to_f(ev(G[2], th))]])
        L = np.linalg.cholesky(Gn)
        core = np.array([[to_f(ev(p[0], th)), to_f(ev(p[1], th))] for p in cpts]) @ L
        Dv = [to_f(ev(zz, th)) for zz in z]
        scale = max(Dv)
        STATS['embeddings'] += 1
        found = search_embedding(core, Dv, m, scale)
        for P, provs, _ in found:
            STATS['num_sizes_ok'] += 1
            # exact re-check
            # map numeric provenance values (indices into Dv) -> exact z
            pts = exact_rebuild(K, G, cpts, z, pat, provs)
            if pts is None:
                # LT provenance: rebuild via trying all exact pair values of Q
                pts = rebuild_lt(K, G, cpts, z, provs, P, L)
            ex = ex_check_config(K, G, pts) if pts is not None else None
            item = dict(root=i, provs=[list(map(str, p)) for p in provs])
            if ex is None:
                item['exact'] = 'REBUILD_FAILED'
                STATS['rebuild_failed'] += 1
            else:
                item['sizes'] = ex.get('sizes')
                item['collinear'] = ex.get('collinear')
                item['concyclic'] = ex.get('concyclic')
                item['crescent'] = (ex.get('sizes') == list(range(1, 9)) and not ex['collinear'] and not ex['concyclic'])
                item['dvals'] = sorted(to_f(ev(v, th)) for v in ex['values']) if 'values' in ex else None
                item['points'] = [[str(c) for c in (p[0].coeffs())] + ['|'] + [str(c) for c in p[1].coeffs()] for p in pts]
                if item['sizes'] == list(range(1, 9)):
                    STATS['exact_sizes_ok'] += 1
                if item['crescent']:
                    STATS['CRESCENT'] += 1
            found_all.append(item)
    rec['found'] = found_all
    return rec


def rebuild_lt(K, G, cpts, z, provs, Pnum, L):
    pts = exact_rebuild(K, G, cpts, z, None, provs[:2])
    if pts is None:
        return None
    p = provs[2]
    if p[0] == 'LT' or p[0] == 'LC':
        t = p[1]
        if p[0] == 'LC':
            x = ex_circum(K, G, pts[t[0]], pts[t[1]], pts[t[2]])
            return pts + [x]
        # exact values of Q pairs; choose those matching numerically
        vals = {}
        for a, b in itertools.combinations(range(8), 2):
            vals[(a, b)] = ex_d2(K, G, pts[a], pts[b])
        target = p[2]
        Qn = Pnum[:8]
        choices = []
        for tv in target:
            cands = [ab for ab in vals if abs(((Qn[ab[0]] - Qn[ab[1]]) ** 2).sum() - tv) <= 1e-6 * max(1, abs(tv))]
            choices.append(vals[cands[0]] if cands else None)
        if None in choices:
            return None
        x = ex_trilat(K, G, pts[t[0]], pts[t[1]], pts[t[2]], *choices)
        return None if x is None else pts + [x]
    pts2 = exact_rebuild(K, G, cpts, z, None, provs)
    return pts2


if __name__ == '__main__':
    main()
