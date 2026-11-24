"""V-b independent coverage search (Task 2).

Design (differs from NOTES/coverage.md (b)):
 * own arrangement representatives, low pairs placed on the LAST points;
 * point-by-point extension: for each new point v, a 5-point base B0 of high partners is joined
   against a lookup table built from U's orbit (Wang-style five-point table, NOT used by builder);
   remaining partners are then assigned one by one with completed-6-subset checks;
 * after each completed point the partial pattern is reduced modulo the stabiliser of the
   arrangement restricted to the present points (orbit dedup at every level; builder keeps all
   labelled states);
 * full local tests (Lemma 2.1/2.2) for every point, capacity (Lemma 5.2), U-membership;
 * parallelogram affine tests (Lemma 2.3) on complete 9-point patterns (and optionally per level).
Usage: python coverage_vb.py U_file out_prefix [--affine-each-level]
"""
import sys, itertools, time, json, collections
sys.path.insert(0, __import__('os').path.dirname(__file__))
from common_vb import *

N = 9
P9 = PI[9]
LOW1, LOW2 = -1, -2

# 0-based representatives (my own). (d1 pair, [d2 pairs])
ARR = {
    'disjoint':        ((0, 6), [(1, 7), (2, 8)]),
    'triangle':        ((7, 8), [(6, 7), (6, 8)]),
    'star':            ((0, 8), [(1, 8), (2, 8)]),
    'chain_d1_middle': ((7, 8), [(0, 7), (1, 8)]),
    'chain_d1_end':    ((0, 7), [(7, 8), (1, 8)]),
    'P3K2_d1_disjoint': ((2, 7), [(0, 8), (1, 8)]),
    'P3K2_d1_meeting': ((0, 8), [(1, 8), (2, 7)]),
}


def arrangement_invariant(d1, d2s):
    edges = [d1] + list(d2s)
    deg = collections.Counter(x for e in edges for x in e)
    shape = tuple(sorted(deg.values()))
    return (shape, tuple(sorted((deg[d1[0]], deg[d1[1]]))), len(deg))


def verify_arrangements():
    """all 36*C(35,2) labelled low structures; classify by an invariant that determines the
    isomorphism type of (graph with 3 edges, marked d1 edge); check 7 classes and that each rep
    hits one class."""
    pairs = PAIRS[9]
    cnt = collections.Counter()
    for d1 in pairs:
        rest = [p for p in pairs if p != d1]
        for d2s in itertools.combinations(rest, 2):
            cnt[arrangement_invariant(d1, d2s)] += 1
    reps = {name: arrangement_invariant(*ARR[name]) for name in ARR}
    return cnt, reps


def stabiliser(d1, d2s, k):
    """restrictions to {0..k-1} of perms of 9 points preserving the low structure and {0..k-1}."""
    d1s = frozenset(d1)
    d2set = {frozenset(e) for e in d2s}
    low_pts = set(d1) | {x for e in d2s for x in e}
    # perms of non-low points are free (within the prefix / outside prefix); enumerate perms of
    # the low points explicitly and of the rest by prefix blocks
    lowl = sorted(low_pts)
    res = set()
    for img in itertools.permutations(range(9), len(lowl)):
        g = dict(zip(lowl, img))
        if frozenset((g[d1[0]], g[d1[1]])) != d1s:
            continue
        if {frozenset((g[a], g[b])) for a, b in d2s} != d2set:
            continue
        # images must lie in low_pts (structure preserved implies this) and prefix preserved
        if any((x < k) != (g[x] < k) for x in lowl):
            continue
        free_in = [x for x in range(k) if x not in low_pts]
        free_tgt = sorted(set(range(k)) - {g[x] for x in lowl if x < k})
        if len(free_in) != len(free_tgt):
            continue
        for fi in itertools.permutations(free_tgt):
            gg = dict(g)
            gg.update(zip(free_in, fi))
            res.add(tuple(gg[x] for x in range(k)))
    return sorted(res)


class Search:
    def __init__(self, Upats, name, affine_each=False):
        self.name = name
        d1, d2s = ARR[name]
        self.d1, self.d2s = d1, d2s
        self.lowlab = {}
        self.lowlab[frozenset(d1)] = LOW1
        for e in d2s:
            self.lowlab[frozenset(e)] = LOW2
        self.affine_each = affine_each
        # U orbit set and five-point lookup table
        self.Uorb = set()
        for u in Upats:
            self.Uorb |= orbit6(list(u))
        self.table = collections.defaultdict(set)
        pi6 = PI[6]
        old = [(i, j) for i in range(5) for j in range(i + 1, 5)]
        for s in self.Uorb:
            seq = [s[pi6[p]] for p in old] + [s[pi6[(i, 5)]] for i in range(5)]
            rs = rename(seq)
            self.table[rs[:10]].add(rs[10:])
        self.stab = {k: stabiliser(d1, d2s, k) for k in range(6, 10)}
        # all-high prefix
        self.low_pts = set(d1) | {x for e in d2s for x in e}

    def is_low(self, a, b):
        return frozenset((a, b)) in self.lowlab

    def all_high(self, pts):
        return not any(self.is_low(a, b) for a, b in itertools.combinations(pts, 2))

    # ---- canonical form of a partial state on k points under the stabiliser ----
    def canon_state(self, lab, k):
        pi = PI[k]
        best = None
        pairs = PAIRS[k]
        for g in self.stab[k]:
            out = [None] * len(pairs)
            for idx, (i, j) in enumerate(pairs):
                out[pi[(g[i], g[j])]] = lab[idx]
            # rename high labels by first occurrence; low labels kept
            m = {}
            r = []
            for x in out:
                if x < 0:
                    r.append(x)
                else:
                    if x not in m:
                        m[x] = len(m)
                    r.append(m[x])
            r = tuple(r)
            if best is None or r < best:
                best = r
        return best

    def in_U(self, lab9, pts):
        """lab9: dict-like accessor c(a,b)"""
        return rename([lab9(a, b) for a, b in itertools.combinations(pts, 2)]) in self.Uorb

    # ---- local tests involving new point v on full pattern of points 0..v ----
    def local_ok_v(self, c, v):
        n = v + 1
        # L2.1(i)
        cnt = collections.Counter(c(x, v) for x in range(v))
        if max(cnt.values()) > 3:
            return False
        for a in range(v):
            l = c(a, v)
            k = sum(1 for x in range(n) if x != a and c(a, x) == l)
            if k > 3:
                return False
        # L2.1(ii),(iii) with a,b < v (x = v is new)
        for a in range(v):
            for b in range(a + 1, v):
                key = (c(a, v), c(b, v))
                k = 0
                e = 0
                for x in range(n):
                    if x == a or x == b:
                        continue
                    cax, cbx = c(a, x), c(b, x)
                    if (cax, cbx) == key:
                        k += 1
                    if cax == cbx:
                        e += 1
                if k > 2 or e > 2:
                    return False
        # a = v
        for b in range(v):
            cnt2 = collections.Counter()
            e = 0
            for x in range(v):
                if x == b:
                    continue
                key = (c(v, x), c(b, x))
                cnt2[key] += 1
                if cnt2[key] > 2:
                    return False
                if key[0] == key[1]:
                    e += 1
                    if e > 2:
                        return False
        # L2.2 quadruples containing v
        for a, b, cc in itertools.combinations(range(v), 3):
            if c(a, b) == c(cc, v) and c(a, cc) == c(b, v) and c(b, cc) == c(a, v):
                return False
        return True

    @staticmethod
    def capacity_ok(lab, final=False):
        cnt = collections.Counter(x for x in lab if x >= 0)
        ms = sorted(cnt.values())
        s = len(ms)
        if s > 6:
            return False
        if final:
            return ms == [3, 4, 5, 6, 7, 8]
        return all(m <= 8 - s + i + 1 for i, m in enumerate(ms))

    def extend(self, lab, v):
        """lab: tuple over lex pairs of 0..v-1. yield label tuples over 0..v."""
        k = v
        piold = PI[k] if k >= 2 else None
        cold = lambda a, b: lab[piold[(a, b)]]
        used = sorted(set(x for x in lab if x >= 0))
        nused = len(used)
        partners = list(range(v))
        lowp = {a: self.lowlab[frozenset((a, v))] for a in partners if self.is_low(a, v)}
        highp = [a for a in partners if a not in lowp]
        B0 = None
        for B in itertools.combinations(highp, 5):
            if self.all_high(B):
                B0 = B
                break
        assert B0 is not None
        old = [(B0[i], B0[j]) for i in range(5) for j in range(i + 1, 5)]
        oldlabs = [cold(a, b) for a, b in old]
        m = {}
        for x in oldlabs:
            if x not in m:
                m[x] = len(m)
        inv = {r: g for g, r in m.items()}
        t = len(m)
        Bset = set(oldlabs)
        notB = [g for g in used if g not in Bset]
        entries = self.table.get(tuple(m[x] for x in oldlabs), ())
        rest = [a for a in highp if a not in B0]
        for ent in entries:
            nnew = max([e - t + 1 for e in ent if e >= t], default=0)
            # injective assignment of new entry labels to global labels not on B0
            for choice in self._new_assign(nnew, notB, nused):
                gl = {}
                for e, a in zip(ent, B0):
                    gl[a] = inv[e] if e < t else choice[e - t]
                total = nused + sum(1 for g in choice if g >= nused)
                yield from self._fill_rest(lab, v, gl, lowp, rest, 0, total, B0)

    @staticmethod
    def _new_assign(nnew, notB, nused):
        res = []

        def rec(j, cur, fresh):
            if j == nnew:
                res.append(tuple(cur))
                return
            for g in notB:
                if g not in cur:
                    rec(j + 1, cur + [g], fresh)
            if fresh < 6:
                rec(j + 1, cur + [fresh], fresh + 1)
        rec(0, [], nused)
        return res

    def _fill_rest(self, lab, v, gl, lowp, rest, idx, total, B0):
        if idx == len(rest):
            yield from self._finish(lab, v, gl, lowp)
            return
        a = rest[idx]
        assigned_high = list(B0) + rest[:idx]
        piold = PI[v]
        for g in list(range(total)) + ([total] if total < 6 else []):
            gl[a] = g
            ok = True
            # count at v of label g <= 3
            if sum(1 for x in gl.values() if x == g) > 3:
                ok = False
            if ok:
                # completed all-high 6-subsets {v, a} + 4 of assigned_high
                for S4 in itertools.combinations(assigned_high, 4):
                    pts = tuple(sorted(S4 + (a,)))
                    if not self.all_high(pts):
                        continue
                    # build lex order over sorted 6 points pts+(v,), v is largest
                    full = []
                    P6 = pts + (v,)
                    for i in range(6):
                        for j in range(i + 1, 6):
                            x, y = P6[i], P6[j]
                            full.append(gl[x] if y == v else lab[piold[(x, y)]])
                    if rename(full) not in self.Uorb:
                        ok = False
                        break
            if ok:
                yield from self._fill_rest(lab, v, gl, lowp, rest, idx + 1, max(total, g + 1), B0)
            del gl[a]

    def _finish(self, lab, v, gl, lowp):
        pin = PI[v + 1]
        piold = PI[v]
        new = [None] * (len(lab) + v)
        for (i, j) in PAIRS[v + 1]:
            if j == v:
                new[pin[(i, j)]] = lowp[i] if i in lowp else gl[i]
            else:
                new[pin[(i, j)]] = lab[piold[(i, j)]]
        c = lambda a, b: new[pin[(a, b)]]
        if not self.capacity_ok(new, final=(v == 8)):
            return
        if not self.local_ok_v(c, v):
            return
        if self.affine_each and affine_reject(new, v + 1) is not None:
            return
        yield tuple(new)

    def images(self, lab, k):
        pi = PI[k]
        pairs = PAIRS[k]
        out_all = set()
        for g in self.stab[k]:
            out = [None] * len(pairs)
            for idx, (i, j) in enumerate(pairs):
                out[pi[(g[i], g[j])]] = lab[idx]
            m = {}
            r = []
            for x in out:
                if x < 0:
                    r.append(x)
                else:
                    if x not in m:
                        m[x] = len(m)
                    r.append(m[x])
            out_all.add(tuple(r))
        return out_all

    def seeds(self, Upats):
        reps = set()
        for u in Upats:
            seen = set()
            for s in orbit6(list(u)):
                if s in seen:
                    continue
                im = self.images(s, 6)
                seen |= im
                reps.add(min(im))
        return reps

    def run(self, Upats, log):
        t0 = time.time()
        level = {6: self.seeds(Upats)}
        if self.affine_each:
            level[6] = {s for s in level[6] if affine_reject(list(s), 6) is None}
        log('%s level6 %d (%.1fs)' % (self.name, len(level[6]), time.time() - t0))
        for v in (6, 7, 8):
            nxt = set()
            raw = 0
            for st in level[v]:
                for new in self.extend(st, v):
                    raw += 1
                    nxt.add(self.canon_state(new, v + 1))
            level[v + 1] = nxt
            log('%s level%d raw %d dedup %d (%.1fs)' % (self.name, v + 1, raw, len(nxt), time.time() - t0))
        return level


def to_sizes(lab9):
    cnt = collections.Counter(lab9)
    return [1 if x == LOW1 else 2 if x == LOW2 else cnt[x] for x in lab9]


def seven_point_canon(level7):
    """canonicalise level-7 states (all points high) up to S7, plus affine filter."""
    perms7 = list(itertools.permutations(range(7)))
    canon = set()
    for f in level7:
        assert min(f) >= 0
        canon.add(min(rename(relabel_points(list(f), 7, p)) for p in perms7))
    aff = {c for c in canon if affine_reject(list(c), 7) is None}
    return canon, aff


if __name__ == '__main__':
    Ufile, out = sys.argv[1], sys.argv[2]
    affine_each = '--affine-each-level' in sys.argv
    only = [a.split('=')[1] for a in sys.argv if a.startswith('--only=')]
    Upats = read_pattern_file(Ufile)
    logf = open(out + '.log', 'a')

    def log(s):
        print(s, flush=True)
        logf.write(s + '\n'); logf.flush()
    log('U lines %d, canonical distinct %d' % (len(Upats), len({canon6(list(u)) for u in Upats})))
    cnt, reps = verify_arrangements()
    log('arrangement classes: %s' % sorted(cnt.values()))
    log('reps -> class sizes: %s' % {n: cnt[r] for n, r in reps.items()})
    assert len(cnt) == 7 and len(set(reps.values())) == 7
    if '--level7' in sys.argv:
        t0 = time.time()
        S = Search(Upats, 'triangle')
        l6 = S.seeds(Upats)
        l7 = set()
        for st in l6:
            for new in S.extend(st, 6):
                l7.add(S.canon_state(new, 7))
        c7, a7 = seven_point_canon(l7)
        log('triangle level7 states %d' % len(l7))
        log('7-point stage: %d canonical (no affine), %d after affine (%.0fs)' % (len(c7), len(a7), time.time() - t0))
        with open(out + '_level7.txt', 'w') as f:
            for c in sorted(c7):
                f.write(''.join(map(str, c)) + (' affine_ok' if c in a7 else ' affine_rejected') + '\n')
    summary = {}
    with open(out + '_found.txt', 'a') as fo:
        for name in ARR:
            if only and name not in only:
                continue
            S = Search(Upats, name, affine_each)
            level = S.run(Upats, log)
            comp = set()
            for st in level[9]:
                comp.add(canon_complete(to_sizes(st), 9))
            aff_ok = []
            for cf in sorted(comp):
                if affine_reject(list(cf[1]), 9) is None:
                    aff_ok.append(cf)
            for cf in sorted(comp):
                fo.write('%s %s\n' % (name, ''.join(map(str, cf[1]))))
            summary[name] = dict(levels={k: len(vv) for k, vv in level.items()},
                                 distinct9=len(comp), after_affine=len(aff_ok))
            log('%s RESULT %s' % (name, json.dumps(summary[name])))
    log('SUMMARY ' + json.dumps(summary))
