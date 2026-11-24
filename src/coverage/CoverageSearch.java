import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import java.util.zip.*;

/**
 * Builder C: independent coverage search (re-implementation of the role of Wang's Prop 5.1 /
 * Computation 5.4, but with our own design).
 *
 * Question decided: does there exist a labelled distance pattern on n points (n = 7, 8, 9) such that
 *   - the three low pairs carry the labels D1 (one pair) and D2 (two pairs) in one of the 7 low-pair
 *     arrangements (representatives derived and checked in arrangements.py),
 *   - every other pair carries one of at most H = n-3 high labels,
 *   - every k-subset (k = n-3) spanning only high pairs has a pattern equivalent to a member of U,
 *   - Lemma 2.1 (i)-(iii) and Lemma 2.2 hold (on all labels, low labels included),
 *   - the final high-label class sizes are exactly {3, ..., n-1}?
 *
 * Design (no five-point lookup, no placement symmetry reduction):
 *   - SEEDS: positions 0..k-1 never span a low pair (in all 7 representatives). Their labelled pattern
 *     (high labels renamed by first occurrence in lex pair order) must lie in the ORBIT of some U member
 *     under all k! point permutations. We enumerate every distinct element of every orbit.
 *   - GROWTH: positions k..n-1 are added one at a time; the pairs (u,v), u = 0..v-1, are assigned in
 *     that order. A low pair gets its fixed low label. A high pair gets an existing high label or the
 *     next unused one (first-occurrence renaming, <= H labels).
 *   - PRUNING after every assignment of pair (u,v): incremental local tests on configurations that
 *     contain the pair (u,v); the capacity test (Wang Lemma 5.2 / Hall); and U-membership of every
 *     all-high k-subset whose last two elements are u < v (exactly those subsets become complete).
 *   - LEAF: high class sizes exactly {3..n-1}; paranoid full re-check of everything.
 *
 * Usage:
 *   java CoverageSearch --n 9 --U data/E1_U.txt.gz [--threads 1] [--arr all|name,...] [--mode auto|orbit|canon]
 *        [--subset N --rseed S] [--nolocal] [--paranoid] [--out file] [--level7 file] [--maxlabels6 H]
 */
public class CoverageSearch {
    static final int L1 = 6, L2 = 7, NL = 8; // label codes: 0..5 high, 6 = delta1, 7 = delta2

    // ---- arrangements (1-based, for n = 9; shifted by n-9 for smaller n). {delta1, delta2, delta2} ----
    static final String[] ARR_NAMES = {"disjoint", "triangle", "star", "chain_d1_middle", "chain_d1_end",
            "P3K2_d1_disjoint", "P3K2_d1_meeting"};
    static final int[][][] ARR = {
            {{4, 7}, {5, 8}, {6, 9}},
            {{8, 9}, {7, 8}, {7, 9}},
            {{6, 9}, {7, 9}, {8, 9}},
            {{7, 9}, {6, 9}, {7, 8}},
            {{6, 9}, {7, 9}, {7, 8}},
            {{5, 8}, {6, 9}, {7, 9}},
            {{6, 9}, {7, 9}, {5, 8}},
    };

    // ---- configuration ----
    final int n, k, H, kp;     // kp = C(k,2)
    final int[][] lowLab;      // lowLab[i][j] = L1/L2 or -1 (high)
    final boolean useLocal, paranoid;
    final int[] pu, pv, plow;  // remaining pair order
    final int[][][] subsets;   // subsets[idx] = list of k-subsets (sorted) completed by pair idx
    final int[][][] lookSubs;  // lookSubs[idx] = all-high k-subsets T = S+{w,v}, u in S subset of {0..u}, u<w<v: only (w,v) unassigned
    final boolean lookahead;
    final boolean[] endOfPoint;
    final int[][] pairIdxK;    // lex pair index on k points
    final int[][] permK;       // all k! perms
    final UMember umem;

    // stats
    final AtomicLong nodes = new AtomicLong(), leaves = new AtomicLong(), seeds = new AtomicLong();
    final AtomicLongArray levelCount;
    final AtomicLong lastPrint = new AtomicLong(System.nanoTime());
    final Set<String> found = ConcurrentHashMap.newKeySet();      // labelled leaf patterns (class-size digits), capped
    static final int FOUND_CAP = 200_000;
    final AtomicBoolean foundOverflow = new AtomicBoolean(false);
    static Map<String, List<String[]>> targetsByInv = null;        // invariant -> [name, canonical]
    final Set<String> targetHits = ConcurrentHashMap.newKeySet();
    final Set<Long> level7 = ConcurrentHashMap.newKeySet();        // optional: canonical 7-point all-high patterns
    final boolean collect7;
    final int stopAfterPoint;  // if >=0, stop the search after completing this point (used for --level7)

    interface UMember { boolean contains(long renamedPacked, int[] lab); }

    CoverageSearch(int n, int arr, UMember umem, int[][] perms, boolean useLocal, boolean paranoid, boolean collect7, int stopAfterPoint, boolean lookahead) {
        this.lookahead = lookahead;
        this.n = n; this.k = n - 3; this.H = n - 3; this.kp = k * (k - 1) / 2;
        this.umem = umem; this.permK = perms; this.useLocal = useLocal; this.paranoid = paranoid;
        this.collect7 = collect7; this.stopAfterPoint = stopAfterPoint;
        lowLab = new int[n][n];
        for (int[] r : lowLab) Arrays.fill(r, -1);
        int sh = 9 - n;
        for (int e = 0; e < 3; e++) {
            int a = ARR[arr][e][0] - 1 - sh, b = ARR[arr][e][1] - 1 - sh;
            if (a < 0 || b < 0 || a >= n || b >= n) throw new IllegalStateException("bad rep");
            lowLab[a][b] = lowLab[b][a] = (e == 0) ? L1 : L2;
        }
        for (int i = 0; i < k; i++) for (int j = i + 1; j < k; j++)
            if (lowLab[i][j] >= 0) throw new IllegalStateException("seed block spans a low pair");
        pairIdxK = new int[k][k];
        { int p = 0; for (int i = 0; i < k; i++) for (int j = i + 1; j < k; j++) { pairIdxK[i][j] = pairIdxK[j][i] = p++; } }
        List<int[]> order = new ArrayList<>();
        for (int v = k; v < n; v++) for (int u = 0; u < v; u++) order.add(new int[]{u, v});
        int m = order.size();
        pu = new int[m]; pv = new int[m]; plow = new int[m]; subsets = new int[m][][]; endOfPoint = new boolean[m];
        lookSubs = new int[m][][];
        for (int idx = 0; idx < m; idx++) {
            int u = order.get(idx)[0], v = order.get(idx)[1];
            pu[idx] = u; pv[idx] = v; plow[idx] = lowLab[u][v]; endOfPoint[idx] = (u == v - 1);
            List<int[]> subs = new ArrayList<>();
            if (k - 2 <= u) {
                for (int[] s : combos(u, k - 2)) {
                    int[] t = new int[k];
                    System.arraycopy(s, 0, t, 0, k - 2);
                    t[k - 2] = u; t[k - 1] = v;
                    boolean allHigh = true;
                    for (int a = 0; a < k && allHigh; a++) for (int b = a + 1; b < k; b++) if (lowLab[t[a]][t[b]] >= 0) { allHigh = false; break; }
                    if (allHigh) subs.add(t);
                }
            }
            subsets[idx] = subs.toArray(new int[0][]);
            List<int[]> look = new ArrayList<>();
            if (k - 3 <= u) {
                for (int w = u + 1; w < v; w++) for (int[] s : combos(u, k - 3)) {
                    int[] t = new int[k];
                    System.arraycopy(s, 0, t, 0, k - 3);
                    t[k - 3] = u; t[k - 2] = w; t[k - 1] = v;
                    boolean allHigh = true;
                    for (int a = 0; a < k && allHigh; a++) for (int b = a + 1; b < k; b++) if (lowLab[t[a]][t[b]] >= 0) { allHigh = false; break; }
                    if (allHigh) look.add(t);
                }
            }
            lookSubs[idx] = look.toArray(new int[0][]);
        }
        levelCount = new AtomicLongArray(n + 1);
    }

    static List<int[]> combos(int m, int r) { // r-subsets of {0..m-1}, sorted
        List<int[]> out = new ArrayList<>();
        combRec(new int[r], 0, 0, m, out);
        return out;
    }
    static void combRec(int[] cur, int pos, int start, int m, List<int[]> out) {
        if (pos == cur.length) { out.add(cur.clone()); return; }
        for (int x = start; x < m; x++) { cur[pos] = x; combRec(cur, pos + 1, x + 1, m, out); }
    }

    // ================= search state (one per thread) =================
    final class State {
        final int[][] c = new int[n][n];
        final int[] cnt = new int[NL];
        int t; // number of high labels used
        final int[] subLab = new int[kp];

        void loadSeed(int[] lab) {
            for (int[] r : c) Arrays.fill(r, -1);
            Arrays.fill(cnt, 0);
            t = 0;
            int p = 0;
            for (int i = 0; i < k; i++) for (int j = i + 1; j < k; j++) {
                int l = lab[p++];
                c[i][j] = c[j][i] = l; cnt[l]++; t = Math.max(t, l + 1);
            }
        }

        long myNodes = 0;
        void dfs(int idx) {
            myNodes++;
            if (idx == pu.length) { leaf(); return; }
            int u = pu[idx], v = pv[idx];
            if (plow[idx] >= 0) {
                set(u, v, plow[idx]);
                if (ok(idx, false)) after(idx);
                unset(u, v);
                return;
            }
            int maxl = Math.min(t, H - 1);
            for (int l = 0; l <= maxl; l++) {
                int oldT = t;
                set(u, v, l);
                if (l == t) t++;
                if (ok(idx, true)) after(idx);
                t = oldT;
                unset(u, v);
            }
        }

        void after(int idx) {
            if (endOfPoint[idx]) {
                int v = pv[idx];
                levelCount.incrementAndGet(v + 1);
                if (collect7 && v == 6) level7.add(canon7());
                if (v == stopAfterPoint) return;
            }
            dfs(idx + 1);
        }

        void set(int u, int v, int l) { c[u][v] = c[v][u] = l; cnt[l]++; }
        void unset(int u, int v) { cnt[c[u][v]]--; c[u][v] = c[v][u] = -1; }

        boolean ok(int idx, boolean high) {
            int u = pu[idx], v = pv[idx];
            if (useLocal && !localInc(c, n, u, v)) { if (paranoid) assertFull(false); return false; }
            if (paranoid && useLocal && !fullLocal(c, n)) throw new IllegalStateException("incremental local test missed a violation");
            if (high && !capacityOk()) return false;
            for (int[] s : subsets[idx]) if (!inU(s)) return false;
            if (lookahead) for (int[] s : lookSubs[idx]) if (!completable(s)) return false;
            return true;
        }

        /** T sorted, all pairs assigned except the last lex pair (T[k-2],T[k-1]). True iff some label for it
         *  (a label already on T, or a label not on T) puts T in U. */
        boolean completable(int[] s) {
            Arrays.fill(ren, -1);
            int nx = 0, p = 0;
            long packed = 0;
            for (int a = 0; a < k; a++) for (int b = a + 1; b < k; b++) {
                if (a == k - 2 && b == k - 1) break;
                int l = c[s[a]][s[b]];
                if (l < 0 || l >= H) throw new IllegalStateException("lookahead subset not all-high/assigned");
                int r = ren[l];
                if (r < 0) { r = nx++; ren[l] = r; }
                subLab[p++] = r;
                packed = (packed << 3) | r;
            }
            if (c[s[k - 2]][s[k - 1]] != -1) throw new IllegalStateException("lookahead pair already assigned");
            int maxr = Math.min(nx, H - 1);
            for (int r = 0; r <= maxr; r++) {
                subLab[kp - 1] = r;
                if (umem.contains((packed << 3) | r, subLab)) return true;
            }
            return false;
        }
        final int[] ren = new int[NL];
        final int[] capm = new int[NL];

        void assertFull(boolean expect) {
            if (fullLocal(c, n) != expect) throw new IllegalStateException("incremental/full local mismatch");
        }

        boolean capacityOk() {
            int[] m = capm;
            int s = 0;
            for (int l = 0; l < H; l++) if (cnt[l] > 0) {
                int x = cnt[l], i = s++;
                while (i > 0 && m[i - 1] > x) { m[i] = m[i - 1]; i--; } // insertion sort
                m[i] = x;
            }
            // available multiplicities b_1..b_H = 3..n-1 ; need m_i <= b_{H-s+i}  (1-based)
            for (int i = 0; i < s; i++) if (m[i] > 3 + (H - s + i)) return false;
            return true;
        }

        boolean inU(int[] s) {
            Arrays.fill(ren, -1);
            int nx = 0, p = 0;
            long packed = 0;
            for (int a = 0; a < k; a++) for (int b = a + 1; b < k; b++) {
                int l = c[s[a]][s[b]];
                if (l < 0 || l >= H) throw new IllegalStateException("subset not all-high/assigned");
                int r = ren[l];
                if (r < 0) { r = nx++; ren[l] = r; }
                subLab[p++] = r;
                packed = (packed << 3) | r;
            }
            return umem.contains(packed, subLab);
        }

        long myLeaves = 0;
        void leaf() {
            int[] m = new int[H];
            for (int l = 0; l < H; l++) m[l] = cnt[l];
            Arrays.sort(m);
            for (int i = 0; i < H; i++) if (m[i] != 3 + i) return;
            if (cnt[L1] != 1 || cnt[L2] != 2) throw new IllegalStateException("low counts");
            if (paranoid) { // full re-check of every condition with independent code paths
                if (useLocal && !fullLocal(c, n)) throw new IllegalStateException("leaf fails full local test");
                for (int[] s : combos(n, k)) {
                    boolean allHigh = true;
                    for (int a = 0; a < k && allHigh; a++) for (int b = a + 1; b < k; b++) if (c[s[a]][s[b]] >= H) { allHigh = false; break; }
                    if (allHigh && !inU(s)) throw new IllegalStateException("leaf has a core outside U");
                }
            }
            myLeaves++;
            char[] sb = new char[n * (n - 1) / 2];
            int q = 0;
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
                int l = c[i][j];
                int size = (l == L1) ? 1 : (l == L2) ? 2 : cnt[l];
                sb[q++] = (char) ('0' + size);
            }
            String str = new String(sb);
            if (targetsByInv != null) {
                List<String[]> cand = targetsByInv.get(invariant(str, n));
                if (cand != null) {
                    String cf = canonComplete(str, n);
                    for (String[] tg : cand) if (tg[1].equals(cf)) targetHits.add(tg[0]);
                }
            }
            if (found.size() < FOUND_CAP) found.add(str); else foundOverflow.set(true);
        }

        long canon7() { // canonical form of the 7-point prefix (all high), min over 7! perms of renamed string
            long best = Long.MAX_VALUE;
            for (int[] pi : PERM7) {
                int[] ren = new int[NL]; Arrays.fill(ren, -1); int nx = 0; long cur = 0;
                for (int i = 0; i < 7; i++) for (int j = i + 1; j < 7; j++) {
                    int l = c[pi[i]][pi[j]]; int r = ren[l]; if (r < 0) { r = nx++; ren[l] = r; }
                    cur = (cur << 3) | r;
                }
                if (cur < best) best = cur;
            }
            return best;
        }
    }

    static int[][] PERM7;

    // ================= local tests =================
    /** Configurations containing the newly assigned pair (u,v) only. */
    static boolean localInc(int[][] c, int n, int u, int v) {
        int L = c[u][v];
        // Lemma 2.1 (i)
        int cu = 0, cv = 0;
        for (int x = 0; x < n; x++) { if (c[u][x] == L) cu++; if (c[v][x] == L) cv++; }
        if (cu > 3 || cv > 3) return false;
        // Lemma 2.1 (ii),(iii): pairs {a,b} with third point x, where {a,x} = {u,v}
        for (int side = 0; side < 2; side++) {
            int a = side == 0 ? u : v, x = side == 0 ? v : u;
            for (int b = 0; b < n; b++) {
                if (b == a || b == x) continue;
                int kb = c[b][x];
                if (kb < 0) continue;
                int ka = c[a][x];
                int same = 0, key = 0;
                for (int y = 0; y < n; y++) {
                    if (y == a || y == b) continue;
                    int p = c[a][y], q = c[b][y];
                    if (p < 0 || q < 0) continue;
                    if (p == ka && q == kb) key++;
                    if (p == q) same++;
                }
                if (key > 2) return false;
                if (ka == kb && same > 2) return false;
            }
        }
        // Lemma 2.2 on 4-sets {u,v,a,b}
        for (int a = 0; a < n; a++) {
            if (a == u || a == v) continue;
            for (int b = a + 1; b < n; b++) {
                if (b == u || b == v) continue;
                int ab = c[a][b], ua = c[u][a], vb = c[v][b], ub = c[u][b], va = c[v][a];
                if (ab < 0 || ua < 0 || vb < 0 || ub < 0 || va < 0) continue;
                if (L == ab && ua == vb && ub == va) return false;
            }
        }
        return true;
    }

    /** Full (non-incremental) local tests on a partial pattern; independent code path used for checking. */
    static boolean fullLocal(int[][] c, int n) {
        for (int a = 0; a < n; a++) {
            int[] cnt = new int[NL];
            for (int x = 0; x < n; x++) if (x != a && c[a][x] >= 0 && ++cnt[c[a][x]] > 3) return false;
        }
        for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) {
            int[] cnt = new int[NL * NL];
            int same = 0;
            for (int x = 0; x < n; x++) {
                if (x == a || x == b) continue;
                int p = c[a][x], q = c[b][x];
                if (p < 0 || q < 0) continue;
                if (++cnt[p * NL + q] > 2) return false;
                if (p == q && ++same > 2) return false;
            }
        }
        for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) for (int d = b + 1; d < n; d++) for (int e = d + 1; e < n; e++) {
            int ab = c[a][b], de = c[d][e], ad = c[a][d], be = c[b][e], ae = c[a][e], bd = c[b][d];
            if (ab < 0 || de < 0 || ad < 0 || be < 0 || ae < 0 || bd < 0) continue;
            if (ab == de && ad == be && ae == bd) return false;
        }
        return true;
    }

    // ================= permutations / U handling =================
    static int[][] perms(int k) {
        List<int[]> out = new ArrayList<>();
        permRec(new int[k], new boolean[k], 0, out);
        return out.toArray(new int[0][]);
    }
    static void permRec(int[] cur, boolean[] used, int pos, List<int[]> out) {
        if (pos == cur.length) { out.add(cur.clone()); return; }
        for (int v = 0; v < cur.length; v++) if (!used[v]) { used[v] = true; cur[pos] = v; permRec(cur, used, pos + 1, out); used[v] = false; }
    }

    /** renamed packed string of lab after permuting points by pi (new i = old pi[i]). */
    static long permuted(int[] lab, int[] pi, int[][] idx, int k) {
        int[] ren = new int[NL]; Arrays.fill(ren, -1); int nx = 0; long cur = 0;
        for (int i = 0; i < k; i++) for (int j = i + 1; j < k; j++) {
            int l = lab[idx[pi[i]][pi[j]]]; int r = ren[l]; if (r < 0) { r = nx++; ren[l] = r; }
            cur = (cur << 3) | r;
        }
        return cur;
    }
    static long canonK(int[] lab, int[][] perms, int[][] idx, int k) {
        long best = Long.MAX_VALUE;
        for (int[] pi : perms) { long x = permuted(lab, pi, idx, k); if (x < best) best = x; }
        return best;
    }
    static int[] unpack(long code, int np) {
        int[] lab = new int[np];
        for (int p = np - 1; p >= 0; p--) { lab[p] = (int) (code & 7); code >>>= 3; }
        return lab;
    }
    static String toStr(long code, int np) {
        StringBuilder sb = new StringBuilder();
        for (int l : unpack(code, np)) sb.append((char) ('1' + l));
        return sb.toString();
    }

    /** Open-addressing set of non-negative longs. */
    static final class LongSet {
        long[] tab; int size; int mask;
        LongSet(int cap) { int c = 1; while (c < cap * 2) c <<= 1; tab = new long[c]; Arrays.fill(tab, -1L); mask = c - 1; }
        static int h(long x) { x ^= x >>> 33; x *= 0xff51afd7ed558ccdL; x ^= x >>> 33; return (int) x; }
        boolean add(long x) {
            if (size * 2 >= tab.length) grow();
            int i = h(x) & mask;
            while (tab[i] != -1L) { if (tab[i] == x) return false; i = (i + 1) & mask; }
            tab[i] = x; size++; return true;
        }
        boolean contains(long x) {
            int i = h(x) & mask;
            while (tab[i] != -1L) { if (tab[i] == x) return true; i = (i + 1) & mask; }
            return false;
        }
        void grow() {
            long[] old = tab; tab = new long[old.length * 2]; Arrays.fill(tab, -1L); mask = tab.length - 1; size = 0;
            for (long x : old) if (x != -1L) add(x);
        }
    }

    // ================= main =================
    public static void main(String[] args) throws Exception {
        Map<String, String> opt = new HashMap<>();
        Set<String> flags = new HashSet<>();
        for (int i = 0; i < args.length; i++) {
            String a = args[i];
            if (a.equals("--nolocal") || a.equals("--paranoid") || a.equals("--nolook")) flags.add(a);
            else opt.put(a, args[++i]);
        }
        int n = Integer.parseInt(opt.getOrDefault("--n", "9"));
        int k = n - 3, H = n - 3, kp = k * (k - 1) / 2;
        int threads = Integer.parseInt(opt.getOrDefault("--threads", "1"));
        if (threads > 2) { System.err.println("CPU cap: threads clamped to 2"); threads = 2; }
        String uFile = opt.get("--U");
        String mode = opt.getOrDefault("--mode", "auto");
        boolean useLocal = !flags.contains("--nolocal"), paranoid = flags.contains("--paranoid");
        String level7File = opt.get("--level7");

        // ---- load U ----
        int[][] permK = perms(k);
        int[][] idxK = new int[k][k];
        { int p = 0; for (int i = 0; i < k; i++) for (int j = i + 1; j < k; j++) { idxK[i][j] = idxK[j][i] = p++; } }
        List<Long> canonU = new ArrayList<>();
        Set<Long> canonSet = new HashSet<>();
        int skippedFmt = 0, skippedLabels = 0, dupes = 0, totalLines = 0;
        List<String> toks = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new InputStreamReader(
                uFile.endsWith(".gz") ? new GZIPInputStream(new FileInputStream(uFile)) : new FileInputStream(uFile), "US-ASCII"))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) continue;
                totalLines++;
                toks.add(line.split("[\\s,;]+")[0]);
            }
        }
        // optional random subset of the file's lines (stand-in U), chosen before canonicalisation
        if (opt.containsKey("--subset")) {
            int m = Integer.parseInt(opt.get("--subset"));
            Collections.shuffle(toks, new Random(Long.parseLong(opt.getOrDefault("--rseed", "1"))));
            if (m < toks.size()) toks = new ArrayList<>(toks.subList(0, m));
        }
        for (String tok : toks) {
            if (tok.length() != kp || !tok.matches("[1-9]+")) { skippedFmt++; continue; }
            int[] lab = new int[kp];
            Map<Character, Integer> ren = new HashMap<>();
            for (int p = 0; p < kp; p++) { char ch = tok.charAt(p); Integer r = ren.get(ch); if (r == null) { r = ren.size(); ren.put(ch, r); } lab[p] = r; }
            if (ren.size() > H) { skippedLabels++; continue; }
            long cf = canonK(lab, permK, idxK, k); // re-canonicalise with OUR canonical form
            if (!canonSet.add(cf)) { dupes++; continue; }
            canonU.add(cf);
        }
        Collections.sort(canonU);
        System.out.printf("U: file=%s lines=%d kept=%d (skipped: bad format %d, >%d labels %d, duplicate classes %d)%n",
                uFile, totalLines, canonU.size(), skippedFmt, H, skippedLabels, dupes);

        // ---- orbits (seeds + membership) ----
        long est = (long) canonU.size() * permK.length;
        if (mode.equals("auto")) mode = est <= 40_000_000L ? "orbit" : "canon";
        final List<long[]> orbitsList = new ArrayList<>();
        LongSet orbitSet = new LongSet(1024);
        long t0 = System.nanoTime();
        for (long cf : canonU) {
            int[] lab = unpack(cf, kp);
            LongSet one = new LongSet(64);
            List<Long> el = new ArrayList<>();
            for (int[] pi : permK) { long x = permuted(lab, pi, idxK, k); if (one.add(x)) el.add(x); }
            long[] arr = new long[el.size()]; for (int i = 0; i < arr.length; i++) arr[i] = el.get(i);
            Arrays.sort(arr);
            orbitsList.add(arr);
            if (mode.equals("orbit")) for (long x : arr) orbitSet.add(x);
        }
        long totalSeeds = 0; for (long[] o : orbitsList) totalSeeds += o.length;
        System.out.printf("mode=%s, labelled seeds (orbit elements) = %d, built in %.1fs%n", mode, totalSeeds, (System.nanoTime() - t0) / 1e9);

        final Set<Long> canonSetF = canonSet;
        final String modeF = mode;
        final LongSet orbitSetF = orbitSet;
        // canon mode: per-thread cache
        final ThreadLocal<HashMap<Long, Boolean>> cache = ThreadLocal.withInitial(HashMap::new);
        UMember umem = (packed, lab) -> {
            if (modeF.equals("orbit")) return orbitSetF.contains(packed);
            HashMap<Long, Boolean> cch = cache.get();
            Boolean b = cch.get(packed);
            if (b == null) { b = canonSetF.contains(canonK(lab, permK, idxK, k)); if (cch.size() < 5_000_000) cch.put(packed, b); }
            return b;
        };

        if (level7File != null) PERM7 = perms(7);
        if (opt.containsKey("--targets")) {
            targetsByInv = new HashMap<>();
            int nt = 0;
            for (String line : java.nio.file.Files.readAllLines(java.nio.file.Paths.get(opt.get("--targets")))) {
                line = line.trim(); if (line.isEmpty() || line.startsWith("#")) continue;
                String[] f = line.split("[ \\t]+");
                if (f[1].length() != n * (n - 1) / 2) continue;
                String cf = canonComplete(f[1], n);
                targetsByInv.computeIfAbsent(invariant(cf, n), z -> new ArrayList<>()).add(new String[]{f[0], cf});
                nt++;
            }
            System.out.printf("targets: %d loaded for n=%d%n", nt, n);
        }
        Set<String> allHits = new TreeSet<>();

        // ---- arrangements ----
        List<Integer> arrs = new ArrayList<>();
        String arrOpt = opt.getOrDefault("--arr", "all");
        if (level7File != null) arrOpt = "star";
        for (int a = 0; a < ARR.length; a++) if (arrOpt.equals("all") || Arrays.asList(arrOpt.split(",")).contains(ARR_NAMES[a])) arrs.add(a);

        PrintWriter out = opt.containsKey("--out") ? new PrintWriter(new FileWriter(opt.get("--out"))) : null;
        long grandLeaves = 0;
        Set<String> allCanon = new TreeSet<>();
        for (int a : arrs) {
            long ts = System.nanoTime();
            CoverageSearch cs = new CoverageSearch(n, a, umem, permK, useLocal, paranoid, level7File != null, level7File != null ? 6 : -1, !flags.contains("--nolook"));
            ExecutorService ex = Executors.newFixedThreadPool(threads);
            List<Future<?>> fs = new ArrayList<>();
            AtomicInteger done = new AtomicInteger();
            final int total = orbitsList.size();
            for (long[] orbit : orbitsList) {
                fs.add(ex.submit(() -> {
                    State st = cs.new State();
                    for (long seed : orbit) {
                        int[] lab = unpack(seed, kp);
                        st.loadSeed(lab);
                        if (useLocal && !fullLocal(st.c, n)) continue; // U members should pass; guard anyway
                        cs.seeds.incrementAndGet();
                        cs.levelCount.incrementAndGet(k);
                        st.dfs(0);
                        cs.nodes.addAndGet(st.myNodes); st.myNodes = 0;
                        cs.leaves.addAndGet(st.myLeaves); st.myLeaves = 0;
                        long now = System.nanoTime();
                        if (now - cs.lastPrint.get() > 30_000_000_000L) { cs.lastPrint.set(now);
                            System.err.printf("  [%s] %d/%d U patterns done, seeds=%d nodes=%d leaves=%d levels=%s%n", ARR_NAMES[a], done.get(), total, cs.seeds.get(), cs.nodes.get(), cs.leaves.get(), cs.levelCount.toString()); }
                    }
                    int d = done.incrementAndGet();
                    if (d % 20000 == 0) System.err.printf("  [%s] %d/%d U patterns, nodes=%d leaves=%d%n", ARR_NAMES[a], d, total, cs.nodes.get(), cs.leaves.get());
                    return null;
                }));
            }
            for (Future<?> f : fs) f.get();
            ex.shutdown();
            double sec = (System.nanoTime() - ts) / 1e9;
            Set<String> canonFound = new TreeSet<>();
            for (String s : cs.found) canonFound.add(canonComplete(s, n));
            allCanon.addAll(canonFound);
            StringBuilder lv = new StringBuilder();
            for (int v = k; v <= n; v++) lv.append(v).append(':').append(cs.levelCount.get(v)).append(' ');
            System.out.printf("arr=%-17s seeds=%d nodes=%d levels[%s] leaves=%d distinctCanonical=%d time=%.1fs%n",
                    ARR_NAMES[a], cs.seeds.get(), cs.nodes.get(), lv.toString().trim(), cs.leaves.get(), canonFound.size(), sec);
            grandLeaves += cs.leaves.get();
            if (targetsByInv != null) { System.out.printf("  target hits in %s: %s%n", ARR_NAMES[a], new TreeSet<>(cs.targetHits)); allHits.addAll(cs.targetHits); }
            if (cs.foundOverflow.get()) System.out.printf("  NOTE: more than %d labelled leaves; stored set truncated (canonical count is a lower bound)%n", FOUND_CAP);
            if (out != null) for (String s : canonFound) out.println(ARR_NAMES[a] + " " + s);
            if (level7File != null) {
                List<Long> l7 = new ArrayList<>(cs.level7); Collections.sort(l7);
                try (PrintWriter w = new PrintWriter(new FileWriter(level7File))) { for (long x : l7) w.println(toStr(x, 21)); }
                System.out.printf("level7: %d distinct canonical 7-point all-high patterns (local+capacity+all 6-subsets in U) -> %s%n", l7.size(), level7File);
            }
        }
        if (out != null) out.close();
        System.out.printf("TOTAL leaves=%d distinctCanonical(all arrangements)=%d%n", grandLeaves, allCanon.size());
        if (targetsByInv != null) System.out.printf("TARGETS HIT: %d %s%n", allHits.size(), allHits);
    }

    /** Permutation-invariant: sorted list of each point's sorted class-size multiset. */
    static String invariant(String s, int n) {
        int[][] c = new int[n][n];
        int p = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { c[i][j] = c[j][i] = s.charAt(p++) - '0'; }
        String[] rows = new String[n];
        for (int i = 0; i < n; i++) {
            char[] r = new char[n - 1]; int q = 0;
            for (int j = 0; j < n; j++) if (j != i) r[q++] = (char) ('0' + c[i][j]);
            Arrays.sort(r); rows[i] = new String(r);
        }
        Arrays.sort(rows);
        return String.join("|", rows);
    }

    /** Canonical form of a complete pattern given as class-size digits (unique names): min over all n! perms. */
    static String canonComplete(String s, int n) {
        int[][] c = new int[n][n];
        int p = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { c[i][j] = c[j][i] = s.charAt(p++) - '0'; }
        char[] best = null, cur = new char[s.length()];
        int[] pi = new int[n];
        for (int i = 0; i < n; i++) pi[i] = i;
        // iterate all permutations (Heap's algorithm is not lexicographic; we just need the min)
        int[] cc = new int[n];
        best = s.toCharArray();
        while (true) {
            int q = 0; boolean less = false, greater = false;
            for (int i = 0; i < n && !greater; i++) for (int j = i + 1; j < n; j++) {
                char ch = (char) ('0' + c[pi[i]][pi[j]]);
                cur[q] = ch;
                if (!less) { if (ch > best[q]) { greater = true; break; } if (ch < best[q]) less = true; }
                q++;
            }
            if (less && !greater) best = cur.clone();
            // next permutation (lexicographic)
            int i = n - 2; while (i >= 0 && pi[i] >= pi[i + 1]) i--;
            if (i < 0) break;
            int j = n - 1; while (pi[j] <= pi[i]) j--;
            int tmp = pi[i]; pi[i] = pi[j]; pi[j] = tmp;
            for (int x = i + 1, y = n - 1; x < y; x++, y--) { tmp = pi[x]; pi[x] = pi[y]; pi[y] = tmp; }
        }
        return new String(best);
    }
}
