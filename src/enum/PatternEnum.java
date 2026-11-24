import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import java.util.zip.*;

/**
 * Builder E, Phase 2: enumerate distance patterns on k = 3..6 points with at most 6 labels,
 * up to (S_k on points) x (label renaming), keeping only patterns passing the local tests
 * (Wang Lemma 2.1 (i),(ii),(iii)) and the four-point exclusion (Lemma 2.2).
 *
 * Encoding: a pattern on k points is a label array over pairs in lex order
 * (1,2),(1,3),...,(1,k),(2,3),...,(k-1,k); labels 0..5, renamed by first occurrence.
 * Packed into a long with 3 bits per pair, first pair most significant, so numeric order
 * on packed longs == lexicographic order on label strings of equal length.
 *
 * Canonical form (our own): the minimum, over ALL k! point permutations, of the
 * first-occurrence-renamed label string. No invariant-based shortcut (k! <= 720 is cheap),
 * so canonical form correctness does not depend on any invariant being correct.
 *
 * Two independent generation methods:
 *   A ("extend"): extend each canonical (k-1)-point rep by point k with every label choice
 *      (existing labels or one new label, at most 6 total), test, canonicalise, dedupe.
 *   B ("orderly"): direct depth-first search over ALL restricted-growth label strings on the
 *      C(k,2) pairs (pairs assigned in colex order so tests prune early), keep a full
 *      labelling iff it passes the tests and its lex string equals its canonical form.
 * Both must give the same set.
 *
 * Usage: java PatternEnum A|B outdir threads
 */
public class PatternEnum {
    static final int MAXL = 6;

    // ---------- pair indexing ----------
    static int[][] pairIndexLex(int k) {           // lex order index for i<j (0-based)
        int[][] idx = new int[k][k];
        int p = 0;
        for (int i = 0; i < k; i++) for (int j = i + 1; j < k; j++) { idx[i][j] = p; idx[j][i] = p; p++; }
        for (int i = 0; i < k; i++) idx[i][i] = -1;
        return idx;
    }

    // ---------- permutations ----------
    static List<int[]> perms(int k) {
        List<int[]> out = new ArrayList<>();
        permRec(new int[k], new boolean[k], 0, k, out);
        return out;
    }
    static void permRec(int[] cur, boolean[] used, int pos, int k, List<int[]> out) {
        if (pos == k) { out.add(cur.clone()); return; }
        for (int v = 0; v < k; v++) if (!used[v]) { used[v] = true; cur[pos] = v; permRec(cur, used, pos + 1, k, out); used[v] = false; }
    }

    /** For each permutation pi (new point i = old point pi[i]) and each new lex pair index,
     *  the old lex pair index. */
    static int[][] permPairMaps(int k) {
        int[][] idx = pairIndexLex(k);
        List<int[]> ps = perms(k);
        int np = k * (k - 1) / 2;
        int[][] maps = new int[ps.size()][np];
        for (int t = 0; t < ps.size(); t++) {
            int[] pi = ps.get(t);
            int p = 0;
            for (int i = 0; i < k; i++) for (int j = i + 1; j < k; j++) maps[t][p++] = idx[pi[i]][pi[j]];
        }
        return maps;
    }

    // ---------- canonical form ----------
    static final class Canon {
        final int k, np;
        final int[][] maps;
        Canon(int k) { this.k = k; this.np = k * (k - 1) / 2; this.maps = permPairMaps(k); }

        /** Canonical packed string: min over all perms of first-occurrence-renamed string. */
        long canon(int[] lab) {
            long best = Long.MAX_VALUE;
            int[] ren = new int[MAXL];
            for (int[] m : maps) {
                Arrays.fill(ren, -1);
                int next = 0;
                long cur = 0;
                // compare progressively against best (best stored as packed; extract digit)
                boolean greater = false, less = (best == Long.MAX_VALUE);
                for (int p = 0; p < np; p++) {
                    int l = lab[m[p]];
                    int r = ren[l];
                    if (r < 0) { r = next++; ren[l] = r; }
                    cur = (cur << 3) | r;
                    if (!less) {
                        int bd = (int) ((best >>> (3 * (np - 1 - p))) & 7);
                        if (r > bd) { greater = true; break; }
                        if (r < bd) less = true;
                    }
                }
                if (!greater && less) best = cur;
            }
            return best;
        }
    }

    static int[] unpack(long code, int np) {
        int[] lab = new int[np];
        for (int p = np - 1; p >= 0; p--) { lab[p] = (int) (code & 7); code >>>= 3; }
        return lab;
    }
    static long pack(int[] lab) { long c = 0; for (int l : lab) c = (c << 3) | l; return c; }
    static String toStr(long code, int np) {
        StringBuilder sb = new StringBuilder();
        for (int l : unpack(code, np)) sb.append((char) ('1' + l));
        return sb.toString();
    }

    // ---------- tests on a (possibly partial) labelling given as matrix, -1 = unassigned ----------
    /** Returns true iff the assigned pairs violate none of: Lemma 2.1 (i),(ii),(iii), Lemma 2.2. */
    static boolean passes(int k, int[][] c) {
        // (i): at each point, at most 3 incident pairs with the same label
        for (int a = 0; a < k; a++) {
            int[] cnt = new int[MAXL];
            for (int x = 0; x < k; x++) if (x != a && c[a][x] >= 0) if (++cnt[c[a][x]] > 3) return false;
        }
        // (ii) and (iii)
        for (int a = 0; a < k; a++) for (int b = a + 1; b < k; b++) {
            int same = 0;
            int[] cnt = new int[MAXL * MAXL];
            for (int x = 0; x < k; x++) {
                if (x == a || x == b) continue;
                int u = c[a][x], v = c[b][x];
                if (u < 0 || v < 0) continue;
                // (ii): ordered label pair (u,v) for the ordered point pair (a,b); includes u==v
                if (++cnt[u * MAXL + v] > 2) return false;
                // (iii)
                if (u == v && ++same > 2) return false;
            }
        }
        // Lemma 2.2: all three perfect matchings of {a,b,c,d} have equal labels on both edges
        for (int a = 0; a < k; a++) for (int b = a + 1; b < k; b++) for (int cc = b + 1; cc < k; cc++) for (int d = cc + 1; d < k; d++) {
            int ab = c[a][b], cd = c[cc][d], ac = c[a][cc], bd = c[b][d], ad = c[a][d], bc = c[b][cc];
            if (ab < 0 || cd < 0 || ac < 0 || bd < 0 || ad < 0 || bc < 0) continue;
            if (ab == cd && ac == bd && ad == bc) return false;
        }
        return true;
    }

    static int[][] toMatrix(int k, int[] lab) {
        int[][] c = new int[k][k];
        int p = 0;
        for (int i = 0; i < k; i++) { c[i][i] = -1; }
        for (int i = 0; i < k; i++) for (int j = i + 1; j < k; j++) { c[i][j] = lab[p]; c[j][i] = lab[p]; p++; }
        return c;
    }

    static int numLabels(int[] lab) { int m = -1; for (int l : lab) m = Math.max(m, l); return m + 1; }

    // ---------- Method A: extension ----------
    static Set<Long> extend(Set<Long> prev, int k, int threads) throws Exception {
        int npPrev = (k - 1) * (k - 2) / 2;
        int np = k * (k - 1) / 2;
        Canon canon = new Canon(k);
        int[][] idxK = pairIndexLex(k);
        int[][] idxP = pairIndexLex(k - 1);
        Set<Long> out = ConcurrentHashMap.newKeySet();
        List<Long> reps = new ArrayList<>(prev);
        Collections.sort(reps);
        ExecutorService ex = Executors.newFixedThreadPool(threads);
        AtomicLong tested = new AtomicLong(), passed = new AtomicLong();
        AtomicInteger done = new AtomicInteger();
        List<Future<?>> fs = new ArrayList<>();
        for (long rep : reps) {
            fs.add(ex.submit(() -> {
                int[] old = unpack(rep, npPrev);
                int t0 = numLabels(old);
                int[][] c = new int[k][k];
                for (int[] row : c) Arrays.fill(row, -1);
                for (int i = 0; i < k - 1; i++) for (int j = i + 1; j < k - 1; j++) { c[i][j] = c[j][i] = old[idxP[i][j]]; }
                int[] lab = new int[np];
                long[] cnt = new long[2];
                extRec(k, c, 0, t0, lab, idxK, canon, out, cnt);
                tested.addAndGet(cnt[0]); passed.addAndGet(cnt[1]);
                int d = done.incrementAndGet();
                if (d % 50 == 0) System.err.printf("  k=%d: %d/%d reps, %d distinct%n", k, d, reps.size(), out.size());
            }));
        }
        for (Future<?> f : fs) f.get();
        ex.shutdown();
        System.err.printf("  k=%d: candidates completed=%d passing=%d distinct=%d%n", k, tested.get(), passed.get(), out.size());
        return out;
    }

    static void extRec(int k, int[][] c, int i, int t, int[] lab, int[][] idxK, Canon canon, Set<Long> out, long[] cnt) {
        int n = k - 1; // new point index
        if (i == n) {
            cnt[0]++;
            if (!passes(k, c)) return;
            cnt[1]++;
            for (int a = 0; a < k; a++) for (int b = a + 1; b < k; b++) lab[idxK[a][b]] = c[a][b];
            out.add(canon.canon(lab));
            return;
        }
        int maxL = Math.min(t, MAXL - 1); // labels 0..t-1 existing, t new (if t < MAXL)
        for (int l = 0; l <= maxL; l++) {
            c[i][n] = c[n][i] = l;
            extRec(k, c, i + 1, (l == t) ? t + 1 : t, lab, idxK, canon, out, cnt);
        }
        c[i][n] = c[n][i] = -1;
    }

    // ---------- Method B: orderly brute force over restricted growth strings ----------
    static long[] orderly(int k, int threads, Set<Long> sink) throws Exception {
        int np = k * (k - 1) / 2;
        // colex pair order: (0,1),(0,2),(1,2),(0,3),(1,3),(2,3),...
        int[][] order = new int[np][2];
        int p = 0;
        for (int j = 1; j < k; j++) for (int i = 0; i < j; i++) { order[p][0] = i; order[p][1] = j; p++; }
        int[][] idxK = pairIndexLex(k);
        Canon canon = new Canon(k);
        // split work on the labels of the first few pairs
        int split = Math.min(np, 6);
        List<int[]> prefixes = new ArrayList<>();
        prefRec(new int[split], 0, 0, split, prefixes);
        ExecutorService ex = Executors.newFixedThreadPool(threads);
        AtomicLong full = new AtomicLong(), canonCount = new AtomicLong();
        List<Future<?>> fs = new ArrayList<>();
        for (int[] pre : prefixes) {
            fs.add(ex.submit(() -> {
                int[][] c = new int[k][k];
                for (int[] row : c) Arrays.fill(row, -1);
                int t = 0;
                for (int q = 0; q < pre.length; q++) {
                    c[order[q][0]][order[q][1]] = c[order[q][1]][order[q][0]] = pre[q];
                    if (pre[q] == t) t++;
                    // test after each completed point (when pair (j-1,j) assigned) -- cheap enough to test every step
                    if (!passes(k, c)) return;
                }
                long[] cnt = new long[2];
                int[] lab = new int[np];
                orderlyRec(k, c, order, pre.length, t, idxK, canon, lab, cnt, sink);
                full.addAndGet(cnt[0]); canonCount.addAndGet(cnt[1]);
            }));
        }
        for (Future<?> f : fs) f.get();
        ex.shutdown();
        return new long[]{full.get(), canonCount.get()};
    }

    static void prefRec(int[] cur, int pos, int t, int len, List<int[]> out) {
        if (pos == len) { out.add(cur.clone()); return; }
        for (int l = 0; l <= Math.min(t, MAXL - 1); l++) { cur[pos] = l; prefRec(cur, pos + 1, l == t ? t + 1 : t, len, out); }
    }

    static void orderlyRec(int k, int[][] c, int[][] order, int pos, int t, int[][] idxK, Canon canon, int[] lab, long[] cnt, Set<Long> sink) {
        int np = order.length;
        if (pos == np) {
            cnt[0]++;
            for (int a = 0; a < k; a++) for (int b = a + 1; b < k; b++) lab[idxK[a][b]] = c[a][b];
            // relabel lex string by first occurrence, compare with canonical form
            int[] ren = new int[MAXL]; Arrays.fill(ren, -1); int nx = 0;
            int[] rl = new int[np];
            for (int q = 0; q < np; q++) { int l = lab[q]; if (ren[l] < 0) ren[l] = nx++; rl[q] = ren[l]; }
            long self = pack(rl);
            if (canon.canon(lab) == self) { cnt[1]++; if (sink != null) sink.add(self); }
            return;
        }
        int a = order[pos][0], b = order[pos][1];
        for (int l = 0; l <= Math.min(t, MAXL - 1); l++) {
            c[a][b] = c[b][a] = l;
            if (passes(k, c)) orderlyRec(k, c, order, pos + 1, l == t ? t + 1 : t, idxK, canon, lab, cnt, sink);
        }
        c[a][b] = c[b][a] = -1;
    }

    // ---------- main ----------
    public static void main(String[] args) throws Exception {
        String mode = args.length > 0 ? args[0] : "A";
        String outdir = args.length > 1 ? args[1] : ".";
        int threads = args.length > 2 ? Integer.parseInt(args[2]) : 4;
        int kmax = args.length > 3 ? Integer.parseInt(args[3]) : 6;
        if (mode.equals("A")) {
            Set<Long> cur = new HashSet<>();
            cur.add(0L); // 2-point pattern, single label
            for (int k = 3; k <= kmax; k++) {
                long t0 = System.nanoTime();
                cur = extend(cur, k, threads);
                double sec = (System.nanoTime() - t0) / 1e9;
                int np = k * (k - 1) / 2;
                int[] byR = new int[MAXL + 1];
                for (long code : cur) byR[numLabels(unpack(code, np))]++;
                System.out.printf("A k=%d count=%d byR=%s time=%.1fs%n", k, cur.size(), Arrays.toString(byR), sec);
                List<Long> sorted = new ArrayList<>(cur);
                Collections.sort(sorted);
                File f = new File(outdir, "patterns" + k + "_local" + (k == 6 ? ".txt.gz" : ".txt"));
                OutputStream os = new FileOutputStream(f);
                if (k == 6) os = new GZIPOutputStream(os);
                try (Writer w = new BufferedWriter(new OutputStreamWriter(os, "US-ASCII"))) {
                    for (long code : sorted) { w.write(toStr(code, np)); w.write('\n'); }
                }
            }
        } else {
            for (int k = 3; k <= kmax; k++) {
                long t0 = System.nanoTime();
                Set<Long> sink = ConcurrentHashMap.newKeySet();
                long[] r = orderly(k, threads, sink);
                double sec = (System.nanoTime() - t0) / 1e9;
                int np = k * (k - 1) / 2;
                int[] byR = new int[MAXL + 1];
                for (long code : sink) byR[numLabels(unpack(code, np))]++;
                System.out.printf("B k=%d labelledPassing=%d canonical=%d distinct=%d byR=%s time=%.1fs%n", k, r[0], r[1], sink.size(), Arrays.toString(byR), sec);
                List<Long> sorted = new ArrayList<>(sink);
                Collections.sort(sorted);
                File f = new File(outdir, "patterns" + k + "_local_B.txt.gz");
                try (Writer w = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(f)), "US-ASCII"))) {
                    for (long code : sorted) { w.write(toStr(code, np)); w.write('\n'); }
                }
            }
        }
    }
}
