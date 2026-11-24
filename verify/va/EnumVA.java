// V-a independent enumerator of 6-point distance patterns (Wang Lemma 2.1 (i)-(iii), Lemma 2.2), <=6 labels.
// Written from the TeX spec only. Single-threaded. Canonical form: min over all k! point orderings of
// the first-occurrence-renamed label string in lex pair order (1,2),(1,3),...,(k-1,k).
import java.io.*;
import java.util.*;

public class EnumVA {
    static int K;
    static int[][] pidx = new int[7][7];
    static int[][] perms;

    static void setK(int k) {
        K = k; int c = 0;
        for (int i = 0; i < k; i++) for (int j = i + 1; j < k; j++) { pidx[i][j] = c; pidx[j][i] = c; c++; }
        List<int[]> ps = new ArrayList<>();
        permute(new int[k], new boolean[k], 0, ps);
        perms = ps.toArray(new int[0][]);
    }
    static void permute(int[] a, boolean[] used, int d, List<int[]> out) {
        if (d == a.length) { out.add(a.clone()); return; }
        for (int i = 0; i < a.length; i++) if (!used[i]) { used[i] = true; a[d] = i; permute(a, used, d + 1, out); used[i] = false; }
    }
    static int lab(int[] p, int i, int j) { return p[pidx[i][j]]; }

    // 0 = unassigned. Returns true if no violation among assigned pairs.
    static boolean ok(int[] p, int k) {
        // L(i)
        for (int a = 0; a < k; a++) {
            int[] cnt = new int[8];
            for (int x = 0; x < k; x++) if (x != a) { int l = lab(p, a, x); if (l > 0 && ++cnt[l] > 3) return false; }
        }
        for (int a = 0; a < k; a++) for (int b = a + 1; b < k; b++) {
            int same = 0; int[][] uv = new int[8][8];
            for (int x = 0; x < k; x++) {
                if (x == a || x == b) continue;
                int u = lab(p, a, x), v = lab(p, b, x);
                if (u == 0 || v == 0) continue;
                if (++uv[u][v] > 2) return false;          // L(ii)
                if (u == v && ++same > 2) return false;      // L(iii)
            }
        }
        // Lemma 2.2: all three perfect matchings of a 4-set equal-labelled
        for (int a = 0; a < k; a++) for (int b = a + 1; b < k; b++) for (int c = b + 1; c < k; c++) for (int d = c + 1; d < k; d++) {
            int ab = lab(p, a, b), cd = lab(p, c, d), ac = lab(p, a, c), bd = lab(p, b, d), ad = lab(p, a, d), bc = lab(p, b, c);
            if (ab > 0 && ab == cd && ac > 0 && ac == bd && ad > 0 && ad == bc) return false;
        }
        return true;
    }

    static String canon(int[] p, int k) {
        int n = k * (k - 1) / 2;
        int[] best = null; int[] cur = new int[n]; int[] map = new int[8];
        for (int[] pi : perms) {
            Arrays.fill(map, 0); int next = 1; boolean worse = false, better = (best == null);
            int c = 0;
            for (int i = 0; i < k && !worse; i++) for (int j = i + 1; j < k; j++) {
                int l = lab(p, pi[i], pi[j]);
                if (map[l] == 0) map[l] = next++;
                int v = map[l]; cur[c] = v;
                if (!better) { if (v < best[c]) better = true; else if (v > best[c]) { worse = true; break; } }
                c++;
            }
            if (!worse && better) { best = cur.clone(); }
        }
        StringBuilder sb = new StringBuilder();
        for (int v : best) sb.append((char) ('0' + v));
        return sb.toString();
    }

    static TreeSet<String> out;
    static void extend(int[] p, int k, int t, int idx) {
        // assign pair (idx, k-1) for idx = 0..k-2
        if (idx == k - 1) { out.add(canon(p, k)); return; }
        int pi = pidx[idx][k - 1];
        int maxl = Math.min(t + 1, 6);
        for (int l = 1; l <= maxl; l++) {
            p[pi] = l;
            if (ok(p, k)) extend(p, k, Math.max(t, l), idx + 1);
            p[pi] = 0;
        }
    }

    public static void main(String[] args) throws Exception {
        String outdir = args[0];
        TreeSet<String> prev = new TreeSet<>(); prev.add("1");
        for (int k = 3; k <= 6; k++) {
            setK(k);
            out = new TreeSet<>();
            int nOld = (k - 1) * (k - 2) / 2;
            setK(k - 1); int[][] oldIdx = new int[7][7]; for (int i=0;i<k-1;i++) for(int j=0;j<k-1;j++) oldIdx[i][j]=pidx[i][j];
            setK(k);
            for (String s : prev) {
                int[] p = new int[k * (k - 1) / 2]; int t = 0;
                for (int i = 0; i < k - 1; i++) for (int j = i + 1; j < k - 1; j++) { int v = s.charAt(oldIdx[i][j]) - '0'; p[pidx[i][j]] = v; t = Math.max(t, v); }
                extend(p, k, t, 0);
            }
            int[] byR = new int[16];
            try (PrintWriter w = new PrintWriter(new FileWriter(outdir + "/va_patterns" + k + "_local.txt"))) {
                for (String s : out) { w.println(s); int r = 0; for (char ch : s.toCharArray()) r = Math.max(r, ch - '0'); byR[r]++; }
            }
            System.out.println("k=" + k + " count=" + out.size() + " byR=" + Arrays.toString(Arrays.copyOfRange(byR, 1, 7)));
            prev = out;
        }
    }
}
