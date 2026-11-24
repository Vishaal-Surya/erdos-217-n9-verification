
## Phase 2/3 — V-a (independent re-implementation, red team)

**Independence:** written only from PROMPT.md, NOTES/*.md (specs/claims) and Wang's TeX (`raw/2611.03688/src/crescent9_v4__1_.tex`, §2–4). I opened nothing in `src/` or `tests/`. All code is in `verify/va/`. Everything ran in **1 process** (Java with `-XX:ActiveProcessorCount=1`), run sequentially on 2026-11-24/25. Total wall-clock about 4.2 h.

Files:
- `verify/va/EnumVA.java`: enumerator, local tests, canonical form.
- `verify/va/va_alg.py`: affine tests, W0, division, exclusion, finite criterion, certificate replay.
- `verify/va/selftest.py`: polynomial identities at random rational point sets.
- `verify/va/run_all.py`, `verify/va/run_e2u.py`: drivers.
- `verify/va/out/`: `va_patterns6_local.txt` (sha256 97947b7e…649640c), `va_patterns6_affine.txt` (sha256 db201dcd…c179c5a1), `classify_retained.jsonl`, `replay_all.jsonl`, `classify_excluded_sample.jsonl`, `classify_e2unresolved.jsonl`, `run_all.log`.

### Claim 1 — enumeration: CONFIRMED
- **Method.** An extension generator in Java (k = 3..6) with full-pattern local tests applied to partial labellings. Every test is a necessary condition, so partial rejection is sound:
  - L(i): at each point, at most 3 pairs share a label.
  - L(ii): for each pair a,b and each ordered label pair (u,v), at most 2 points x.
  - L(iii): for each pair a,b, at most 2 points x with c_ax = c_bx.
  - Lemma 2.2: no 4-set has all three perfect matchings equal-labelled.
  - At most 6 labels.
- **Canonical form.** Brute-force minimum over all k! point orderings of the first-occurrence-renamed lex string. This is the same convention as the data files, so no mapping was needed.
- **Affine tests.** My own Python/flint implementation of Lemma 2.3:
  - R holds the parallelogram rows from every ordered quadruple, plus the centroid rows.
  - T1: rank R_{I^c} < rank R for some triple I.
  - T2: e_l or e_l − e_m lies in the row space of A, where A comes from Σρ_i(D_ij − D_i1) = 0 over a basis of the row space of R.
- **Result.** k = 3/4/5/6 gives 3 / 22 / 932 / **625,472** patterns. By r at k = 6: 275 / 18,938 / 166,912 / 439,347.
  - Affine survivors: **344,835**, by r 4 / 2,191 / 60,895 / 281,745.
  - Rejections: 122,907 support, 92,231 zero, 65,499 equal. These are identical to builder E's breakdown.
- **Set comparison** with `data/patterns6_local.txt.gz` and `data/patterns6_affine.txt.gz`: **0 differences in either direction** for both files.

### Claim 2 — exclusion certificates: CONFIRMED (all 340,198)
- **Replay checker.** I wrote my own code for each piece, exact over Q (python-flint `fmpz_mat` rank, rref and nullspace):
  - D, and L_ij = D_1i + D_1j − D_ij.
  - All 100 3×3 minors, expanded in the 15 pair variables and then substituted. They give 55 distinct minors, as the paper states.
  - W0 = span{z_l·minor} + span{(row-space basis of A)·cubic monomial}.
  - Division by f: C_f = left kernel of [f·m_i ; W], restricted to the q-part. W is then enlarged by m·q for every monomial m with deg m = deg f.
  - Factors: z_l, z_l − z_t, and tri_ijk = 4D_ijD_ik − (D_ij + D_ik − D_jk)².
  - Targets: a monomial e_m, or det4 = det(D on 4 points), tested by rank([W; v]) = rank W.
  - Name conventions inferred from the data: z indices and tri indices are 1-based, det4 indices are 0-based.
- **Self-test.** At 20 random rational 6-point sets:
  - every minor vanishes;
  - det4 = −4(det B)² for B with rows (1, x, y, x² + y²);
  - tri_ijk = 4·det(p_j − p_i, p_k − p_i)².
- **Coverage.** I replayed **every one of the 340,198 certificates**, not a sample. 330,975 have an empty trace (the target already lies in W0). The rest have traces of 1–9 steps.
- **Result.** **340,198 / 340,198 verify.** In every step the dimension of my replayed W matched the certificate's recorded (before, after) dimensions exactly. 0 failures. No det4 target is identically zero.
- **Coordinator's priority set.** `data/E2_unresolved_excluded.txt`, the 4,923 patterns where Engine 2 was UNRESOLVED: **4,923 / 4,923 certificates verify** with dimensions matching. In addition, my certificate-free classifier (Claim 3) independently EXCLUDES **4,923 / 4,923** of them.
- **Partition.** The certificate patterns, `E1_F` and `E1_U` are pairwise disjoint, have no duplicates, and their union equals my affine set exactly.

**Soundness proof (written).** Let p_1..p_6 be a general-position realization, and let z* in R^r be its label values. These are positive and pairwise distinct (definition of realization). Let I* = {g in Q[z]_4 : g(z*) = 0}. I* is a Q-subspace.
1. **W0 ⊆ I*.**
   - Translate so that p_1 = 0. Then L(z*)_ij = |p_i|² + |p_j|² − |p_i − p_j|² = 2 p_i·p_j, so L(z*) = 2PPᵀ, where P (5×2) has the vectors p_i as rows. So rank L(z*) ≤ 2, every 3×3 minor vanishes at z*, and so does z_l·minor.
   - Every row ρ of R is a true affine relation at a general-position realization (Lemma 2.2 and the centroid argument; rows sum to 0). So Σρ_i|p_i − p_j|² = Σρ_i|p_i|² − 2p_j·Σρ_i p_i + |p_j|²Σρ_i is independent of j, and every row of A vanishes at z*. Hence (A-row)·m is in I*.
2. **Division keeps W ⊆ I*.**
   - If W ⊆ I*, q in C_f and f(z*) ≠ 0, then f(z*)q(z*) = 0 forces q(z*) = 0, so m·q is in I*.
   - f(z*) ≠ 0 in each case:
     - z_l > 0;
     - z_l ≠ z_t, because the labels are distinct;
     - tri_ijk(z*) = 4|u|²|v|² − 4(u·v)² = 4 det(u,v)², with u = p_j − p_i and v = p_k − p_i. This is nonzero because no 3 points are collinear. It is symmetric in the choice of apex (Heron), so the labelling convention is immaterial.
   - By induction over the trace, the final W ⊆ I*. The certificate needs no claim that the trace is a full fixpoint: any sequence of sound steps works.
3. **Targets are nonzero at z*, a contradiction.**
   - A monomial z^α is a product of positive numbers, so it is nonzero.
   - For points a,b,c,d, D_abcd(z*) = B J Bᵀ with det J = −4 (checked symbolically in the self-test), so det = −4(det B)².
   - If det B = 0, some (λ0, λ1, λ2, λ3) ≠ 0 has λ0 + λ1x + λ2y + λ3(x² + y²) = 0 at all 4 points.
     - λ3 = 0 would put 3 (indeed 4) points on a line (or force λ = 0), which is excluded.
     - λ3 ≠ 0 gives a circle (radius² > 0, because distinct points lie on it) through 4 points, which is excluded.
   - So det4(z*) ≠ 0.
   - Hence a target in W ⊆ I* is impossible, and the pattern has no general-position realization.
4. **Verifying the certificate.** Verifying it is exactly re-deriving "target in W" over Q, which is what the replay does. Rank over Q of integer matrices is exact in flint.

### Claim 3 — retained patterns: CONFIRMED
- **Method.** A full certificate-free pipeline: affine → W0 → complete passes over all factors until a pass leaves dim W unchanged (the least fixpoint, which is order-independent) → exclusion test.
  - The monomial test is a rref row with one nonzero entry. This is equivalent to "contains a monomial".
  - The det4 test runs over the 15 quadruples.
  - Finite criterion: rank of W restricted to the monomials not divisible by z_a equals C(r+2, 4). F if r in {5,6} and some a passes. The smallest such a is recorded.
- **Results.**
  - On all 4,637 retained patterns: **0 excluded**. Builder F → F for all **4,450**. Builder U → U for all **187**, by r **3 / 9 / 175** (r = 4/5/6).
  - For every F pattern, my smallest a and my h = dim V/W equal the builder's recorded a and h (0 mismatches).
  - On 6,000 uniformly random excluded patterns (seed 217), plus all 7,577 excluded patterns whose certificate needs ≥ 2 division steps (the stress set; 13,442 distinct in total): **13,442 / 13,442 EXCLUDED**.
  - On the 4,923 E2-unresolved patterns: 4,923 / 4,923 EXCLUDED (see Claim 2).
- **Remaining limit.** My classifier was not run on all 340,198 excluded patterns. That would have taken about 5 CPU-h at 1 worker, and the certificate replay above already covers all of them.

**Overall V-a verdict.**
- Phase 2 counts (625,472 / 344,835): **CONFIRMED**, identical sets.
- All 340,198 quartic exclusions: **CONFIRMED** by exact replay, with a written soundness proof.
- F/U split 4,450 / 187 (3/9/175): **CONFIRMED** by an independent classifier.
- No discrepancy of any kind was found.
