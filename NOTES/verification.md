# Verification log (verifier / red team)

## Phase 1 — checker

**Independence:** written from PROMPT.md §1 only. I did not open anything in `src/` or `tests/`.

### What I implemented: `verify/checker_v.py`
- One generic algorithm with two exact arithmetic modes:
  - `Q3Field`: an element a+b√3 is a pair of `Fraction`s (a, b). +, −, ×, ÷ are hand-written. Division multiplies by the conjugate over the norm a²−3b². Zero test: a=b=0, which is valid because √3 is irrational, so the norm is 0 only at 0. There is also an exact sign test: when a and b have opposite signs, compare a² with 3b².
  - `QField`: plain `Fraction`s (pure-rational mode).
- Checks, in order:
  1. Duplicate points are rejected. The spec does not cover them; I treat a zero distance as invalid.
  2. No 3 collinear: an exact cross product.
  3. No 4 concyclic, tested two independent ways, neither of them the 4×4 lifted determinant:
     - (i) the circumcentre O of A,B,C by Cramer's rule with exact division, then test |OD|² = |OA|²;
     - (ii) cross-ratio realness, Im[(A−C)(B−D)·conj((A−D)(B−C))] = 0, which has no division.

     Once no 3 points are collinear, the two tests must agree, and the code raises an error if they don't. They never disagreed in any run.
  4. Squared distances are grouped by exact equality. The checker requires exactly n−1 classes and sorted class sizes equal to [1..n−1].
- `verify/run_checks.py` runs the test cases. `verify/lattice_search.py` runs the lattice search.

### Results: `run_checks.py`, 18/18 PASS
- **ACCEPT:** Palásti n=8 exactly as printed in 1509.07220 Fig. 1. Class sizes are [1..7], so **the printed coordinates are a valid crescent configuration**. The checker also accepts the configuration after a 30° rotation plus an irrational shift, after a 60° rotation, and after a mirror image. Those transformed copies exercise the full Q(√3) arithmetic. It also accepts isosceles triangles in both modes, and a rational n=4 example, (0,0),(0,1),(1,1),(1,2).
- **REJECT:**
  - Palásti with a duplicate point.
  - Palásti with its last point replaced by (3√3,0), which is collinear with points 1 and 2 (reason: collinear).
  - Palásti with its last point moved to (√3,15/7) (13 distances).
  - Palásti's first 7 points (7 distances, not 6).
  - Palásti with point 2 replaced by the antipode on the circumcircle of points 0,1,3 (reason: concyclic (0,1,2,3)).
  - The equilateral triangle.
  - The unit square (concyclic).
  - A rational 3-collinear set.
  - A rational duplicate.
  - (0,0),(0,2),(1,1),(3,1): in general position with exactly 3 distances, but class sizes [2,2,2] (reason: class sizes).
- **Mode agreement:** Q mode and Q(√3) mode gave identical verdicts and reasons on all 12,650 four-point subsets of the 5×5 integer grid.

### Lattice search: radius-3 hexagonal patch, 37 points (`lattice_search.py`)
Lattice point (i,j) ↦ (i + j/2, (√3/2)·j). The search is an exhaustive backtrack with sound pruning: no 3 collinear, at most n−1 distinct distances, and the k-th largest class size at most n−k. Every leaf is re-checked with the exact Q(√3) checker. Solutions are reduced by the 12 lattice symmetries and translation.

| n | raw solutions | classes mod symmetry | time |
|---|---|---|---|
| 5 (radius 2, 19 pts) | 132 | 5 | 0.2 s |
| 6 | 1344 | 26 | 5.6 s |
| **7** | 432 | **8** | 5.9 s |
| 8 | 36 | **1**, which is Palásti's configuration: its canonical form equals that of the printed coordinates after rotation | 4.5 s |
| 9 | 0 | 0 | 6.7 s |

The n=9 result on this 37-point patch is consistent with the published 91-point result (1509.07220, Remark 3.1). The full n=7 and n=8 lists are in `verify/lattice_n7_r3.txt` and `verify/lattice_n8_r3.txt`.

Exact Cartesian coordinates for two of the n=7 classes:
- Class 0: (0,0), (1/2,√3/2), (−1/2,−3√3/2), (0,−2√3), (1/2,−3√3/2), (7/2,−√3/2), (4,−√3)
- Class 2: (0,0), (1/2,√3/2), (0,−√3), (3/2,−√3/2), (5/2,√3/2), (3/2,−3√3/2), (5/2,−√3/2)

### Definition notes and ambiguities
- **"No 4 concyclic" when 3 points are collinear.** Under the definition this case never arises, because collinear triples are rejected first. If someone instead used the generalized "circle or line" test (cross ratio real), 4 collinear points would also count as "concyclic", and either reading rejects them. A builder whose code tests concyclicity *before* collinearity, with a lifted determinant, will still reject. It may report a different reason, but the verdict is the same.
- **"Some distance occurs exactly i times for each i"** together with "exactly n−1 distances" is equivalent to "the class sizes are exactly {1..n−1}". The counts are forced, since they sum to C(n,2).
- **Duplicates** are not addressed in the spec. I reject them.
- **Coordinate fields.** Triangular-lattice configurations fit in Q(√3), and all squared distances are rational there. A candidate for n ≥ 9 found by algebraic methods may need a larger number field. This checker cannot handle that case: it would need exact algebraic-number arithmetic, for example minimal polynomials with isolating intervals.

Cross-check vs builder checker: DONE 2026-11-24 (orchestrator ran `verify/crosscheck_checkers.py`): all 57 saved examples (builder data/known_examples.json n=3..8 incl. Palásti; verifier's 8 n=7 classes) accepted by BOTH checkers; ~1.2k one-point-perturbed variants: 0 disagreements (1034 joint rejects). Both exhaustive radius-3 lattice searches agree: k=6:26, 7:8, 8:1, 9:0. **Phase 1 (rule 3) SIGNED OFF.**


## Phase 3/4 — V-b (red team + independent re-implementation), 2026-11-25

**Independence:** I did not open anything in `src/` or `tests/`. Inputs used: NOTES/*.md, Wang's TeX (spec only), the data files `data/E1_U.txt.gz`, `data/E1_F.txt.gz`, `data/F_recon_final.jsonl.gz` (only the fields q, z, G, points, status, which I re-verified exactly), `data/F_extend.jsonl.gz` (only the near-miss coordinates), and `data/coverage_E1_found.txt` / `coverage_E1_level7.txt` (only for comparison). All code is in `verify/vb/`. Every run used 1 worker process.

### 1. Red team of the proofs (NOTES/coverage.md, NOTES/extension.md)

| # | Concern | Finding | Verdict |
|---|---|---|---|
| R1 | Arrangement list (7 cases) | Re-derived myself: classified all 36·C(35,2) = 21,420 labelled (δ1, {δ2,δ2}) structures by a complete isomorphism invariant (degree multiset of the 3-edge graph, degrees of the δ1 endpoints, number of endpoints). Exactly 7 classes, sizes 252, 1512, 1512, 3024, 3780, 3780, 7560, same as builder's orbit sizes. My own 7 representatives (low pairs on points 1–3 and 7–9, different from builder's and Wang's) hit each class once. | CONFIRMED |
| R2 | "Every all-high 6-subset is a core, so it is in F ⊔ U" (Step 2) | Correct: a core is any 6-subset with D(S) ⊆ {δ3..δ8} (Wang Def 3.1), it inherits GP, has ≤ 6 distances. Case split "some core in F" / "all cores in U" is exhaustive. F ∩ U overlap would be harmless. | CONFIRMED (conditional on (A)) |
| R3 | Seeds / first-occurrence numbering (Steps 3–5), lookahead N4 | Logic checked. U-membership uses every ordering of every U member, so canonical-form choice cannot lose anything. The lookahead tries a superset of the true last label. Empirically: my search (no lookahead, different branching) gives the identical 487-pattern set. | CONFIRMED |
| R4 | Pruning rules N1a–d, N2 capacity, N5 | Each re-proved (Lemma 2.1 (i)–(iii), Lemma 2.2 via SSS isometry, Lemma 5.2 Hall-type). N2 formula m_i ≤ 8−s+i is Lemma 5.2 with b = 3..8. Low labels: counts fixed by construction. | CONFIRMED |
| R5 | N6 affine post-filter only at 9 points | Sound: rows of R/A for a prefix are rows of the full R/A (indices of prefix points are a subset), so a prefix rejection implies the full rejection; omitting tests can only add survivors. T1 (support-≤3 affine relation ⇒ coincidence/collinearity) and T2 (Σρ_i D_ij independent of j) re-derived. | CONFIRMED |
| R6 | Symmetry breaking | Builder uses none beyond the arrangement (redundant but complete). I used a different one (orbit dedup under the arrangement stabiliser at every level) and got the same survivors. | CONFIRMED |
| R7 | Lemma H (≥ 2 high-linked outside points) | Low pairs from distinct outside points to S are distinct pairs; two non-high-linked points need ≥ 4 low pairs > 3. | CONFIRMED |
| R8 | Lemma H(b), m = 0, only triangles (1,2,3),(4,5,6) | A high-linked x has ≤ 1 low partner, which lies in at most one of the two disjoint triangles; the other triangle carries D(S) = H values. Lemma T gives x ∈ K². | CONFIRMED |
| R9 | |D(S)| = 5 case (Lemma H(c)) | ≥ 5 high partners, values in D(S) ∪ {d0}: ≥ 3 in D(S) (trilateration, all 20 triangles) or ≥ 3 at d0 (circumcentre). Correct. Note it is not tested end-to-end by a known configuration (builder's own caveat); my independent search uses a different m = 1 argument (see 3) and agrees. | CONFIRMED (proof); test coverage PARTIAL |
| R10 | Lemma L (last point) with covering family 𝒯(8,5) of 9 triangles | Proof correct: either an absent value with multiplicity ≥ 3 (all its pairs contain w ⇒ circumcentre) or ≥ 5 partners with D(Q) values ⇒ some triangle of 𝒯 inside them. Requires Q GP (pruned earlier). | CONFIRMED |
| R11 | Capacity Lemma C (m = 0 / m = 1 split with k, d0) | Necessary conditions; my extension does not use it (uses a generic bipartite matching sizes→{1..8}, core values ≥ 3) and gets the same final answer. | CONFIRMED |
| R12 | "One abstract K_q search covers all PD embeddings" | Correct (all decisions are zero tests in K_q; embeddings are injective). My run instead searches each PD real embedding separately and agrees. | CONFIRMED |
| R13 | Upstream assumption (A) = Wang Cor 4.5 for the computed F/U | Not provable by these notes; it rests on the Phase-2/E1 computations. Independent spot check below (item 4). Also: 6 F patterns reconstructed only by the E1 route; the reconstruction step (factor list complete) is trusted here (Hermite-rank certificate by builder), not re-derived by me. | PARTIAL (spot check passes; no proof) |
| R14 | Real-root / PD decisions for non-realized factors | I recomputed real roots (arb) and PD at each real root for all 4,492 factors: see item 3. | see 3 |

**Red-team conclusion:** I found no logical gap that lets a 9-point crescent configuration escape both cases. The whole argument still depends on (A), i.e. that the E1 split F ⊔ U contains every realizable GP 6-point pattern with ≤ 6 distances (R13). The spot check in item 4 found no counterexample to (A).

### 2. Independent coverage search (Case 2): CONFIRMED

Code: `verify/vb/coverage_vb.py`, `common_vb.py`, `compare_cov_vb.py`. Logs: `verify/vb/cov_E1.log`, `cov_E1_affeach.log`, `compare_cov_vb.log`. Survivors: `verify/vb/cov_E1_found.txt`, `cov_E1_level7.txt`.

How my design differs from builder §(b):
- I use my own arrangement representatives. δ1/δ2 pairs (1-based): disjoint 17 / 28, 39. Triangle 89 / 78, 79. Star 19 / 29, 39. Chain, δ1 middle 89 / 18, 29. Chain, δ1 end 18 / 89, 29. P3+K2, δ1 disjoint 38 / 19, 29. P3+K2, δ1 meeting 19 / 29, 38.
- The pattern grows one point at a time. For each new point, a base of 5 high partners is joined against a five-point lookup table built from all orderings of U. This is Wang's device, and the builder did not use it. The remaining partners are then assigned one by one, with a U-membership check on every completed all-high 6-subset.
- After every point, states are reduced modulo the stabiliser of the arrangement restricted to the points present. The builder keeps all labelled states.
- The full local tests (Lemma 2.1 (i)–(iii), Lemma 2.2) are applied at every point, then capacity (Lemma 5.2, capacities 3..8), and at 9 points the class sizes must be exactly 3..8.
- The affine tests are my own implementation of Lemma 2.3 with parallelogram rows only: T1 is a row-space vector of support ≤ 3, T2 is z_l = 0 or z_l = z_m forced by Az = 0. Ranks are exact integer ranks (flint fmpz_mat).
- Canonical forms are my own: canon6 is the minimum over 720 orderings of the first-occurrence string; complete 9-point patterns are canonicalised by min over permutations within cells of a point invariant.

Results (U = E1_U, 187 lines, all canonically distinct). The level counts are orbit-reduced, so they cannot be compared with the builder's labelled counts.

| arrangement | states at level 6 → 7 → 8 → 9 | distinct 9-pt patterns before affine | builder | after affine (end) | Wang-style: affine after every point, level 8 → 9 |
|---|---|---|---|---|---|
| disjoint | 7389 → 28295 → 36077 → 222 | **222** | 222 | **0** | 2125 → 0 |
| triangle | 187 → 69 → 142 → 3 | **3** | 3 | **0** | 39 → 0 |
| star | 7389 → 930 → 0 → 0 | **0** | 0 | **0** | 0 → 0 |
| chain, δ1 middle | 2018 → 196 → 451 → 12 | **12** | 12 | **0** | 119 → 0 |
| chain, δ1 end | 3711 → 376 → 838 → 20 | **20** | 20 | **0** | 226 → 0 |
| P3+K2, δ1 disjoint | 7389 → 930 → 1065 → 68 | **68** | 68 | **0** | 282 → 0 |
| P3+K2, δ1 meeting | 14626 → 1830 → 2094 → 162 | **162** | 162 | **0** | 557 → 0 |
| **total** | | **487** | 487 | **0** | 0 |

- **The set of 487 patterns is identical to `data/coverage_E1_found.txt`** after canonicalising both with my own canonical form, arrangement by arrangement.
- My affine code kills all 487: 384 by T1 and 103 by T2, the same split as the builder's.
- **7-point stage:** 13 canonical 7-point all-high patterns, all of whose 6-subsets are in U and which pass local + capacity. **The set is identical to `data/coverage_E1_level7.txt`.** All 13 also pass the 7-point affine tests, so this is consistent with Wang's 13, which he computes with the affine tests applied.
- Runtime: about 4 min (end-affine run) plus about 5 min (per-level affine run), on 1 worker.

**Verdict: CONFIRMED.** Conditional on (A), no 9-point crescent configuration has all its cores in U.

### 3. Independent extension check (Case 1): CONFIRMED (full run, not a sample)

Code: `verify/vb/extend_vb.py`, `summarize_ext_vb.py`. Output: `verify/vb/ext_full.jsonl` (one record per factor), `ext_full.stdout`, `ext_summary.txt`. Runtime 967 s on 1 worker, covering all 4,450 patterns, 4,492 factors and 7,829 embeddings. No reduction to a sample was needed.

How the design differs from extension.md:
- **No Lemma H, no covering family 𝒯, no Lemma C.** I use my own case split instead.
  - **m = 0 (|D(S)| = 6):** every outside point has ≥ 3 high pairs to S, and H = D(S). So every outside point is a trilateration point of some triangle of S, using any of the 20 triangles, with values in D(S)³. All three outside points are drawn from this one set, under the exact budget "exactly 3 of the 36 pairs carry a value outside D(S)", since these are exactly the low pairs.
  - **m = 1 (|D(S)| = 5):** at most one outside point has ≥ 2 low pairs to S, so at least two outside points y, z have ≥ 3 partners with values in D(S) (trilateration) or ≥ 3 partners at d0 (circumcentre). The third point w either has ≥ 5 partners in Q = S ∪ {y,z} with values in D(Q), at least 3 of them in S, which gives trilateration on the 20 core triangles with D(Q)³; or every pair carrying d0 (≥ 3 of them) contains w, which makes w a circumcentre of one of the 56 triangles of Q. Pairs {y,z} are filtered by a generic injective size assignment D(Q) → {1..8} (bipartite matching, core values ≥ 3) and by "≤ 3 values outside D(S)".
  - No GP pruning is done before the 9-point stage, so my funnel is wider than the builder's.
- **I search each PD real embedding separately.** The builder runs one abstract K_q search instead. I recompute the real roots of q myself (arb, conjugate-isolation test) and test PD of G at each root.
- **Screening runs in float64 Euclidean coordinates** (Cholesky of G). Every hit is rebuilt and re-checked **exactly in K_q** by my own arithmetic (fmpq_poly mod q): 36 distances, 84 collinearity determinants and 126 four-point determinants with rows (1, x, y, pᵀGp).
- Numerical decisions use tolerance 1e-8 relative to the largest core value. Over the whole run, every value accepted as equal was within 2.4e-12 and every value rejected was off by ≥ 1.85e-6, so the gap is about 6 orders of magnitude and no decision was borderline. Every numeric hit was confirmed exactly: 879 of 879, with 0 rebuild failures.

| quantity | V-b | builder |
|---|---|---|
| factors / REALIZED / no real root / G never PD | 4,492 / 3,879 / 374 / 239 | 4,492 / 3,879 / 374 / 239 |
| real-root count per factor | equal for all 4,492 factors | |
| PD real embeddings | **7,829** (per-factor equal, 0 undecided) | 7,829 |
| realized cores re-verified exactly (15 distances = z-values, GP, r labels) | 3,879 / 3,879 | — |
| m = 0 / m = 1 realized factors | 3,720 / 159 | r=6: 3,720 / r=5: 159 patterns |
| trilateration trials (core) / distinct candidates | 33.3 M / 429,734 | (different method) |
| candidate pairs passing budget/capacity | 273,967 | — |
| last-point trials (m=1) / candidates | 215 M / 2.92 M | — |
| 9-point sets with exact class sizes 1..8 (no GP pruning; counted per core × embedding) | 879 (195 of them with no 3 collinear; 23 distinct δ-signatures) | 5 (with GP pruning at 7 and 8 points) |
| **9-point crescent configurations** | **0** | **0** |

The comparison of the funnels is limited. Only the final count (0) and the realized-core and embedding counts are directly comparable, and those agree exactly. My 879 is larger than the builder's 5 because I do not prune GP at 7 and 8 points.

Both builder near-miss configurations appear among my exact hits. Config 2 (δ-signature 0.065…0.935, 1) is reached from 14 cores, including 111232433454125, 111232442554362 and 112232145661354. Config 1 is reached from core 111231145613562. Two of the builder's near-miss keys produce no hit in my search: config 1 via 111231134516425 and config 2 via 112233241524646. I checked these exactly: in both, the 6-set S carries a δ2 pair (the class sizes met by S are {2,4,5,6,7,8}), so S is **not** a core of that 9-set. My core-strict search correctly does not produce them, and the builder's search is simply more permissive. This is not a discrepancy in completeness.

**Near-miss sets, exact re-check (`verify/vb/nearmiss_vb.py`, `nearmiss_vb.json`, my own K_q code). CONFIRMED.**
- All 5 builder records have class sizes exactly [1..8] and 0 collinear triples.
- Concyclic quadruples:
  - Config 1, q = u² − 6u + 25/3: 3 concyclic quadruples, (0,1,6,8), (0,4,7,8) and (2,3,7,8). Both real roots are PD, and δ1..δ8/δ8 ≈ 4.633, 1.367, 0.5505, 3, 1.8165, 0.1835, 2, 1 at one root.
  - Config 2 (4 records, fields q = u² − 17u/5 + 1, u² − 8u − 5, u² + u − 5, u² − 4u − 16/3): 2 concyclic quadruples each, with δ/δ8 ≈ 1.7913, 1.6748, 0.7913, 0.2087, 0.4417, 0.1165, 0.5583, 1, matching NOTES/extension.md.
- `checker_v.py` cannot be used here: it supports only Q and Q(√3), and these fields are not Q(√3). The check was therefore done with my independent K_q implementation.

**Verdict: CONFIRMED, and now rigorous; see 3b.** No realizable F-core extends to a 9-point crescent, by a different reduction and a different search. The float64 run described above was the first pass. Following the coordinator's follow-up (rule 1), it has been superseded by the certified re-run in 3b, which reproduces every count exactly.

#### 3b. Rigorous re-run, ball arithmetic plus exact number field (2026-11-25)

Code: `verify/vb/extend_rig_vb.py`. Output: `verify/vb/ext_rig.jsonl` (3,879 records, one per realized factor) and `ext_rig.stdout`. Same case split and search order as item 3, but **no float64 is used in any decision**:
- **Ball arithmetic.** All quantities are python-flint `arb` balls computed from the exact K_q data at the embedding's root. The root is isolated by `complex_roots` at the working precision. The run starts at 128 bits and escalates per item to 256 and then 512 bits.
- **Rejecting a candidate** requires the ball to exclude 0: the trilateration residual, or a "different value" decision where two balls are disjoint. A 9-point set is rejected when its 36 distance balls form more than 8 overlap components; balls in different components are disjoint, so this certifies more than 8 distinct values.
- **Anything the balls do not decide** goes to the exact K_q test: `ex_trilat`, exact distance keys (reduced polynomials mod q), exact collinearity, and the full `ex_check_config`.
- **Explained triples.** A triple is treated as an exact solution without a new test when an already-constructed exact point P satisfies its three distance equations, checked by exact key equality. The linear system has a unique solution, so P is that solution.
- **PD of G** is decided by balls. An undecided case would abort the run; none occurred.

| decision | by balls | by exact K_q test |
|---|---|---|
| trilateration residual (248,632,740 trials) | 236,786,999 rejected at 128 bits (0 needed 256 or 512 bits) | 3,242,211 explicit exact tests + 8,603,530 exactly explained by a known exact point |
| candidate-pair distance vs D(S) (m=0) | 313,571 | 765,073 (incl. all m=1 pair distances) |
| circumcentre degeneracy | — | 18,411 collinear triples decided exactly |
| final 9-point test (3,029,989 sets) | 2,937,013 rejected (>8 certified distinct values) | 91,976 full exact checks |
| PD of G at 7,829 real roots | 7,829 (0 undecided) | — |

**Results.**
- The funnel is identical to the float run: 429,734 distinct candidates, 273,967 pairs, 104,801 budget triples, 2,924,188 last-point candidates, and 879 exact 9-point sets with class sizes 1..8.
- **0 crescent configurations.**
- Runtime 1,979 s (33 min) on 1 worker.

The extension check is therefore an **independent exact certificate**: every rejection is either a ball that provably excludes 0 or an exact computation in K_q.

### 4. Spot check of (A) = Wang Cor 4.5 on the computed lists: PARTIAL (no counterexample)

Code: `verify/vb/cor45_vb.py`, `cor45b_vb.py`. Output: `cor45_vb.json`, `cor45b_vb.json`.

Exact GP 6-point sets with ≤ 6 distances were taken from four sources:
- all 6-subsets of my 8 radius-3 lattice n=7 crescents and of Palásti n=8;
- 3,512 two-concentric-C3-orbit sets in Q(√3);
- 300k random 6-subsets of the radius-3 triangular patch (2,081 qualified);
- **all 8,523 GP ≤6-distance 6-subsets of the 879 exact 9-point sets from item 3.** These lie in non-lattice number fields of degree up to 12, a regime the Phase-2 soundness test (lattice forms only) did not exercise.

Results:
- Lattice and C3 sources: 127 distinct patterns (56 in F, 71 in U), 0 missing from F ⊔ U.
- 9-point-set subsets: 252 distinct patterns (186 in F, 66 in U), 0 missing.

**Verdict: PARTIAL.** This is evidence, not a proof. (A) still rests on the upstream E1 computation, and the reconstruction completeness for F, the factor lists, is trusted (builder's Hermite-rank certificate), not re-derived by me.

### V-b overall
- Coverage search: **CONFIRMED** (487 / 13 / 0 reproduced exactly; identical sets).
- Extension: **CONFIRMED** (0 crescents; 3,879 / 7,829 reproduced; near-misses verified exactly). Since 3b this is rigorous: arb balls plus exact K_q, with no float decisions.
- Red team: **no gap found**.
- (A): **PARTIAL**, spot-checked with no counterexample, and it remains the only unproven link.

## Phase 2/3 — V-a (independent enumeration, certificate replay, classifier), 2026-11-25
Collected by the orchestrator from V-a's outputs in `verify/va/out/`. V-a's code is in `verify/va/` and never read `src/` or `tests/`. V-a's own narrative write-up is pending; the numbers below come straight from its output files.
- **Enumeration** (own Java enumerator `EnumVA.java`, own canonical form and affine tests):
  - counts 3 / 22 / 932 / 625,472 / 344,835;
  - the set comparison with the builder's lists shows 0 differences: `va_patterns6_local.txt` equals `data/patterns6_local.txt.gz`, and the affine stage gives mine−data 0, data−mine 0;
  - affine rejections by type: support 122,907, zero column 92,231, equal columns 65,499;
  - **CONFIRMED.**
- **Certificate replay** (own exact replay over Q, `replay_all.jsonl`): **all 340,198 certificates verified (340,198 True, 0 False)**, including all 4,923 patterns on which E2 was initially UNRESOLVED. Runtime 6,914 s on 1 worker. **CONFIRMED.**
- **Independent quartic-span classifier on the 4,637 retained patterns** (`classify_retained.jsonl`): agrees with the builder on every pattern, F 4,450 and U 187 (r=4: 3, r=5: 9, r=6: 175). It excludes none of them. **CONFIRMED.**

## Phase 3 — Engine 2 (different elimination engine) over all 344,835 patterns, 2026-11-25
- **First pass** (limits 1e4/1e5/1e5, 120 s timeout): of the E1-excluded patterns, 335,275 EXCLUDED and 4,923 UNRESOLVED. **0 came out ZERODIM/POSDIM**, so E2 never contradicts an E1 exclusion.
- **Retry of the 4,923** (limits 1e5/1e6/1e6, 600 s timeout; `data/E2_retry/out.jsonl`, sha256 404e6f63…db779): **4,923 / 4,923 EXCLUDED**.
- ⇒ **Every one of the 340,198 E1 exclusions is independently certified by a Gröbner/saturation computation** (GB = {1}, or a product of the nonzero factors lies in the ideal).
- On the retained patterns: F → ZERODIM 4,422 + UNRESOLVED 28; U → POSDIM 178, ZERODIM 2, UNRESOLVED 7. Unresolved retained patterns are harmless, because retaining is always the safe direction.

## Factor-count discrepancy 4,118 vs 4,119: RESOLVED
- The E1 multiplication-matrix route (builder D's cross-check `data/F_recon_e1check.jsonl`, all 4,450 patterns, status ok and match True) gives 4,493 factors and **4,119 with a real root, exactly Wang's number**. The E2 route gives 4,492 and 4,118.
- The only difference is pattern **112233145565546**. E1 gives 2 factors, both with real roots, 1 realized. E2 gives only the realized factor, q = u² − 5u + 37/16. The extra E1 factor is a component removed by E2's stronger saturation; its non-degeneracy factors include 4-point det(D) and the triangle areas.
- Both routes give **3,879 realized cores**, so the extension input is identical and the n=9 result is unaffected.

## Overall sign-off for goal A
Every headline computation has two independent implementations that agree:
- checker (Phase 1);
- enumeration, with 3 generators;
- exclusions: E1 certificates replayed by V-a, plus E2 Gröbner certification of all 340,198;
- F/U classification: E1 and V-a classifier;
- coverage search: builder C and V-b, the same 487 → 0 and the same 13 seven-point patterns;
- reconstruction and extension: builder D, and V-b's rigorous arb plus exact implementation, giving the same 4,492 / 3,879 / 7,829 → 0.

The red team (V-b) found no gap in the written arguments.

**Remaining caveat (Cor 4.5):** "the computed F ⊔ U contains every realizable core" is not a separate theorem. It follows from (i) every exclusion being a sound necessary-condition certificate (proved on paper: Wang Lemmas 2.1–2.3, 4.1–4.2, re-derived by the reader, V-a and V-b), plus (ii) the enumeration being complete. Both have now been checked by two independent implementations. The empirical spot checks are consistent with this: 145 + 127 + 252 realizable canonical patterns, including non-lattice ones over fields of degree ≤12, all lie in F ⊔ U.

**SIGNED OFF (orchestrator, collecting V-a/V-b outputs): A1 conditions met.**

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

**Sign-off addendum (orchestrator, after V-a's own report):** V-a's own write-up (the "Phase 2/3 — V-a" section added by V-a) confirms the numbers cited in the sign-off above:
- enumeration set-identical;
- all 340,198 certificates verified, including the 4,923 E2-initially-unresolved patterns, which V-a's certificate-free classifier also excludes;
- the classifier agrees on all 4,637 retained patterns, with the same a and h for every F pattern;
- the classifier also excludes all 13,442 excluded patterns it was run on (6,000 random plus all 7,577 with ≥2 division steps).

Sign-off stands. No job is running.
