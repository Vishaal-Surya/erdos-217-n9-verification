# Coverage search (Builder C): independent re-implementation of Wang Prop 5.1 / Computation 5.4

Date: 2026-11-24. Spec: arXiv:2611.03688 §5 (read as a spec only), NOTES/wang_method.md §6, NOTES/phase2_counts.md.
Code: `src/coverage/` (`CoverageSearch.java`, `arrangements.py`, `targets.py`, `planted.py`). Tests: `tests/test_coverage.py`.
Data/logs: `data/coverage_*`.

**Status (2026-11-24, real run done, §5):** with U = `data/E1_U.txt.gz` (187 patterns, sha256 cbd6781f…4214a), the search keeps 45,216 labelled leaves, which are 487 distinct complete 9-point patterns. **All 487 fail the exact parallelogram affine tests (N6), so 0 survive in all 7 arrangements.** This agrees with Wang's Computation 5.4.

---

## Part (a): the coverage proof

### Notation

- P = {p_1, …, p_9} ⊂ ℝ² is a 9-point crescent configuration: general position, and the 36 pairs split into 8 distance classes of sizes 1, …, 8. δ_k denotes the squared distance of the class of size k.
- The **pattern** of P is the map c: pairs → classes. Two patterns are **equivalent** if they differ by a permutation of the points and a bijective renaming of the labels.
- **Low pairs**: the 3 pairs with distance δ_1 or δ_2. All other pairs are **high**.
- A **core** (Wang Def 3.1) is a 6-subset S ⊆ P whose 15 pairs are all high.
- F and U are the two families of 6-point patterns produced upstream (E1 builder / Wang Procedure 4.4).

**Upstream assumption (A).** This is Wang's Corollary 4.5, and the U file must be built so that it holds. Every 6-point pattern realised by a 6-point set in general position with ≤ 6 distinct distances is equivalent to a member of F ⊔ U. Our search uses only U, and it re-canonicalises every entry itself. This file does not prove (A). (A) is the job of the generation, affine and quartic stages and of the F/U split, which are checked elsewhere (NOTES/phase2_counts.md, NOTES/algebra_engine.md, verifier).

**Theorem C.** Assume (A). Suppose every leaf of `CoverageSearch --n 9 --U <U>` over all 7 arrangements fails the affine post-filter N6 (`affine_n.py`), for instance because there are no leaves. Then every 9-point crescent configuration contains a core whose pattern is equivalent to a member of F.

**Proof.** Let P be a 9-point crescent configuration with no core in F. We show the search has a leaf.

**Step 1 (arrangement).** The three low pairs are three distinct pairs e_1 (carrying δ_1) and {e_2, e_2'} (carrying δ_2). There are 36·C(35,2) = 21,420 such labelled structures on 9 points. `arrangements.py` enumerates all of them and computes their orbits under S_9 by brute-force canonical forms. There are **exactly 7 orbits**, of sizes 252, 1512, 1512, 3024, 3780, 3780 and 7560. Our 7 representatives (table below) hit each orbit exactly once, and so do Wang's. Matched by name, they lie in the same orbits as Wang's (checked). So we can renumber the points of P by a permutation σ, so that its low pairs are exactly those of one representative, with δ_1 on the first listed pair. Renumbering points maps cores to cores and preserves equivalence classes of patterns. So the renumbered P still has no core in F.

| name | δ1 pair | δ2 pairs | all-high prefix | orbit size |
|---|---|---|---|---|
| disjoint | 47 | 58, 69 | 1..6 | 3780 |
| triangle | 89 | 78, 79 | 1..7 | 252 |
| star | 69 | 79, 89 | 1..8 | 1512 |
| chain, δ1 middle | 79 | 69, 78 | 1..7 | 1512 |
| chain, δ1 end | 69 | 79, 78 | 1..7 | 3024 |
| P3+K2, δ1 disjoint | 58 | 69, 79 | 1..7 | 3780 |
| P3+K2, δ1 meeting | 69 | 79, 58 | 1..7 | 7560 |

(The same check for n=8 and n=7, with the representatives shifted down by 9−n, also gives exactly 7 orbits. These are used in the tests.)

**Step 2 (cores are in U).** Every 6-subset T with no low pair is a core: its pairs are all high by definition. Its points are in general position (a subset of a general-position set) and it has ≤ 6 distinct distances, namely those among δ_3..δ_8. By (A) its pattern is in F ⊔ U up to equivalence, and by hypothesis it is not in F. So **every all-high 6-subset of P has a pattern equivalent to a member of U.**

**Step 3 (the encoding of P that the search must produce).** Fix the pair order O. It lists the 15 pairs of {1..6} in lex order, then for v = 7, 8, 9 the pairs (1,v), (2,v), …, (v−1,v). Encode P as follows:
- the δ_1 pair gets label L1 and the δ_2 pairs get L2;
- each high distance gets a high label 0..5, numbered by its first occurrence along O.

So the high labels along O form a restricted-growth string. This is one fixed labelled pattern c_P.

**Step 4 (the seed is enumerated).** In every representative, {1..6} contains no low pair (column "all-high prefix"). So c_P restricted to {1..6} is a core pattern. By Step 2 it is equivalent to some u ∈ U: there are a permutation π and a label bijection with c_P(i,j) ↔ u(π(i),π(j)). The lex first-occurrence renaming of a labelled pattern depends only on its equality structure. So the lex-renamed string of c_P restricted to {1..6} equals the lex-renamed string of u∘π, which is an element of the **orbit** of u. The program enumerates every distinct element of every orbit (all 720 permutations, deduplicated) as a seed. The first 15 pairs of O are exactly the lex pairs of {1..6}, so global first-occurrence numbering agrees with the seed's own renaming. Hence the seed equal to c_P on {1..6} is one of the seeds.

**Step 5 (the branch of c_P is never cut).** The search assigns the remaining 21 pairs in the order O. A low pair gets its prescribed low label, which is the only possibility. A high pair (u,v) gets one of the t high labels already used, or the label t if t < 6. In c_P a high pair either repeats an earlier high label or introduces a new one, and by first-occurrence numbering a new one is exactly t, with t ≤ 5 because there are 6 high distances. So the branch that follows c_P exists. Part (b) shows that every pruning test is a necessary condition, i.e. c_P passes it at every prefix. At the leaf, c_P has high class sizes {3,…,8}, the leaf condition. So the search outputs at least one leaf for the arrangement of P, and that leaf (or its stored canonical representative; N6 is invariant under renumbering points) passes N6. ∎

Contrapositive: no N6-passing leaf in any of the 7 arrangements implies no 9-point crescent configuration has all its cores in U, so every one has a core in F (under (A)).

**What Theorem C relies on.**
- (A), and the U file being exactly U. The search re-canonicalises each line and reports how many lines it skipped (bad format, more than 6 labels, duplicate classes). A line with more than 6 labels cannot be a core pattern, so skipping it is harmless. A bad-format count > 0 on the real file must be investigated.
- The proofs of the necessary conditions below.
- Correct code. This is evidenced by the tests in §4. Two code paths are cross-checked: incremental vs full local tests, and orbit-set vs canonical-form membership.

---

## Part (b): search design and soundness of every pruning and symmetry step

### Design (independent of Wang's)

- **No seven-point intermediate stage.** Wang's step (i) builds 13 seven-point patterns first. We do not; we grow directly from 6-point seeds.
- **No placement-group symmetry reduction.** We use *all* orbit elements of every U-member as seeds.
- **No five-point lookup table.**
- **Own arrangement representatives.** Five of our seven differ from Wang's. They were derived by our own orbit enumeration and then matched to his orbits.
- **Own data structures.** Patterns are packed 3 bits per pair. An open-addressing hash set holds all labelled orbit elements of U (≈ |U|·720 longs). A canonical-form mode (min over all 720 permutations of the renamed string) is available, and the tests check that it gives the same results.
- **Lookahead (our own device, see N4).** After each assignment we test every all-high 6-subset with **exactly one** unassigned pair, asking whether some label for that pair puts it in U. This is a completion test on single 6-sets, generated from U's orbit set at run time. It is *not* Wang's table of five-point restrictions with injective label mapping. Disabling it (`--nolook`) gives identical results (tested).

### Symmetry breaking and why nothing is lost

1. **Point renumbering to an arrangement representative.** This is Step 1 of the proof. Every structure lies in exactly one of the 7 orbits, verified by exhaustive enumeration over all 21,420 labelled structures.
2. **First-occurrence numbering of high labels along O.** The map from labelled patterns (with high labels as abstract classes) to restricted-growth strings along O is a bijection. Every pattern is reached once, by Step 5.
3. **Nothing else.** In particular, positions 1..6 are not reduced modulo the stabiliser of the arrangement. Wang's placement trick is deliberately not used. The search is therefore redundant (each 9-point pattern can appear up to |Stab(arrangement)| times), but it is trivially complete.

### Pruning rules: each is a necessary condition for c_P at every prefix

Every rule only inspects assigned pairs. In c_P, equal labels ⇔ equal distances (low labels included), and P is in general position.

- **N1a: Lemma 2.1(i).** At a point a, at most 3 assigned pairs share a label. Otherwise 4 points lie on the circle about a with that radius, so 4 points are concyclic.
- **N1b: Lemma 2.1(ii).** For points a ≠ b and labels (u,v), at most 2 points x have c_ax = u and c_bx = v. Such x lie on two circles with distinct centres a and b, which meet in at most 2 points. (u = v is allowed.)
- **N1c: Lemma 2.1(iii).** For a ≠ b, at most 2 points x have c_ax = c_bx (the common label may vary with x). Such x lie on the perpendicular bisector of ab, so 3 of them would be collinear.
- **N1d: Lemma 2.2.** No 4 points a, b, c, d have c_ab = c_cd, c_ac = c_bd and c_ad = c_bc. Our own proof: by SSS there is an isometry φ with a↦d, b↦c, c↦b, since the triangles abc and dcb have equal corresponding sides; this holds for degenerate triangles too. φ swaps b and c, so it fixes their midpoint m. An isometry fixing m and swapping b and c is either the reflection in the perpendicular bisector ℓ of bc, or the half-turn about m.
  - *Reflection:* d is the mirror image of a in ℓ, and d ≠ a because the points are distinct. The circle through a, b, d has its centre on ℓ, since ℓ bisects both ad and bc perpendicularly. So that circle also passes through c, and a, b, c, d are concyclic. If instead a, b, d are collinear, 3 points are collinear. Either way general position fails.
  - *Half-turn:* a + d = b + c, so a b d c is a parallelogram whose diagonals are ad and bc. The third equality |ad| = |bc| makes it a rectangle, hence concyclic, or a degenerate collinear figure. Contradiction.
  
  All four local tests are monotone (adding pairs only adds configurations), so testing only the configurations that contain the newly assigned pair is equivalent to testing everything (`localInc` vs `fullLocal`; `--paranoid` asserts that they agree at every node).
- **N2: Capacity.** The 6 high classes have final sizes exactly {3,…,8}. So the map (high label ↦ final size) is an injection into {3..8}, and the current count of each label is ≤ its final size. By the elementary Hall-type lemma (Wang Lemma 5.2; our proof: with current counts sorted m_1 ≤ … ≤ m_s and capacities sorted b_1 < … < b_q, an injection with m_i ≤ b_π(i) exists iff m_i ≤ b_{q−s+i} for all i), the following is necessary: s ≤ 6 and m_i ≤ 3 + (6 − s) + (i − 1). Proof of the "only if" direction: m_i, …, m_s are s−i+1 counts, each of which needs a capacity ≥ m_i, but only s−i capacities exceed b_{q−s+i} when m_i > b_{q−s+i}.
- **N3: U-membership on completed cores.** Right after pair (u,v) is assigned, every 6-subset T = S ∪ {u,v} with S ⊆ {1..u−1} has all 15 pairs assigned. These are exactly the subsets that become complete at this step. If T is all-high (a static property, given the arrangement), its lex-renamed string must lie in the orbit set of U, by Step 2 and the argument of Step 4. Every all-high 6-subset is checked exactly once: at the step of its last pair (u,v) with u < v its two largest elements, or as the seed when T = {1..6}.
- **N4: One-missing-pair lookahead.** Right after (u,v) is assigned, take any all-high T = S ∪ {u, w, v} with S ⊆ {1..u−1}, |S| = 4 and u < w < v. Its only unassigned pair is (w,v), which is the last pair of T in lex order. In c_P, T is a core, so its renamed string is in U's orbit set. Lex renaming is prefix-stable: the first 14 renamed entries depend only on the 14 assigned labels. The last entry is either a label already among those 14 (renamed value r < nx) or a label not among them, renamed nx. The second case needs nx ≤ 5, since a core has ≤ 6 labels. The program tries r = 0..min(nx, 5), a superset of the true value. So if no option is in U, c_P is not on this branch.
- **N5: Leaf.** High class sizes must be exactly {3..8}. The low counts are 1 and 2 by construction (asserted).

Everything outside N1–N5 is bookkeeping. The `--paranoid` flag re-checks every leaf against all C(9,6) subsets with a separate loop, and checks incremental vs full local tests at every node.

- **N6: Affine post-filter on complete leaves** (`src/coverage/affine_n.py`, exact ranks with python-flint, cross-checked with Fractions via `--check`). These are Wang's parallelogram-only affine tests, which he applies after each point. We apply them once, to the complete 9-point pattern. This loses nothing, because the rows of R and A for any prefix are rows of the full R and A, so a prefix rejection implies a full rejection.
  - R has a row e_a + e_d − e_b − e_c for each ordered quadruple with c_ab = c_cd and c_ac = c_bd. By the proof of N1d, the reflection case violates general position, so a + d = b + c, i.e. ρ·p = 0 with Σρ_i = 0.
  - **T1:** a nonzero vector of rowspace(R) supported on ≤ 3 points is an affine dependency among ≤ 3 points. Support 2 means two points coincide; support 3 means three points are collinear. So we reject if rank R[:, I^c] < rank R for some triple I.
  - **T2:** Σρ_i|p_i − p_j|² = Σρ_i|p_i|² − 2p_j·Σρ_i p_i + |p_j|²Σρ_i, which does not depend on j. This gives linear equations A z = 0 in the 8 squared distances, and the true z has all z_l > 0 and pairwise distinct. So we reject if e_l or e_l − e_m lies in rowspace(A).
  
  The centroid rows (`--centroid`) are also sound, but they are not needed.

Omitting a necessary test can only increase the number of survivors, never lose a configuration. That is why N6 can be applied as a post-filter.

---

## Part (c): tests (all pass)

`tests/test_coverage.py` runs 9 tests in about 70 s with 1 worker:
1. For n = 7, 8, 9: exactly 7 orbits, the representatives hit each orbit exactly once, and the seed block spans no low pair.
2. n=8, U = the 5-point cores of the known n=8 examples. **Palásti's pattern is found**, together with the triangular-lattice n=8 example, which has the same pattern.
3. n=8 negative control: dropping one of Palásti's cores from U makes Palásti disappear. The U test really prunes.
4. n=7, U = cores of the 8 lattice examples. **All 8 known n=7 lattice examples are found.** They form 5 distinct patterns.
5. Orbit mode vs canonical mode vs `--nolook` give identical leaf counts and distinct-canonical counts.
6. n=9 planted test for the disjoint arrangement, with `--paranoid`.
7. The target patterns are crescent (class sizes 1..n−1), and the Python canonical forms are fixed points.

Longer runs, with logs in `data/`:

| Run | Result |
|---|---|
| n=7, U = cores of lattice examples, all arrangements (`coverage_n7_coresU_found.txt`) | 8/8 targets hit. 444,532 labelled leaves (≥ 40,618 canonical patterns; the disjoint count is a lower bound). This shows that 4-point cores carry little information at n=7. |
| n=8, U = cores of Palásti, all arrangements (`coverage_n8_coresU_found*.txt`) | Palásti hit (P3K2_d1_meeting). 16,992 leaves, 721 canonical. Identical output with and without lookahead (`cmp`), 85 s vs 30 s. |
| n=7, U = **all 19 local 4-point patterns with ≤ 4 labels** (the n=7 analogue of "U = all survivors"; `coverage_n7_localU.log`) | **8/8 targets hit** (star: tri7_3; chain_d1_end: tri7_0,1,6,7; P3K2_d1_meeting: tri7_2,4,5). 1,076,978,404 labelled leaves (≥ 501,644 canonical, a lower bound because the stored set is capped at 200k per arrangement; targets were matched on *every* leaf via an invariant and then the canonical form), 4.67·10⁹ nodes, 2,708 s with 1 thread |
| n=9 planted (seed 1, `--paranoid`), one planted local-passing pattern per arrangement, U = its cores (`coverage_planted_s1.log`) | **7/7 found.** The incremental and full local tests agreed at every node (≈ 50 M nodes). |
| n=9 planted + 300 random affine patterns as noise in U (seed 2; `coverage_planted_s2_noise300.log`) | **4/4 found** (disjoint 432 s, star 140 s, chain_d1_end 107 s, P3K2_d1_meeting 105 s; ≈ 2.2·10⁵ seeds each). The noise produced extra leaves only in disjoint/star, as variants of the planted pattern (distinctCanonical 3/15/2/4) |
| Python canonical forms of the targets vs the Java `--out` canonical forms | identical for all 10 targets |

Cross-implementation: arrangement orbits (Python) vs representatives (Java), Python canonical (`src/canon.py`, `targets.py`) vs Java `canonComplete`, and Java `fullLocal` vs `localInc`.

---

## Part (d): runtime

- Throughput is about 0.3–0.5 M nodes/s per thread, with about 200 nodes per seed when U is a random affine subset at n=9.
- Seeds = Σ_{u∈U} |orbit(u)| ≤ 720·|U|.
- Stand-in U, random 1000-line subset of `patterns6_affine` (716,220 seeds): 300 random lines of `patterns6_affine` (seed 7, 214,560 seeds; `coverage_standin_affine300.log`): **0 leaves in all 7 arrangements**. Disjoint: 2.1·10⁸ nodes in 430 s, levels 6:214,560 → 7:603,120 → 8:5,376 → 9:0. The other six arrangements: 4.8·10⁷ nodes each in 86–107 s, with no 7-point survivor. Total ≈ 17 min with 1 thread
- **Real U (187 patterns, 87,260 seeds):** measured at about 10 min wall-clock with 2 threads, 6.2·10⁸ nodes (disjoint 366 s, the others 30–44 s each). Canonical-form mode is several times slower. The affine post-filter takes 4 s.
- **Full affine list as U (344,835 patterns → 2.5·10⁸ seeds):** not run. The extrapolation is ≈ 200 h with 1 thread (~100 h with 2), and the orbit set would need about 4 GB, so it would have to run in `--mode canon`. This is not feasible under the CPU cap and not needed: it is only a stress test.

---

## §5. Real run

U = `data/E1_U.txt.gz`: sha256 cbd6781f94f2e2f9de47202ee7b60965e7c5f69e977789accee356f376e4214a, 187 lines. All 187 were kept (0 bad format, 0 with more than 6 labels, 0 duplicate classes after our re-canonicalisation), giving 87,260 labelled seeds. The run used 2 threads, orbit-set membership and lookahead, with a wall-clock of about 10 min. Log: `data/coverage_E1_run.log`. Survivors: `data/coverage_E1_found.txt`.

| arrangement | nodes | level counts (points 6→7→8→9) | labelled leaves | distinct 9-point patterns | after N6 affine | time |
|---|---|---|---|---|---|---|
| disjoint | 385,012,658 | 87,260 → 671,620 → 1,339,872 → 21,312 | 21,312 | 222 | **0** | 366 s |
| triangle | 41,878,992 | 87,260 → 43,540 → 100,080 → 4,320 | 4,320 | 3 | **0** | 41 s |
| star | 29,248,716 | 87,260 → 43,540 → 0 → 0 | 0 | 0 | **0** | 30 s |
| chain, δ1 middle | 41,737,440 | 87,260 → 43,540 → 100,080 → 2,880 | 2,880 | 12 | **0** | 44 s |
| chain, δ1 end | 41,734,080 | 87,260 → 43,540 → 100,080 → 2,400 | 2,400 | 20 | **0** | 44 s |
| P3+K2, δ1 disjoint | 41,820,360 | 87,260 → 43,540 → 100,080 → 6,528 | 6,528 | 68 | **0** | 43 s |
| P3+K2, δ1 meeting | 41,830,728 | 87,260 → 43,540 → 100,080 → 7,776 | 7,776 | 162 | **0** | 43 s |
| **total** | 623,262,974 | | **45,216** | **487** | **0** | |

- **The search alone (N1–N5: local, capacity and U-membership) does not finish the job.** It leaves 487 complete patterns. Wang's procedure also applies the parallelogram affine tests at each point.
- **N6 kills all 487** (`src/coverage/affine_n.py`, `data/coverage_E1_affine.txt`): 384 by T1 (a short row-space vector, i.e. a forced collinearity or coincidence) and 103 by T2 (forced z_l = 0 or z_l = z_m). The flint and Fraction rank engines agree on every one (`--check`).
- **Final: 0 complete nine-point patterns in each of the 7 arrangements**, matching Wang's Computation 5.4. By Theorem C, and conditional on (A) for this U, every 9-point crescent configuration has a core in F.
- 45,216 labelled leaves is below the 200k storage cap, so the stored survivor set is complete.
- Note: the 87,260 seeds are fewer than 187·720 = 134,640 because U-members with nontrivial automorphisms have smaller orbits.

Seven-point comparable (Wang reports 13). Every non-disjoint arrangement reaches level 7 with exactly 43,540 labelled 7-point all-high patterns, each with all of its 6-subsets in U and passing local + capacity, but with no affine tests. These 43,540 labelled patterns canonicalise (min over 7! permutations) to **exactly 13 distinct 7-point patterns** (`data/coverage_E1_level7.txt`, `--level7`, 28 s). **This equals Wang's reported 13**, even though we do not apply affine tests at this stage.

Cross-checks: (i) `--paranoid --nolook` re-run of the triangle arrangement (`data/coverage_E1_run_paranoid_spot.log`, 799 s): identical level counts 87,260 → 43,540 → 100,080 → 4,320 and 3 distinct patterns. Incremental and full local tests agreed at all 1.96·10⁸ nodes, and every leaf was re-checked against all C(9,6) subsets. (ii) An independent Python re-check of all 487 survivors and the 13 seven-point patterns with `src/canon.py` (`src/coverage/verify_survivors.py`): 0 inconsistencies, and the Java and Python canonical forms agree. (iii) A full canonical-form-mode re-run was **not completed**: it is roughly 30× slower (about 35 min for the triangle arrangement alone) and was stopped to stay inside the CPU budget. Orbit mode and canonical mode were shown to be equivalent in the n=8 pytest.

Reproduce: `bash src/coverage/run_real.sh data/E1_U.txt.gz`, then `.venv/Scripts/python src/coverage/affine_n.py data/coverage_E1_found.txt data/coverage_E1_affine.txt --check`.

