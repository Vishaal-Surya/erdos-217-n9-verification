# Wang, arXiv:2611.03688 — method spec (reader notes)

Read 2026-11-24 from the arXiv TeX source. Every statement below comes from the source; my own checks are tagged.
Tags: **[verified]** = I checked the logic or number myself. **[unclear]** = underspecified, or I could not check it. **[suspicious]** = a reason for doubt is given.

## 0. Provenance

- Abstract page: title "The nonexistence of nine-point crescent configurations", author Yifei Wang (Nankai University, LPMC). v1 only, submitted Tue 27 Oct 2026 01:19 UTC (29 KB). Subjects math.MG, math.CO. MSC 52C10, 05C30, 14Q10, 68R05. **The paper exists and matches the description.** [verified]
- Files in `raw/2611.03688/`: `abs.html`, `eprint.bin` (gzip tarball, sha256 ee40307a…eb715), and `src/` containing `crescent9_v4__1_.tex` (1567 lines, sha256 e912b71f…d5ad), `crescent9_v2.bib`, `00README.json`. The tarball was extracted with Python tarfile using `filter='data'`. Nothing in it was executed.
- The file names (`v4__1_` for the TeX, `v2` for the bib) point to several private drafts. The gzip header mtime is 4 Nov 2026, but this is probably arXiv repackaging and is not evidence of a v2.
- **No code, data or certificates are released.** No availability statement appears anywhere. [verified: absent from the source]
- Acknowledgment: "OpenAI's GPT-6-Astra model assisted in the search for the constructions presented in this paper." The paper presents no constructions, so the wording is odd. [suspicious: AI-assisted, wording mismatched]
- Corrections to PROMPT §2:
  - (a) The method deletes **one endpoint of each low pair**, not "the endpoints".
  - (b) The runtime is not garbled in the source: "approximately four hours on a machine with 16 GB of RAM configured for 24 threads", wall-clock.

## 1. Definitions

- **General position:** no 3 points collinear, no 4 concyclic. **Crescent configuration:** a set of n points in general position with exactly n−1 distances, one of multiplicity k for each 1≤k≤n−1. Multiplicity counts unordered pairs, and the distances are not ordered by size. [verified, matches the problem]
- Squared distances are used throughout. δ_k is the squared distance occurring k times (k=1..8 for n=9). The δ_k are pairwise distinct and positive, with no ordering assumed.
- **Low pairs** carry δ1 or δ2. There are exactly 3 of them: one δ1 pair and two δ2 pairs. **High pairs** carry δ3…δ8.
- **Distance pattern** on n indices: a partition of the C(n,2) pairs into r nonempty classes (labels). Equal label means equal squared distance. Distinct labels mean distinct squared distances. **Equivalence** = permuting the points plus bijectively renaming the labels. **Realization** = distinct points p_i with |p_i−p_j|² = z_{c_ij} for some pairwise-distinct positive z.
- Formal matrix D: D_ii = 0, D_ij = z_{c_ij}.
- **High-distance core** (Def 3.1): a 6-point subset S ⊂ P with D(S) ⊆ {δ3,…,δ8}. This is any such 6-subset, not only the ones produced by the deletion.
- "Multiplicity classes": the paper never relates a core's internal class sizes to 3..8. Inside a core the multiplicities are free, and a high distance may occur once in S. [verified]

## 2. Reduction to the 6-point core (Lemma 3.2)

- Delete **one endpoint of each of the 3 low pairs**, which removes at most 3 points. Every pair among the remaining ≥6 points is high. Any 6 of the remaining points form a core. [verified]
- **≤6 distances:** D(S) ⊆ {δ3..δ8}, so S has at most 6 distinct distances (r ≤ 6). [verified]
- **Disjoint low pairs:** exactly 3 points are deleted and exactly 6 remain. There are 2³ = 8 choices, so a 9-set has exactly 8 cores. [verified]
- **Low pairs share an endpoint:** deleting the shared endpoint plus one endpoint of the third pair removes ≤2 points. This leaves a 7-point (or larger) set whose pairs are all high. For the star arrangement (19, 29, 39), points 1..8 are all high. Any 6 of these points form a core, so the core is **not unique**.
- The core is always exactly 6 points by definition. The 7-point all-high set is used only in the coverage step (§6). [verified]

## 3. Pattern enumeration (Procedure 3.3, Table 1)

- **Object:** distance patterns on 6 points with ≤6 labels, of arbitrary class sizes, taken up to equivalence (S_6 acting on points × renaming of labels). [verified as stated]
- **Generation:** start from the 2-point pattern. For n=3..6, extend each canonical (n−1)-point representative by a point n. Label (1,n)…(n−1,n) in order, each with an existing label or the next new label (t+1, only if t<6). Partial assignments are rejected on:
  - Lemma 2.1 (i): at most 3 pairs at a point share a label. Four equidistant neighbours would be concyclic.
  - Lemma 2.1 (ii): for points a,b and labels (u,v), at most 2 points x have c_ax=u and c_bx=v. Two circles with distinct centres meet in at most 2 points.
  - Lemma 2.1 (iii): for points a,b, at most 2 points x have c_ax = c_bx. Three such x would all lie on the perpendicular bisector, so they would be collinear.
  - Lemma 2.2: c_ab=c_cd, c_ac=c_bd **and** c_ad=c_bc is impossible.
  
  All four tests are necessary conditions. [verified: proofs checked, including the parallelogram computation s=−u, t=±v, where t=v gives an isosceles trapezoid (concyclic) or collinear points, and t=−v gives a+d=b+c]
- **Canonical form:** each point gets a signature s_i = sorted multiset of (m_ℓ, d_iℓ), where m_ℓ is the class size and d_iℓ the number of label-ℓ pairs at point i. Group the points by signature and order the groups. Over all orderings within the groups, take the lexicographically least first-occurrence-renamed label sequence in pair order (1,2),(1,3),…,(n−1,n). This is a valid canonical form because the signatures are invariant under equivalence. [verified]
- Completeness (Prop 3.4): delete a point, induct, and use the fact that restrictions of general-position sets stay in general position with ≤6 labels. [verified]
- **Affine tests at 6 points (Lemma 2.3):**
  - Build R from rows p_a+p_d−p_b−p_c, one for each ordered quadruple with c_ab=c_cd and c_ac=c_bd (the parallelogram relation, Lemma 2.2).
  - Add a row p_a+p_b+p_c−3p_d whenever abc is equilateral and d is equidistant from a, b, c (d is the centroid).
  - **Reject** if the row space contains a nonzero vector with ≤3 nonzero entries. Test: for some triple I, rank R_{I^c} < rank R. [verified]
  - **Linear distance equations:** for each row ρ and each j=2..n, Σ_i ρ_i (D_ij − D_i1) = 0. This gives Az=0. Take a basis N of ker A. **Reject** if some column of N is zero (forces z_i=0) or two columns are equal (forces z_i=z_j). [verified: Σρ_i|p_i−p_j|² is independent of j]
- **Counts (Table 1):**

  | Stage | Patterns |
  |---|---|
  | 3-point | 3 |
  | 4-point | 22 |
  | 5-point | 932 |
  | 6-point (local + parallelogram tests) | 625,472 |
  | after the 6-point affine tests | 344,835 |
  | after the quartic tests (F ⊔ U) | 4,637 |
  | F | 4,450 |
  | U | 187 (3 four-label, 9 five-label, 175 six-label) |

  - **I reproduced 3, 22 and 932 by brute force** (all set partitions with ≤6 blocks, filtered by the local and parallelogram tests on every subset, canonicalised by full permutation). (Scratch script not included; superseded by src/enum and verify/va.) [verified]
  - 625,472 / 344,835 / 4,637 / 4,450 / 187 were not checked. [unclear; this is Phase 2 work]
  - 4450 + 187 = 4637 and 3 + 9 + 175 = 187. [verified]

## 4. Algebraic conditions (§4)

For each 6-point pattern that passes, with r ≤ 6 labels:

- **Gram rank-2 condition:** L_ij = D_1i + D_1j − D_ij for 2≤i,j≤6, a 5×5 matrix. Its value is 2·Gram(p_i − p_1), so rank L ≤ 2.
  - All 3×3 minors (principal and nonprincipal) are homogeneous cubics in z.
  - Up to the transpose symmetry there are 55 of them: 100 minors in total, 10 fixed by the symmetry, (100−10)/2 + 10 = 55. [verified]
- **Working space:** V = Q[z]_4 (homogeneous quartics), with dim V = C(r+3,4) ≤ C(9,4) = 126. [verified]
- **Initial W:** the span of {z_ℓ · minor} ∪ {linear relation from Az=0 · cubic monomial}.
- **Division (Lemma 4.1):** the nonzero factors f are
  - z_ℓ,
  - z_ℓ − z_t,
  - 4D_ij D_ik − (D_ij + D_ik − D_jk)² for i<j<k, which equals 4·det(p_j − p_i, p_k − p_i)² > 0. [verified]
  
  Set C_f = {q ∈ Q[z]_{4−deg f} : f q ∈ W} and W += span{m·q : q ∈ C_f, m a monomial of degree deg f}. Repeat full passes to a fixpoint (≤126 enlargements). Sound because f ≠ 0 at a realization. [verified]
- **Exclusions (Lemma 4.2):** reject the pattern if W contains
  - (a) a degree-4 monomial, or
  - (b) det(D restricted to any 4 points). Here (D_ij) = B·J·Bᵀ with rows of B equal to (1, x, y, x²+y²) and J = [[0,0,0,1],[0,−2,0,0],[0,0,−2,0],[1,0,0,0]], so det = −4(det B)², which is nonzero for 4 points in general position. [verified, including det J = −4]
- **Finite criterion:** V = W + z_a·Q[z]_3 for some label a, i.e. the monomials divisible by z_a span V/W. Let h = dim V/W.
- **Degree reduction (Lemma 4.3):** with z_a = 1, every polynomial reduces modulo ⟨W|_{z_a=1}⟩ to a combination of b_1…b_h (deg ≤ 3). So the quotient ring is ≤ h-dimensional and there are finitely many normalized tuples. [verified]

## 5. F/U split (Procedure 4.4)

- **F** = retained patterns with r ∈ {5,6} that satisfy the finite criterion for some a. The smallest such a and the final W are stored.
- **U** = all other retained patterns: every pattern with r ≤ 4 (even if it is finite), plus the 5- and 6-label patterns that fail the criterion for every a.
- The restriction r ∈ {5,6} exists because the extension lemma (§7) needs |D(S)| ≥ 5. [verified as a design rationale]
- **U is not classified geometrically.** It is handled only through the overlap/coverage argument of §6.
- Corollary 4.5: every core pattern is equivalent to a member of F ⊔ U. This follows because every test above is necessary. [verified, conditional on a correct implementation]

## 6. Completeness argument (Prop 5.1: "every 9-point crescent has a core in F")

**Argument (paraphrase).** Suppose no core is in F. Then every core is in U. Every 6-subset with only high pairs is a core, so every such restriction of the global 9-point pattern is equivalent to a U-member. Search for any 9-point pattern with this property:

- **Low-pair arrangements:** 3 edges with no isolated vertex give 5 graphs (triangle, star K_{1,3}, path P4, P3+K2, 3K2). Marking the δ1 edge gives 7 cases:

  | Arrangement | δ1 | δ2 |
  |---|---|---|
  | disjoint | 47 | 58, 69 |
  | triangle | 89 | 18, 19 |
  | star | 19 | 29, 39 |
  | chain, δ1 middle | 89 | 18, 29 |
  | chain, δ1 end | 18 | 89, 29 |
  | P3+K2, δ1 disjoint | 39 | 18, 28 |
  | P3+K2, δ1 meeting | 18 | 28, 39 |

  [verified: the enumeration is complete and the representatives are correct; in the disjoint row points 1..6 are all high, and in every other row points 1..7 are all high]
- **(i) Seven-point patterns:** extend each U-member by one point (same labeling rule). Keep it only if all seven 6-subsets are equivalent to U-members, it passes the local and parallelogram tests, the affine tests with parallelogram relations only (row-space and linear-distance tests), and the capacity test with capacities {3..8}. Canonicalise. **Reported result: 13 seven-point patterns.** [unclear: not checked]
- **(ii) Placement:** m=6 for the disjoint row (place every U pattern), m=7 otherwise (place every one of the 13).
  - Symmetry reduction: positions 1..m are grouped by identical relations (δ1/δ2/high) to every later position, and only increasing placements within a group are kept. This preserves the low pairs, so it is sound. [verified]
  - Duplicates are then removed after renaming high labels by first occurrence.
- **(iii) Extension:** add positions m+1..9.
  - Low pairs are prescribed. The δ1/δ2 labels are fixed and cannot be reused.
  - Every other new pair, in increasing order of its old endpoint, gets an existing high label or the next unused one, with ≤6 high labels in total.
  - After each assignment: the local tests, the capacity test (Lemma 5.2: sorted counts m_i ≤ b_{q−s+i}, proof checked [verified]), and a U-membership test for every fully assigned all-high 6-subset.
  - After all pairs to the new point are assigned: the affine tests (parallelogram rows only).
  - At 9 points: low counts must be 1 and 2, and high counts must be exactly {3,4,5,6,7,8}.
- **Five-point lookup acceleration:**
  - The table is indexed by ordered 5-point patterns and lists, over all orderings of all U patterns, the five labels joining them to the sixth point.
  - When adding a point: choose 5 old points whose mutual pairs and pairs to the new point are all high, and look the pattern up.
  - Labels on the 10 old pairs map through the existing correspondence. Remaining entry labels map injectively to high labels absent from those 10 pairs, with ≤6 overall.
  - Sound, because the 6-set is a core and must be in U. [verified logically]
- **Reported result (Computation 5.4): 0 complete nine-point patterns in all 7 arrangements.**

**Assessment of each step.**

1. Any 6-subset with only high pairs is a core. This holds by the definition. [verified]
2. "All cores in U" follows from "no core in F" via Corollary 4.5. [verified]
3. The 7 arrangements are exhaustive. [verified]
4. The non-disjoint cases have a 7-point all-high subset. [verified]
5. Every test in (i)–(iii) is a necessary condition for a real general-position 9-set: the local tests, the parallelogram exclusion, parallelogram affine rows (including the linear-distance test with low labels included in z), capacity, and U-membership. [verified]
6. The placement symmetry reduction is sound. [verified]
7. The label-assignment branching covers every labeling, with high labels renamed by first occurrence and low labels fixed. [verified]
8. **The whole weight rests on the two unreleased computations:** the correctness of the 187-pattern list U (and so the correctness of every exclusion and of the finite criterion upstream), and the 13 / 0 results. If any genuinely realizable core pattern was wrongly excluded upstream (bug in the minors, the affine tests, the division or the canonical form), it is in neither F nor U. The coverage search would then wrongly prune 9-point patterns containing it. **The argument is logically sound but entirely computational. No proof is independent of the code.** [suspicious only in the sense of being unverifiable; no logical gap found]
9. The logic does not need U-patterns to be non-realizable or infinite: U is "everything else". This makes the split robust. A finite but misclassified pattern costs nothing logically. [verified]
10. The dependency figure (Fig. 1) cites "Proposition~\ref{sec:coverage}", a section label, and contains typos ("Definiton", "lie in"). These are cosmetic. [suspicious: sloppiness only]

## 7. Core reconstruction and extension

**Multiplication matrices (§6).**

- For a pattern in F with z_a = 1: the basis B_1..B_h of V/W uses monomials divisible by z_a, chosen greedily by degree after normalization and then by lex order on the degree-4 exponent vector. b_i = B_i|_{z_a=1}.
- The coefficient column of a monomial m (deg ≤ 4) is the coordinates of [z_a^{4−deg m}·m] in that basis.
- e = column of 1. M_a = I. For j ≠ a, column k of M_j is the column of z_j·b_k.
- The evaluation functional ℓ(v) = Σ v_k b_k(z) satisfies ℓ(e) = 1 and ℓ(M_j v) = z_j ℓ(v). [verified]
- **Lemma 6.1:**
  - Take M = Σα_j M_j such that f = det(uI − M) is squarefree and T = (e, Me, …, M^{h−1}e) is invertible.
  - Define g_j by its coefficient vector T^{−1}M_j e.
  - Then every realization has λ = Σα_j z_j with f(λ) = 0 and z_j = g_j(λ). So there are at most h tuples. [verified]
  - Note: z_j = g_j(λ) needs only g_j(M)e = M_j e (true by construction) plus ℓ. Commutativity of the M_j is used only for g_j(M) = M_j. Soundness does not depend on it.
- **Choice of M (App. A.2):** for c = 0, 1, 2, run M ← cM + M_j over j ≠ a in increasing order. Take the first squarefree candidate, then check T and g_j(M) = M_j.
  - With c = 0 the result is just the last M_j. With c = 1 it is the sum. With c = 2 it is the weights 2^k.
  - If all three fail, the case would be "unresolved". The paper says all pass. [unclear: no fallback is described]
- **Reported:** all of F passes, with order h ≤ 12, and **4119 irreducible factors** are retained in total (factors with a real root, with z_j ≠ 0 and z_j ≠ z_k in K_q). [unclear]

**Reconstruction (§6.1).**

- Work in K_q = Q[u]/(q). Use the frame p1 = (0,0), p2 = (1,0), p3 = (0,1) with metric G = [[D12, (D12+D13−D23)/2], [·, D13]].
- Discard if det G = 0. Otherwise require some real root with D12 > 0 and det G > 0, checked by exact isolating intervals.
- Solve the 2×2 system (coefficient matrix 2G) for p4, p5, p6. Verify all 15 distances.
- Check general position: det(p_j − p_i, p_k − p_i) ≠ 0, and the 4×4 det with rows (1, pᵀ, pᵀGp) ≠ 0.
- [verified]

**Extension (§7).**

- **Trilateration (Lemma 7.1):** for noncollinear a, b, c and squared distances α, β, γ in K, the unique candidate solves 2(b−a)ᵀGx = bᵀGb − aᵀGa + α − β and 2(c−a)ᵀGx = cᵀGc − aᵀGa + α − γ. Accept it if ‖x−a‖²_G = α. Circumcenter: both differences set to 0. [verified]
- **High-linked point:** an outside point with ≤1 low pair to S.
- **Lemma 7.2:** ≥2 of the 3 outside points are high-linked. Otherwise ≥4 low pairs would exist. [verified]
  - If |D(S)| = 6: all of its ≥5 high pairs carry core values. Using the triangles (1,2,3) and (4,5,6) is enough, since the single low pair spoils at most one of them. That is 2(6³+1) = 434 trials. [verified]
  - If |D(S)| = 5: either ≥3 pairs carry core values (trilateration), or ≥3 carry the missing high value (circumcenter). Use all 20 triangles: 20·(5³+1) = 2520 trials. [verified]
- **Lemma 7.3 (last point):** for any 8-subset Q, either some value absent from D(Q) has multiplicity ≥3 (then x is a circumcenter of 3 of its partners), or only δ1/δ2 are absent (≤3 pairs, so ≥5 old values and trilateration). The construction uses all 56 triangles × D(Q)³ plus circumcenters. [verified]
- Consequence: all points lie in K². One search in K_q covers every embedding with G positive definite, because the tests are only zero/equality tests in a field. [verified]
- **Partial tests (7- and 8-point sets):** there must exist an injective μ: D(Q) → {1..8} with:
  - m_Q(d) ≤ μ(d),
  - μ ≥ 3 on D(S),
  - each selected outside point having ≥5 pairs to S with μ ≥ 3.
  
  How this is checked:
  - With 6 core distances: core counts must fit {3..8}, noncore counts must fit {1,2}, and there must be ≥5 core values from each outside point.
  - With 5 core distances: for each k ∈ {3..8} (the missing distance's multiplicity), core counts must fit {3..8}\{k}. Then either the missing value is absent (noncore counts fit {1,2}), or it is some d0 with m_Q(d0) ≤ k (other noncore counts fit {1,2}). Each outside point needs ≥5 high-valued pairs, with the same k and d0 for all outside points.
  
  [verified: necessary conditions]
- **Final (9-point) test:** general position, 8 distinct values with multiplicities exactly 1..8, and every core value occurring ≥3 times. **Reported: no configuration survives.** [unclear]

## 8. Reported computations and availability

- **Table 1 counts.** Seven-point residual patterns: 13. Retained 9-point patterns: 0 in each of the 7 arrangements. Max matrix order: 12. Retained irreducible factors: 4119.
- **Not reported:** the number of factors surviving reconstruction, the number of extension candidates, the degree distribution of q, and the per-stage timings.
- **Runtime:** ~4 h wall-clock, 16 GB RAM, 24 threads (single run).
- Implementation details: Appendix A.1 expands the 55 minors once in 15 pair variables, then substitutes labels. No language, libraries or CAS are named.
- **No code, data, U-list, F-list or certificates are released, and there is no availability statement.** [verified absent]

## Most likely failure points (ranked)

1. **Upstream soundness bugs that silently drop realizable core patterns.** These would sit in the affine, Gram-minor, division or quartic-exclusion code, or in canonicalisation, and would push a pattern out of F ⊔ U. Coverage (§6) then uses U-membership as a *pruning* test, so a missing pattern prunes real 9-point patterns. This is the single point where a bug turns into a false theorem, and none of it can be checked without code. It is reproducible: recount 625,472 → 344,835 → 4,637 (4,450 / 187) independently.
2. **The residual search (Computation 5.4: 13 seven-point patterns, 0 completions).** This is a complex branch-and-prune with a five-point lookup table, injective label mapping and placement symmetry reduction. The logic checks out, but off-by-one or symmetry bugs are typical here. It needs independent re-implementation, with the U list rebuilt from scratch.
3. **The finite-core extension (Computation 7.6).** The search covers 4119 factors, number-field arithmetic in K_q, candidate generation, and the capacity tests with the subtle 5-distance case (choice of k and d0). The risks are a mishandled |D(S)| = 5 case, using only triangles (1,2,3)/(4,5,6) when |D(S)| = 6 (sound only if the high-linked logic is implemented exactly), and the duplicate and core-point removal.
4. **Missing M-choice fallback.** If some pattern in F had failed the squarefree/T tests, the case would be unresolved. The paper asserts that none failed, but this is unverifiable.
5. **Presentation red flags.** These do not affect the logic: an AI-assisted acknowledgment that mentions "constructions" the paper does not contain, a mis-referenced figure, typos, draft-numbered source files, a single author and no peer review.

**No logical gap was found in the written argument.** The reduction, the F/U split and the coverage logic all check out on paper. The open risk is entirely computational and unreleased.
