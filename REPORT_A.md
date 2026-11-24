# REPORT A — Erdős Problem #217, n = 9: independent computer-assisted verification

**Terminal state: A1 (independent confirmation).** We independently reproduced, stage by stage, the computer-assisted proof in arXiv:2611.03688 (Y. Wang, 27 Oct 2026) that **no 9-point crescent configuration exists in the plane**. Two independent implementations agree on every headline computation, and the computation that eliminates the algebraic cases also ran on a second, different engine (Gröbner bases with saturation). We found no gap and no counterexample.

Date: 2026-11-25. Folder: `open-problems/attack-217/`. Nothing has been posted or sent anywhere.

---

## 1. Statement

A *crescent configuration* is a set of n points in ℝ² with no 3 collinear and no 4 concyclic, determining exactly n−1 distinct distances, where for each i = 1,…,n−1 some distance occurs exactly i times.

**Theorem (confirmed).** There is no planar crescent configuration with n = 9.

This is the result of a computer-assisted proof. The proof structure is Wang's (§2): we did not find a new reduction. Our contribution is a complete, independent, exact re-execution of it, with:
- a second elimination engine;
- an independently designed coverage search and extension;
- an exact replay of every certificate;
- release-ready code and data.

**What is NOT claimed:**
- nothing about n ≥ 10;
- no new human-readable proof;
- no claim that C+Q = 1 is impossible (see §6).

## 2. The argument (Wang's reduction, re-derived)

Squared distances are δ1,…,δ8, where δk occurs k times. The *low pairs* are the δ1 pair and the two δ2 pairs (3 pairs). All other pairs are *high*.

1. **Cores.** Deleting one endpoint of each low pair leaves ≥ 6 points whose pairs are all high. Every 6-subset of them is a *core*: a 6-point set with ≤ 6 distinct distances (the values δ3..δ8). [Wang Lemma 3.2; re-derived by the reader and by V-b]
2. **Enumeration.** Every core realises one of the 6-point distance patterns with ≤ 6 labels that pass the local necessary conditions:
   - ≤ 3 equal distances at a point (else 4 points are concyclic);
   - two circles meet in ≤ 2 points;
   - ≤ 2 points lie on a perpendicular bisector;
   - no "parallelogram" equalities;
   - the exact affine tests.
3. **Algebraic elimination.** For each pattern, the Gram matrix of the 5 difference vectors has rank ≤ 2. Its 3×3 minors, together with the linear relations, generate an ideal. Division by the provably nonzero factors (z_l, z_l − z_t, triangle areas) and membership of a monomial or a 4-point concyclicity determinant prove non-realizability. The survivors split into F (finitely many normalised solutions, r ∈ {5,6} labels) and U (the rest).
4. **Case F.** Some core of the 9-point set is in F. Then its coordinates are reconstructed exactly in number fields K_q, and the 3 remaining points are forced, so they can be enumerated exactly:
   - trilateration or circumcentres for the ≥ 2 "high-linked" outside points;
   - the last-point lemma for the 9th point.
5. **Case U.** No core of the 9-point set is in F. Then every all-high 6-subset of it is in U. A combinatorial search over the 7 possible arrangements of the 3 low pairs, with the local, capacity and exact affine tests, shows that no 9-point distance pattern has this property.

Every pruning step is a necessary condition. The proofs are in `NOTES/coverage.md` (coverage, necessity of rules N1–N6) and `NOTES/extension.md` (Lemmas T, 7.2, 7.3, capacity). V-b's red-team review, which re-derived the 7 arrangements from 21,420 low-pair structures, is in `NOTES/verification.md`.

## 3. Results, stage by stage

| Stage | Wang | Ours (implementation 1) | Independent check (implementation 2) |
|---|---|---|---|
| 3/4/5-point patterns | 3 / 22 / 932 | 3 / 22 / 932 | V-a enumerator: identical |
| 6-point, local tests | 625,472 | 625,472 (Java A ≡ Java B brute force over 446M labelled) | V-a: identical set |
| after affine tests | 344,835 | 344,835 (flint ≡ Fraction) | V-a: identical set |
| excluded by quartic/division (E1) | 340,198 | 340,198, each with a certificate replayed over Q | V-a replay: 340,198/340,198 ✔; **E2 Gröbner: 340,198/340,198 EXCLUDED** |
| retained = F + U | 4,637 = 4,450 + 187 | 4,450 + 187 | V-a classifier: identical on all 4,637 |
| U by labels r = 4/5/6 | 3 / 9 / 175 | 3 / 9 / 175 | identical |
| irreducible factors with a real root | 4,119 | 4,119 (E1 route); 4,118 (E2 route, see note) | — |
| realised cores (factors / embeddings) | n/a | 3,879 / 7,829 | V-b: identical per factor |
| 9-point configurations from F-cores | 0 | 0 (funnel 74,382 → 4,489 → 1,468 → 163,807 → 0) | V-b, arb balls + exact: 0 |
| 7-point stage of the U-search | 13 | 13 | V-b: same 13 |
| complete 9-point patterns with all cores in U | 0 | 487 before the 9-point affine test, 0 after | V-b: same 487, same 384/103 kill split, 0 |

**The 4,118 vs 4,119 factor count is resolved.** The two routes differ on a single pattern, 112233145565546, where the E1 multiplication-matrix route has one extra factor with a real root. The extra factor is non-realizable and is removed by E2's stronger saturation. Both routes give 3,879 realised cores, so nothing downstream changes.

## 4. Independence and verification

- **Rule 3 (checker first).**
  - Two exact checkers: number-field/flint, and Q(√3)/rational with circumcentre and cross-ratio concyclicity.
  - Both accept Palásti's n = 8 set and reject broken variants.
  - They agree on 57 examples and about 1.2k perturbations.
  - Exhaustive radius-3 lattice search: k = 6/7/8/9 → 26 / 8 / 1 / 0, identical in both. The published solutions of sallerk/erdos-notes reduce to the same k = 4..8 counts.
- **Rule 4 (two implementations).**
  - The verifiers V-a and V-b wrote their code without reading `src/`.
  - The extension engines differ in design: the builder uses the high-linked lemma and a covering family; V-b uses neither and splits by real embedding.
  - The coverage searches differ in design: the builder grows patterns with first-occurrence label breaking and an orbit-set; V-b uses its own representatives, a five-point lookup and stabiliser orbits.
- **Exactness.** Every exclusion has an exact certificate. Every "keep" decision is the safe direction. Real-root and sign questions use Sturm counts and arb balls with exact fallback; no float decides anything.
- **Soundness tests.**
  - Realisable patterns from lattices, Palásti, C3 orbits and random triangular sets are never excluded: 145 + 127 canonical classes. A further 252 non-lattice patterns over fields of degree ≤ 12, taken from the 879 exact 9-point sets, are also all in F ⊔ U.
  - 581/581 known F-realisations are recovered exactly.
  - Palásti n = 8 is recovered by the extension code.
  - Planted 9-point patterns are found in 7/7 arrangements.
  - Known n = 7 and n = 8 examples are found by the coverage code.
- **Prior-work re-check, 2026-11-25** (`NOTES/prior_work.md`):
  - No one else has independently confirmed or refuted n = 9.
  - Wang's paper is still v1, with 0 citations and no code.
  - #217 is still listed as open (teorth/erdosproblems, formal-conjectures `217.lean`).

## 5. Caveats and limitations (stated honestly)

1. **Same reduction as Wang.** Our independence lies in the implementation and the engines, not in the mathematics of the reduction. A conceptual error in the reduction itself would affect both. We re-derived every lemma (reader, V-a, V-b) and found none.
2. **Corollary 4.5** ("every realisable core is in the computed F ⊔ U") is not proved separately. It follows from two facts:
   - (i) every exclusion step is a sound necessary condition: Lemmas 2.1–2.3 and 4.1–4.2, proved on paper and re-derived by three of us;
   - (ii) the enumeration is complete.
   Both facts were checked by two independent implementations, and every exclusion carries a certificate checked by two independent methods. The empirical spot checks (§4) are consistent with this, but they are evidence, not proof.
3. **Six F patterns** (out of 4,450) hit the Gröbner limit and were reconstructed only by the multiplication-matrix route, which is complete by Wang's Lemma 6.1. V-b's extension check took its cores from those reconstructions. Its 4,492 factors match the E2-route count, but V-b did not independently re-reconstruct those six patterns.
4. **Computer-assisted proof.** Trust rests on the code, which is released (`release/`), plus agreement between independent implementations. The proof is not formally verified: there is no Lean proof.
5. **One builder did both elimination engines.** E1 and E2 were written by the same builder (A). Their independence comes from the method (linear-algebra quartic spans vs Gröbner saturation in flint), and V-a independently replayed E1's certificates.

## 6. Side finding: near-miss configurations (`NOTES/near_misses.md`)

The forum thread for #217 (IheartTesla, 4 Sep 2026) asked how small C + Q can be: C collinear triples plus Q concyclic quadruples, among 9-point sets with multiplicities exactly 1..8. Their best was (C, Q) = (1, 2) on the triangular lattice, which we confirmed exactly.

- **We found (C, Q) = (0, 2) over Q(√21)**, confirmed by three independent exact codes. The squared distances, scaled so that δ8 = 1, are:

  | δ1 | δ2 | δ3 | δ4 | δ5 | δ6 | δ7 | δ8 |
  |---|---|---|---|---|---|---|---|
  | (−1+√21)/2 | (3+3√21)/10 | (−3+√21)/2 | (5−√21)/2 | (9−√21)/10 | (−4+√21)/5 | (1+√21)/10 | 1 |

  The exact coordinates are in the note and in `verify/orch/nearmiss_config2.py`.
- Among all 879 exact 9-point sets with multiplicities 1..8 produced by the F-core extension, the minimum C + Q is 2. This is not a global search, so C + Q = 1 is not ruled out.

## 7. Reproduction

Requirements: Python ≥ 3.12 with `sympy python-flint numpy networkx pytest`, and Java 21. Peak RAM was about 8 GB, for Java enumeration. All runtimes were measured on a 20-core Windows 11 machine, and most later runs were capped at ≤ 8 worker processes.

| Step | Command (run from the repo root) | Runtime |
|---|---|---|
| checker tests | `python -m pytest tests/test_checker.py -q` | 2 s |
| enumeration | `javac -d src/enum/build src/enum/PatternEnum.java && java -Xmx8g -cp src/enum/build PatternEnum A data 4 6` | 3 s (16 thr) |
| affine tests | `python src/affine.py data/patterns6_local.txt.gz data/patterns6_affine.txt.gz 4 flint` | 11 s |
| E1 (all 344,835) | `python src/algebra_full.py --workers 4 --chunk 5000`, then extract F/U (see `NOTES/algebra_engine.md`) | ≈ 2.5–3 h, 4 workers |
| E2 (all) | `bash src/run_e2_all.sh`, then the retry on unresolved patterns with `ALGEBRA_E2_LIMITS=100000,1000000,1000000` | ≈ 2.5 h + 2 h, 4 workers |
| F reconstruction + extension | `python src/extend_run_recon.py …`, `python src/extend_run.py …` (see `NOTES/extension.md`) | ≈ 7 min + 1 min, 2 workers |
| coverage (case U) | `bash src/coverage/run_real.sh data/E1_U.txt.gz`, then `python src/coverage/affine_n.py data/coverage_E1_found.txt data/coverage_E1_affine.txt --check` | ≈ 10 min, 2 threads |
| V-a | `python verify/va/run_all.py` | ≈ 3.6 h, 1 worker |
| V-b | `python verify/vb/coverage_vb.py`, `python verify/vb/extend_rig_vb.py data/F_recon_final.jsonl.gz verify/vb/ext_rig` | ≈ 10 min + 33 min |

The total compute was about 16 CPU-hours of wall-clock jobs, within the 24 h budget. SHA-256 hashes of the data products are in `data/MANIFEST.md`. Key hashes:
- `E1_F.txt.gz`: 6a6b6e13…f9fe
- `E1_U.txt.gz`: cbd6781f…214a
- `E1_excluded_certs.jsonl.gz`: 7f57361d…e8fe
- `E2_retry/out.jsonl`: 404e6f63…b779

## 8. Relation to arXiv:2611.03688

- **Every numerical claim we could test reproduces exactly:** Table 1, the 13 seven-point patterns, the 4,119 factors, the maximum order 12, and 0 completions.
- **Wang released no code or data.** The completeness of his U-search and the soundness of his elimination were the main unverifiable points. Both are now independently reproduced with exact certificates.
- **Minor corrections to the paper's presentation:**
  - one endpoint of each low pair is deleted (not "the endpoints", as some summaries say);
  - the figure has a mis-reference;
  - the acknowledgement mentions AI-assisted "constructions" that the paper does not contain.
  None of these affects correctness.

**Confidence: high** that n = 9 is impossible:
- every stage reproduces Wang's counts exactly;
- every exclusion carries two independent certificates;
- the two case searches are independent implementations that agree to the pattern;
- extensive known-answer tests pass.

The residual risk is a shared conceptual error in the reduction. We re-derived all of its lemmas and found none.

## 9. Files

- `NOTES/`: prior_work, wang_method, phase2_counts, algebra_engine, coverage, extension, verification, near_misses.
- `src/`, `tests/`: the builders' code.
- `verify/`: the verifiers' independent code (`checker_v.py`, `va/`, `vb/`, `orch/`).
- `data/`: enumerations, verdicts, certificates and MANIFEST.
- `release/`: the publishable package, a local git repository with one commit and not pushed.
- `OUTBOX/`: an optional draft note to the author, NOT sent.
