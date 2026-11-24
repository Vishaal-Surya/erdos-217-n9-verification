# Builder A: algebraic elimination / certificate engine (status 2026-11-24)

Code: `src/algebra.py` (engines), `src/algebra_testgen.py` (exact realizable generators),
`src/algebra_batch.py` (killable workers + timeout, MAX_WORKERS=2 because of the CPU cap),
`src/algebra_compare.py` (E1 vs E2 agreement + timing), `tests/test_algebra.py`.

## Engine 1 (Wang §4, re-implemented)
- Lemma 2.3: R = parallelogram rows + centroid rows. Rejected if some triple I has rank R_{I^c} < rank R. Az=0 is built from every row and every j. Rejected if e_l or e_l−e_t lies in rowspace(A).
- W0 = z_l·(all 3×3 minors of L, 55 distinct) + (rref rows of A)·cubic monomials, inside V = Q[z]_4.
- Division fixpoint by z_l, z_l−z_t and the triangle factors. C_f is the left kernel of [f·m_i ; W].
  The closure does not depend on the order of the divisions (it is the least closed superspace).
- Exclusion: W has an rref row with a single nonzero entry (a monomial), or some 4-point det(D) reduces to 0 mod W.
- Finite criterion: rank of W restricted to the monomials not divisible by z_a equals C(r+2,4). F requires r∈{5,6}; everything else retained is U.
- **Exactness.** All verdicts are exact over Q (flint fmpz_mat: fraction-free rref and integer nullspace).
  A mod-p run (p = 2^62−57, nmod_mat) is only a hint. If it excludes, its enlargement trace is replayed over Q, and if the replay does not exclude, the full Q fixpoint runs. If it keeps, the full Q fixpoint always runs. Any mod-p/Q mismatch is flagged as `modp_discrepancy`.
- **Certificate** for a quartic exclusion: `{"trace": [(factor name, dim before, dim after), ...], "target": [monomial|det4, ...]}`. `verify_certificate()` replays it over Q.
  Affine exclusions carry their reason, e.g. `affine-support(I)`, `lin-zero zl`, `lin-equal`.

## Engine 2 (independent: Gröbner over Q, flint `fmpz_mpoly_vec.buchberger_naive`, degrevlex)
- Linear relations are eliminated by exact rref substitution, and one free z is set to 1. Generators: the 3×3 minors.
- Saturation uses a single Rabinowitsch variable t with P = NF(∏ factors). The product covers z_l, z_l−z_t, the 20 triangle factors and the 15 4-point dets, and is reduced modulo GB(I) one factor at a time. This leaves the ideal I+<1−tP> unchanged.
- EXCLUDED when any of these holds, each with its certificate kind:
  - `factor-zero`: a nonzero factor is identically 0 under the linear relations.
  - `GB1/I`: GB(I) = {1}.
  - `NF0`: the product of the factors lies in I.
  - `GB1/sat`: GB(I+<1−tP>) = {1}.
- ZERODIM: every variable has a pure-power leading monomial. Each z is then reported with its minimal polynomial in Q[z,t]/J, i.e. the eliminant of the saturation, computed by exact Q division and fmpq_mat minpoly.
- POSDIM otherwise.
- UNRESOLVED when the Buchberger limits are hit (10000, 100000, 100000) or the per-pattern timeout fires. UNRESOLVED is never counted as an exclusion.
- Mathematically, E1 exclusion ⇒ E2 EXCLUDED (W ⊆ saturation), so "E1 excludes and E2 reports ZERODIM/POSDIM" would be a bug.

## Finite-core reconstruction (Wang §6, §6.1): implemented
- `f_univariate`: basis B_i, M_j, commutativity check, the c=0,1,2 choice of M, f, T, g_j, and a check of g_j(M)=M_j.
- `reconstruct_F`: number field K_q = Q[u]/q. The real roots are rigorous: an exact Sturm count plus arb isolation. Positive-definiteness of G is checked, p4..p6 are solved, all 15 distances are verified, and general position is checked.
- Test on the 10 F patterns found in a sample of 1500:
  - 7 contain a factor with an exact general-position realization (e.g. 111231234315412 has a rational z ratio of 13/7).
  - 2 fail positive-definiteness, and 1 has no real root.
  - The E2 eliminants match E1's f whenever the normalisation agrees (for example 111232433553364: z6 → x^3−168x+54 in both engines).

## Results so far
- Engine 1 on 1500 random patterns from `data/patterns6_local.txt.gz` (other builder): 675 affine-excluded (45%, the paper has 45%), 815 quartic-excluded, 10 retained (all F; 0.67%, the paper has 0.74%).
  No mod-p discrepancy. Time per pattern: median 0.023 s, p99 0.69 s, max 4.6 s, mean 0.053 s.
  Estimate for 345k affine-passing patterns: about 8–9 CPU-h, about 4.5 h on 2 workers.
  f_univariate+reconstruction takes about 4 s per F pattern (unoptimised), so about 5 CPU-h for 4.45k F patterns.
- Engine 2: 0.01–0.3 s typical, with occasional patterns of 2 s+ at the larger limits.
  In the first smoke test of 40 patterns: every E1 exclusion was matched by E2 EXCLUDED, except 5 that hit the original tight limits (since raised); at the raised limits all 5 were EXCLUDED. No inconsistency was seen.
- Realizable sets (first generator run: 137 canonical classes from Palásti subsets, triangular and square lattice subsets, two concentric equilateral triangles, and 90°-rotation orbits): all pass the local tests. Engine-level soundness on these had **not yet been run when this note was written**; see the pending items.

## Pending (runs were in progress at handback)
- `data/algebra_realizable.json`: a larger cached realizable set, being built by the generator (slow, pure Fraction arithmetic).
- `data/algebra_compare_run1.{jsonl,_summary.json,log}`: E1 vs E2 on 2000 random local-passing strings + 2000 sampled patterns, 1 worker, 60 s timeout.
- `pytest tests/test_algebra.py` has not yet been run end to end. It needs the cache file. Without the file it builds it, which takes 30+ min.
