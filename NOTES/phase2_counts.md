# Phase 2 — independent reproduction of Wang's core-pattern counts (Builder E)

Date: 2026-11-24. Spec: arXiv:2611.03688 TeX (Lemmas 2.1–2.3 = `lem:local`, `lem:parallelogram`, `lem:affine`; Procedure 3.3 = `proc:generation`; Table 1). Read as data only.

## Result

| Stage | Wang (Table 1) | Builder E | Agree |
|---|---|---|---|
| 3-point patterns | 3 | 3 | yes |
| 4-point patterns | 22 | 22 | yes |
| 5-point patterns | 932 | 932 | yes |
| 6-point patterns (local + Lemma 2.2) | 625,472 | 625,472 | yes |
| after the 6-point affine tests | 344,835 | 344,835 | yes |

**No discrepancy.** Every count matches exactly, from two independent generators (Java A and B) and two independent exact rank engines (flint, Fraction).

Counts by the number of labels r. The paper does not give this breakdown.

| k | r=1 | r=2 | r=3 | r=4 | r=5 | r=6 |
|---|---|---|---|---|---|---|
| 3 | 1 | 1 | 1 | – | – | – |
| 4 | 0 | 4 | 8 | 7 | 2 | 1 |
| 5 | 0 | 2 | 49 | 260 | 385 | 236 |
| 6, local | 0 | 0 | 275 | 18,938 | 166,912 | 439,347 |
| 6, after affine | 0 | 0 | 4 | 2,191 | 60,895 | 281,745 |

Affine rejections (first failing test, in the order T1 then T2-zero then T2-equal):
- 122,907 by the row-space support-≤3 test (T1).
- 92,231 by a forced z_l = 0.
- 65,499 by a forced z_l = z_m.
- Total rejected: 280,637 = 625,472 − 344,835.

## Exact filter definitions as implemented

The pattern on k points is a label c_ij for each pair i<j. The string is written in lex pair order (1,2),(1,3),…,(1,k),(2,3),…,(k−1,k), with labels as digits 1..6 renamed by first occurrence.

**Local tests.** A pattern is kept only if it passes all of these. Each test is monotone, meaning a violation inside any subset is a violation of the whole pattern. So "every sub-pattern passes" is the same as "the full pattern passes", and incremental generation with partial rejection, as Wang does it, gives the same set.
- **L(i), Lemma 2.1(i):** for every point a and label ℓ, #{x : c_ax = ℓ} ≤ 3.
- **L(ii), Lemma 2.1(ii):** for every a<b and every ordered label pair (u,v), #{x ∉ {a,b} : c_ax=u, c_bx=v} ≤ 2.
  - u=v is included. This is immaterial: the u=v case is implied by L(iii). The TeX says "prescribed labels u,v" with no distinctness required.
  - Ordered (u,v) relative to the ordered pair (a,b) is the geometrically correct reading: two circles with centres a,b and radii u,v.
- **L(iii), Lemma 2.1(iii):** for every a<b, #{x ∉ {a,b} : c_ax = c_bx} ≤ 2. The common label may vary with x, as the TeX says explicitly.
- **P, Lemma 2.2:** for no 4-set {a,b,c,d} do all three perfect matchings carry equal labels, i.e. c_ab=c_cd, c_ac=c_bd and c_ad=c_bc. The labels of different matchings may coincide; the K4-all-equal case is rejected here.
- At most 6 labels in total.

**Equivalence / canonical form (our own, not Wang's).** The canonical form is the minimum over **all k!** point permutations of the first-occurrence-renamed lex string. No invariant shortcut is used. Wang's signature-group ordering gives a different representative of each class, but the class counts must agree, and they do.

**Affine tests (6 points, Lemma 2.3).** These use exact integer matrices with rank over ℚ.
- **R** contains:
  - the rows e_a+e_d−e_b−e_c for every ordered quadruple (a,b,c,d) of distinct points with c_ab=c_cd and c_ac=c_bd (the parallelogram a+d=b+c);
  - the rows e_a+e_b+e_c−3e_d whenever c_ab=c_ac=c_bc and c_ad=c_bd=c_cd (centroid).
  - Duplicate rows are dropped.
- **T1:** reject if some triple I has rank R[:,I^c] < rank R. This holds exactly when rowspace(R) contains a nonzero vector supported in I. Supports of size 2 are covered by the triples that contain them.
- **T2:** A is the coefficient matrix in z of Σ_i ρ_i(D_ij − D_i1) = 0, taken over each row ρ of R and each j=2..6, with D_ii=0 and D_ij=z_{c_ij}.
  - Reject if e_l ∈ rowspace(A) for some l. This is "z_l=0 forced", the same as a zero column of a kernel basis N.
  - Reject if e_l − e_m ∈ rowspace(A) for some l<m. This is "z_l=z_m forced", the same as equal columns l,m of N.
  - Membership is tested as rank([A; v]) = rank(A).
- **Basis independence of "equal columns of N".** The rows of any basis N of ker A span ker A. So column l of N is zero iff every kernel vector has z_l=0. Columns l and m are equal iff every kernel vector has z_l=z_m, i.e. iff z_l−z_m ∈ (ker A)^⊥ = rowspace(A). So the test does not depend on which basis is chosen, and our basis-free form is equivalent to Wang's.

## Implementations

- `src/enum/PatternEnum.java` gives two independent generators that share only the test function and the canonicaliser.
  - **A (extension):** extend each canonical (k−1)-point rep by point k with every label choice (existing labels plus one new one, at most 6 in total), test, canonicalise and dedupe. This is Wang's scheme.
  - **B (orderly brute force):** depth-first search over all restricted-growth label strings on the C(k,2) pairs, in colex pair order with test pruning. A full labelling is kept iff its own lex string equals its canonical form.
  - The two outputs are **byte-identical** for k=3..6. B visits 446,390,770 labelled passing 6-point patterns, and exactly 625,472 of them are canonical.
- `src/canon.py` is an independent Python canonicaliser and local-test implementation.
- `src/affine.py` runs the affine tests with two rank engines: python-flint `fmpz_mat.rank` (production) and pure-Python `Fraction` Gaussian elimination (cross-check).
- `tests/test_phase2.py` (pytest, 10 tests, about 30 s):
  - Python brute force for k=3,4,5 gives the same sets as the Java files.
  - The Java canonical forms are fixed points of the Python canonicaliser and survive random permutation and relabelling (sample of 1,500).
  - The flint and Fraction engines agree on a sample.
  - Survivors are a subset of the local list.
  - Hand examples for L(i), L(iii) and P.
  - On every stored realizable example, R annihilates the exact coordinates, A·z=0 holds at the true squared distances, and the pattern is in the affine list.
  - Palásti's 6-subsets survive.
- `tests/soundness_run.py` runs the soundness test described below.

Reproduce:
```
JDK="$JAVA_HOME/bin"   # any Java 21 JDK
"$JDK/javac" -d src/enum/build src/enum/PatternEnum.java
"$JDK/java" -Xmx8g -cp src/enum/build PatternEnum A data 16 6      # writes data/patterns{3,4,5}_local.txt, patterns6_local.txt.gz
"$JDK/java" -Xmx8g -cp src/enum/build PatternEnum B <outdir> 16 6  # cross-check generator
.venv/Scripts/python src/affine.py data/patterns6_local.txt.gz data/patterns6_affine.txt.gz 16 [flint|fraction]
.venv/Scripts/python tests/soundness_run.py 8 16
.venv/Scripts/python -m pytest tests/test_phase2.py -q
```

## Runtimes (20-core Windows box)

| Step | Wall-clock | Workers |
|---|---|---|
| Java A, k=3..6 (writes the local lists) | 2.9 s | 16 threads |
| Java B, orderly brute force k=6 (cross-check) | 260 s | 16 threads |
| affine.py, flint engine, 625,472 patterns | 10.9 s | 16 procs |
| affine.py, Fraction engine, full cross-check | 135 s | 5 procs |
| soundness W=8, 16 lattice forms + Palásti | 33 s | 12 procs |
| soundness W=10, 4 main lattice forms | 66 s | 4 procs |
| pytest tests/test_phase2.py | about 25 s | 1 |

The runs above were made before the coordinator's CPU cap. All scripts now default to 4 workers.

Full cross-checks:
- Java A vs Java B: the two patterns6 lists are **byte-identical** (B: 446,390,770 labelled passing 6-point patterns → 625,472 canonical).
- flint vs Fraction affine output: **byte-identical**, 344,835 survivors (`data/affine_fraction_crosscheck.log`).

## Soundness test (mandatory)

Realizable 6-point configurations come from exact point sets. Every one satisfies two conditions, both checked in integer arithmetic:
- **General position:** no collinear triple (integer cross product) and no concyclic quadruple (4×4 squared-distance determinant ≠ 0, Wang eq. (circlefactorization)).
- **At most 6 distinct distances.**

Sources:
- **Lattices ℤ² with the Euclidean metric given by a positive definite integral form** Q = a dx² + b dx dy + c dy². Every such lattice embeds isometrically in ℝ², so its point sets are genuine.
  - (1,0,1) is the square lattice.
  - (1,1,1) is the triangular lattice ((2a+b, b√3)/2, squared distance a²+ab+b²).
  - (3,0,1) is the lattice {(u√3/2, v/2)}, which contains Palásti's 8 points.
  - 13 further forms were also used.
  - Translations are normalised by putting the lex-min point at (0,0). Points range over a box [0,W]×[−W,W].
- **Palásti's 8-point set.** It verifies as a crescent configuration (multiplicities 1..7, general position). We took all of its 6-subsets with at most 6 distances in general position: 15 sets.

| Run | Raw realizable 6-sets | Distinct canonical patterns | Missing from local | Missing from affine |
|---|---|---|---|---|
| W=8, 16 forms + Palásti | 44,445 | 202 (r=4: 1, r=5: 19, r=6: 182) | 0 | 0 |
| W=10, forms (1,1,1),(1,0,1),(3,0,1),(1,0,3) + Palásti | 93,245 | 194 | 0 | 0 |
| union | 137,690 (overlapping) | **204** (r=4: 1, r=5: 19, r=6: 184) | **0** | **0** |

**Result: no misses.** Every realizable pattern found is in both `patterns6_local` and `patterns6_affine`.

A pytest check confirms more for every stored example:
- R annihilates the exact coordinates.
- A·z = 0 at the true squared distances.
- The Python canonical form equals the stored one and lies in the affine list.

This is a direct check that R and A are built correctly, not only that the filter passes.

There is a limitation: the realizable sample is small, 204 classes out of 344,835 survivors, and all of it is lattice-type, rich in parallelograms. That is exactly the regime that exercises the affine filters. Outputs are `data/soundness_phase2_*.json` and `data/soundness_realizable6_*.txt.gz` (pattern, source, points).

## Hashes

See `data/MANIFEST.md`.

## Caveats

- These counts reproduce Wang's generation and affine stages only. The quartic/division stage (→ 4,637 = F 4,450 + U 187) is not covered here.
- The soundness test can only exercise patterns that are realizable on lattices and on Palásti's set. It cannot prove that the filters are sound. Soundness rests on the proofs of Lemmas 2.1–2.3, which were checked in NOTES/wang_method.md. The test is evidence against implementation bugs, both ours and, by the count match, Wang's.
