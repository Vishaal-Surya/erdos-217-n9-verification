# Near-miss 9-point configurations (finding, 2026-11-25)

**Question (forum thread 217, IheartTesla, 4 Sep 2026):** among 9-point sets with distance multiplicities exactly 1..8, how low can C+Q go (C = collinear triples, Q = concyclic quadruples)? Their best was (C,Q) = (1,2) on the triangular lattice, and they asked whether 2 or 1 is reachable.

## Results
1. **Forum config confirmed exactly** (`verify/orch/iheart_check.py`, integer arithmetic): lattice points (4,-2),(4,-3),(2,1),(3,-1),(2,2),(-1,0),(0,0),(1,-1),(3,1) have class sizes 1..8, C = 1 {(1,2,3)}, Q = 2 {(0,2,3,4),(4,6,7,8)}.
2. **C+Q = 2 is attained with C = 0** ("config 2", over Q(√21)). It was found by builder D's extension of F-core 111232433454125, with q = u² − 17/5·u + 1, u = (17 − 3√21)/10.
   - It was confirmed by THREE independent exact implementations: builder D (K_q arithmetic), verifier V-b (own number-field code, `verify/vb/nearmiss_vb.json`), and the orchestrator (`verify/orch/nearmiss_config2.py`, sympy in Q(√21), G-metric lift).
   - Class sizes are 1..8, 0 collinear triples, and exactly 2 concyclic quadruples: (0,1,7,8) and (1,2,6,7) in the point order of `verify/vb/ext_full.jsonl`.
   - Exact squared distances (scaled so δ8 = 1):

     | multiplicity | squared distance |
     |---|---|
     | 1 | (−1+√21)/2 |
     | 2 | (3+3√21)/10 |
     | 3 | (−3+√21)/2 |
     | 4 | (5−√21)/2 |
     | 5 | (9−√21)/10 |
     | 6 | (−4+√21)/5 |
     | 7 | (1+√21)/10 |
     | 8 | 1 |

   - Euclidean coordinates (numeric; exact in the frame plus Cholesky of G, see script): (0,0), (0.747166,0), (0.451555,0.595278), (0.669195,0.332319), (0,0.664637), (0.373583,0.927597), (0.077971,−0.332319), (−0.217640,0.262960), (−0.077971,0.996956).
3. **No C+Q = 1 among our exact sets.** Over all 879 exact 9-point sets with sizes 1..8 produced by V-b's full extension run, the (C,Q) minimum is C+Q = 2: 23 records with (0,2), 1 with (1,2), and none with C+Q ≤ 1. Several (0,2) records come from different cores; some share distance values over Q(√21). The 2nd near-miss family (builder D's config 1, over Q(√6), q = u² − 6u + 25/3) has (0,3).
   - **Scope caveat:** these sets are only those that contain a 6-point all-high core in Wang's family F, i.e. extensions of finite cores. This is NOT a search over all 9-point sets, so "C+Q = 1 is impossible" is NOT claimed. (If the n=9 nonexistence result holds, C+Q = 0 is impossible.)
