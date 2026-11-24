# LOG (append-only)

## 2026-11-24 — cycle 1 (Phase 0 start)
- Read PROMPT.md. Created NOTES/ src/ tests/ data/ raw/ OUTBOX/, STATE.md, LOG.md.
- Created venv `.venv` with sympy, python-flint, numpy, networkx, pytest (pip, rule 7 OK).
- Briefed scout (-> NOTES/prior_work.md) and reader (-> NOTES/wang_method.md) in parallel.

## 2026-11-24 — cycle 2 (Phase 0 scout done; Phase 1 start)
- Scout -> NOTES/prior_work.md. n=9 NOT independently confirmed/refuted; no public code. Goal A stands (no re-scope).
- PROMPT table corrections: Burt "900 hours" (not CPU-hours); Wang runtime "~4 h, 16 GB, 24 threads" (not garbled); formal-conjectures has 217.lean (open, added 2026-10-23). Independence caveat: Wang co-authors 1610.07836v3; Chen same institution; Wang acknowledgment mentions an LLM assisting "search for constructions" (odd for nonexistence paper).
- Palásti n=8 coords from Burt Fig.1 match PROMPT; scout's ad-hoc check passes (to be redone by official checker). n=7 coords not found online.
- Launched builder(i) [checker + tests] and verifier [independent checker from spec] in parallel.

## 2026-11-24 — cycle 3 (reader done; Phase 2 launched)
- Reader -> NOTES/wang_method.md; TeX source in raw/2611.03688/src/. Paper exists as described; no code/data. Correction: deletes ONE endpoint of each of the 3 low pairs. Counts: 3/22/932/625,472 -> 344,835 (affine) -> 4,637 = F 4,450 + U 187; 13 seven-pt residual patterns; 0 nine-pt completions; 4,119 factors. Reader reproduced 3/22/932. No logical gap on paper; risk is all in unreleased computation (upstream drops of realizable patterns are the critical failure mode).
- Tooling: no C compiler; Java 21 + Node available; WSL Ubuntu has no gcc/msolve. Plan: Java for enumeration, python-flint for exact algebra.
- Launched builder E (6-pt enumerator + affine tests, reproduce counts, realizable-pattern soundness suite -> NOTES/phase2_counts.md) and builder A (Engine 1 = Wang quartic-span re-impl, Engine 2 = Gröbner/saturation, independent engine -> NOTES/algebra_engine.md).
- Independent route chosen (Phase 3 option b+c): different elimination engine (Gröbner w/ saturation) + independent re-implementation of coverage search.

## 2026-11-24 — cycle 4 (verifier checker done)
- Verifier independent checker verify/checker_v.py (Q(√3) pairs + rational; concyclic via circumcentre AND cross-ratio, never 4x4 det): 18/18 checks PASS. Palásti n=8 (as printed in 1509.07220) ACCEPTED exactly (sizes 1..7), also after rotations/mirror; broken variants REJECTED for the right reasons.
- Exhaustive radius-3 hex patch (37 pts): n=6 26 configs, n=7 8 configs (432 sols), n=8 1 (Palásti), n=9 none. Examples in verify/lattice_n7_r3.txt. NEW VERIFIED FACT.
- Pending: cross-check vs builder src/checker.py.

## 2026-11-24 19:23 — CPU cap applied (user request via parent session)
- User: CPU at 100 °C, throttling. ~18 busy python procs were running. Killed (attack-217 processes only): builder E affine.py Fraction cross-check (8-pool; flint result data/patterns6_affine.txt.gz already saved), builder E soundness_run W12 (8-pool; W8 run already done), builder(i) lattice_search R4 (optional). Left running: gen.py (1 proc).
- Standing rule recorded in STATE.md: max 8 concurrent worker processes total. Budgets sent to agents: E ≤5, A ≤2, (i) ≤1.

## 2026-11-24 — cycle 5 (Phase 1 closed)
- Builder(i) src/checker.py (number field Q[t]/(m), flint): pytest 57 passed. Accepts Palásti n=8, n=3,4 examples; rejects broken. Exhaustive radius-3: k=3..9 -> 24,36,19,26,8,1,0 (agrees with verifier on 6..9). data/known_examples.json.
- Cross-check verify/crosscheck_checkers.py: 57/57 accepted by both, 0 disagreements on perturbed variants. Rule 3 satisfied; Phase 1 DONE.
- 19:2x: Builder E had relaunched affine.py Fraction cross-check with a 12-worker pool at 19:23:33 (before the cap message reached it). Killed tree PID 43012. Sent E a hard-cap message: ≤5 workers, explicit worker counts, script defaults 4 and never os.cpu_count().

## 2026-11-24 — cycle 6 (builder A partial handback)
- src/algebra.py: Engine 1 (Wang quartic-span, exact over Q; mod-p hint only, exclusions replayed over Q with certificates + verify_certificate()) and Engine 2 (flint Buchberger over Q, saturation; EXCLUDED/ZERODIM/POSDIM/UNRESOLVED). Reconstruction (Sturm + arb isolation) implemented: 10 F test patterns -> 7 exact GP realizations, 2 fail PD, 1 no real root.
- E1 on 1500 sampled local patterns: 45% affine-excl (paper 45%), retained 0.67% (paper 0.74%). E1 median 0.023 s, p99 0.69 s -> ~9 CPU-h for 345k. E1 vs E2: no inconsistency on ~80 patterns (theory: E1 ideal ⊆ E2 saturated ideal).
- OPEN: soundness on realizable sets not yet run (cache data/algebra_realizable.json building; compare run data/algebra_compare_run1 running, 1 worker each).

## 2026-11-24 — cycle 7 (Phase 2 counts reproduced)
- Builder E: 3 / 22 / 932 / 625,472 / 344,835 — ALL MATCH Wang Table 1 exactly. By r (6-pt): local 275/18,938/166,912/439,347 (r=3..6); affine 4/2,191/60,895/281,745. Two Java generators (extension + brute force over ~446M labelled) byte-identical; flint vs Fraction affine engines byte-identical. Soundness: ~137,700 exact GP point sets -> 204 distinct patterns, all retained (weak: lattice-type only). NOTES/phase2_counts.md, hashes in data/MANIFEST.md.
- Next: E1 full run over 344,835 (target 4,637 = 4,450 F + 187 U); E2 independent confirmation; coverage search builder C.
- cycle 8: messaged builder A (E1 full run, ≤5 workers); launched builder C (independent coverage search, ≤2 workers).

## 2026-11-24 22:45 — cycle 9 (E1 full run collected; resumed after usage-limit pause)
- Builder A hit usage limit; its E1 background run completed all 69 chunks. Orchestrator extracted: 344,835 records = affine set exactly; EXCLUDED 340,198 (all quartic stage, every certificate replayed over Q OK, 0 mod-p discrepancies); F 4,450 (r5 210, r6 4,240); U 187 (r4 3, r5 9, r6 175). **EXACT MATCH with Wang Table 1 (4,637 / 4,450 / 187 = 3/9/175).**
- Soundness: 2,024 labelled realizable patterns (145 canonical classes; Palásti, tri/square lattices, two-triangle, C4-orbit sets): E1 verdicts U/F only (1 per-pattern UNRESOLVED), E2 POSDIM/ZERODIM/UNRESOLVED only — 0 EXCLUDED. All 145 canonical classes ∈ E1 F∪U.
- Launched src/run_e2_all.sh: Engine 2 (Gröbner, independent engine) over all 344,835 in 35 resumable chunks (4 workers, 120 s timeout; driver clamp now env ALGEBRA_MAX_WORKERS, hard max 5). Resumed builder C with real U (≤2 threads). Total busy ≤ 6.
- cycle 10: launched builder D (F-core reconstruction via E2 eliminants + extension, ≤2 workers).

## 2026-11-24 23:2x — cycle 11 (builder D: F-branch closed, 0 found)
- 4,450 F patterns reconstructed: 4,444 via E2 Gröbner eliminants (independent engine), 6 via Wang multiplication matrices (Gröbner limit) -> those 6 single-engine for now. 4,492 irreducible factors; 4,118 with real root (Wang: 4,119 — 1 apart, OPEN, E1 cross-check running ~4 h, 507/4,450 matched so far); 3,879 realizable-core factors (3,845 patterns, 7,829 PD embeddings). Max h=12 (= Wang).
- Extension (Lemmas 7.1–7.3 re-derived, NOTES/extension.md): 74,382 -> 4,489 (7-pt) -> 1,468 (8-pt) -> 163,807 ninth-pt candidates -> 0 nine-point crescents. All-triangles rerun identical. 0 UNRESOLVED.
- Near misses: 5 nine-point sets (2 up to similarity, quadratic fields) with class sizes exactly 1..8 and no 3 collinear, failing only on a concyclic quadruple. Recorded as a finding.
- Soundness: 581/581 known F realizations recovered; Palásti n=8 recovered by the extension code from its all-high 6-subsets; capacity test = brute force on 1,500 instances. Gap: |D(S)|=5 high-linked path not exercised end to end by a known config.
- Builder C real coverage run partially seen: disjoint 222 / triangle 3 / star 0 distinct 9-pt patterns before the 9-pt affine stage; 0 after. Awaiting report.

## 2026-11-24 — cycle 12 (builder C coverage done; verifiers launched)
- Builder C real run (U = E1_U, 187; 87,260 seeds; 6.23e8 nodes; ~10 min, 2 threads): 487 distinct complete 9-pt patterns before affine (disjoint 222, triangle 3, star 0, chain-mid 12, chain-end 20, P3K2-disj 68, P3K2-meet 162); ALL 487 fail exact parallelogram affine tests (384 row-space, 103 linear-distance; flint = Fraction) -> 0. Seven-point stage = 13 patterns (= Wang). Tests: n=7 all 8 found with U=all local 4-pt patterns; n=8 Palásti found; planted n=9 7/7 and 4/4 with noise. Coverage proof + necessity of N1–N6 in NOTES/coverage.md (relies on Corollary 4.5 for computed lists).
- Launched verifier V-a (independent enum + independent certificate replay + independent quartic classifier on 4,637 retained) and V-b (red-team coverage/extension proofs, independent coverage search, independent extension on reconstructed cores, near-miss check). 1 worker each. CPU: E2 4 + D-crosscheck 2 + Va 1 + Vb 1 = 8.

## 2026-11-25 — cycle 13 (V-b report)
- V-b (no src/ access): red team found no gap; arrangements (7 classes, re-derived from 21,420 structures), high-linked lemma, |D(S)|=5 case, last-point lemma, capacity, symmetry breaking all CONFIRMED. Cor 4.5 (computed F⊔U ⊇ all realizable cores) PARTIAL: 127 lattice/C3/random-tri patterns + 252 non-lattice patterns (from 8,523 GP 6-subsets of 879 exact 9-pt sets, fields of degree ≤12) all ∈ F⊔U.
- Independent coverage search: 487 identical per arrangement; 384/103 affine kill split identical; 13 seven-pt identical; per-step affine also 0. CONFIRMED.
- Independent extension, full run (967 s): 4,492 factors / 3,879 cores / 7,829 embeddings identical per factor; 0 crescents. CAVEAT: float64 candidate rejection (margin 2.4e-12 vs 1.85e-6) -> asked V-b to redo with arb balls + exact fallback.
- Near misses CONFIRMED exactly: 2 nine-point sets with sizes 1..8, no 3 collinear, failing only by concyclic quadruples (config 1 over q = u² − 6u + 25/3).

## 2026-11-25 — cycle 14 (E2 full run collected)
- E2 (flint Buchberger, saturation; independent elimination engine) over all 344,835: E1-excluded -> E2 EXCLUDED 335,275, UNRESOLVED 4,923, ZERODIM/POSDIM 0 (no disagreement). E1-F -> ZERODIM 4,422 + UNRESOLVED 28. E1-U -> POSDIM 178, ZERODIM 2, UNRESOLVED 7. Cert kinds: GB(I)={1} 197,046; product-of-factors-in-I remainder.
- data/E2_unresolved_excluded.txt (4,923). Retrying with limits 1e5/1e6/1e6, timeout 600 s, 4 workers -> data/E2_retry/out.jsonl. V-a asked to replay all 4,923 certificates independently.
- 2026-11-25 cycle 15: final prior-work re-check (rule 5): nothing changed; Wang still v1, 0 citations; #217 open in teorth + formal-conjectures; no one else confirmed/refuted n=9. NOTES/prior_work.md.

## 2026-11-25 — cycle 16 (forum thread input + near-miss finding)
- Parent session relayed user screenshots of erdosproblems forum thread 217: sallerk lattice search (diam² ≤ 400, no n=9; code github sallerk/erdos-notes p217), IheartTesla (C,Q)=(1,2) lattice near-miss, Alfaiz note. Added to NOTES/prior_work.md.
- Verified exactly: IheartTesla config (C,Q) = (1,2) ✔. Our config 2 over Q(√21): sizes 1..8, (C,Q) = (0,2), confirmed by 3 independent exact codes (builder D, V-b, orchestrator sympy). Min C+Q over V-b's 879 exact 9-pt sets = 2; none ≤ 1 (scope: F-core extensions only). NOTES/near_misses.md, verify/orch/*.py.
- cycle 17: scout read sallerk p217 (read-only): their published solutions reduced to our radius-3 patch match k=4..8 = 36/19/26/8/1 exactly; their exhaustive diam²≤400 n=9 = 0 subsumes our patch result. NOTES/prior_work.md.
- cycle 18: V-b rigorous rerun (arb balls 128-bit + exact NF fallback, no float decisions; 1,979 s): funnel identical (429,734 candidates, 879 exact 9-pt sets sizes 1..8), 0 crescents. Case 1 (F-branch) now has two independent exact implementations.

## 2026-11-25 — SCOPE CHANGE (user via parent session)
- Finish goal A only; do NOT start B (n=10). After REPORT_A.md: build attack-217/release/ (README, src/tests/verify, small data + MANIFEST, near_misses.md, MIT LICENSE, no venv/abs paths/personal info), git init + ONE local commit there, NO remote/push (parent confirms push with user). Then STATE terminal.

## 2026-11-25 03:5x — cycle 20 (all verification in; A1 reached)
- E2 retry: 4,923/4,923 EXCLUDED -> E2 independently certifies all 340,198 exclusions.
- V-a: enumeration 625,472/344,835 set-identical; all 340,198 certificates replayed True (6,914 s); classifier agrees on all 4,637 (F 4,450; U 3/9/175).
- Factor discrepancy resolved: E1 route 4,493 / 4,119 real (= Wang), E2 route 4,492 / 4,118; differs only on 112233145565546 (extra non-realizable component removed by E2 saturation); realized cores 3,879 both.
- NOTES/verification.md signed off. Terminal state A1. Writing REPORT_A.md + release/.

## 2026-11-25 — cycle 21 (CLOSE: terminal A1)
- REPORT_A.md written; OUTBOX/draft_note_to_author.md (not sent). STATE.md set TERMINAL (A1). Goal B not started (user scope change).
- release/ built: src, tests, verify, data (files >1 MB gzipped, total ~41 MB, none >45 MB), NOTES, REPORT_A, README, MIT LICENSE, .gitignore, SHA256SUMS. Excluded: .venv, raw/ (arXiv source not redistributed), OUTBOX, build artefacts. Scrubbed absolute JDK paths (now JAVA_HOME/PATH). No personal info found (grep incl. inside .gz). git init + one local commit with neutral identity; no remote, no push.
- cycle 21b: V-a final report: all CONFIRMED (+ classifier excludes 13,442/13,442 sampled excluded). Addendum in NOTES/verification.md; release refreshed (amended single local commit).
