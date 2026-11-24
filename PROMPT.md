# GOAL PROMPT — Erdős Problem #217: crescent configurations (n = 9 verification, then n = 10)

You are the **orchestrator** of a small research team (you + subagents). Work autonomously toward the goal below, across as many turns and context windows as needed, until you reach one of the **terminal states**. You may stop only at a terminal state or when a **stop rule** fires. Your working folder is `open-problems/attack-217/`. Everything you learn, decide, compute or reject goes in files there, never only in chat.

---

## 1. The problem (exact)

A set of n points in ℝ² is a **crescent configuration** if
1. **general position:** no 3 points are collinear and no 4 are concyclic, and
2. it determines **exactly n−1 distinct distances**, and for each i = 1,…,n−1 some distance occurs **exactly i times** (the C(n,2) = 1+2+…+(n−1) pairs are split into classes of sizes 1,…,n−1).

Erdős asked for which n these exist (Problem #217, erdosproblems.com/217). Erdős conjectured that they fail for all large n.

## 2. What is already known (verified 2026-11-24; re-check in Phase 0)

| Fact | Source |
|---|---|
| Planar constructions exist for n = 3…8 (Pomerance n=5; Liu/Palásti n=7; Palásti n=8, on the triangular lattice) | arXiv:1509.07220 (Burt et al.); Palásti 1987/1989 |
| Palásti n=8 coordinates as printed (UNVERIFIED): (0,1), (√3,0), (2√3,0), (5√3/2,5/2), (3√3/2,9/2), (√3/2,7/2), (3√3/2,7/2), (√3,2) | arXiv:1509.07220, Fig. 1 caption |
| Exhaustive search of a 91-point hexagonal patch of the triangular lattice found no n=9 configuration (~900 CPU-hours) | arXiv:1509.07220, Remark 3.1 |
| Classification of n=4,5,6 by distance-equality patterns (5: 31 realizable strict classes; 6: 6 realizable). Its "strictness" filter can DISCARD realizable classes, so it is not a complete classification | arXiv:1610.07836 (v3, Oct 2026) |
| **CLAIM: no planar 9-point crescent configuration exists.** Exact computation by Yifei Wang. Method: delete the endpoints of the multiplicity-1 and multiplicity-2 distance classes to get a 6-point "high-distance core" with ≤ 6 distinct distances; enumerate core equality patterns; derive necessary polynomial conditions (rank-2 Gram 3×3 minors, affine relations, quartic forms of dimension ≤ 126); split the patterns into a finite-candidate family F and an unclassified residual U; reconstruct the cores via multiplication matrices of order ≤ 12; extend each by ≤ 3 points with exact multiplicity and general-position tests. ~4 h on 24 threads (garbled in the source). **Single author, v1 only, no code released, not peer reviewed.** The completeness step ("every configuration has at least one core in F") is the critical point to check. | arXiv:2611.03688 (27 Oct 2026) |
| Paper says it does NOT settle n ≥ 10: deleting a point does not preserve the multiplicities | same |
| 3D: crescent configurations exist in ℝ³ for every n (helix construction); this does not decide the plane | arXiv:2610.38368 (Chen 2026) |
| Community DB (teorth/erdosproblems, synced 2026-11-23) still lists #217 as open; erdosproblems.com froze statuses and proof claims on 2026-10-06 | github.com/teorth/erdosproblems |

## 3. Goal and terminal states

**Primary goal (A): independently settle whether a planar 9-point crescent configuration exists, with an exact, reproducible certificate, and so confirm or refute arXiv:2611.03688.**

Terminal states for A (reaching one completes A):
- **A1 — Independent confirmation.** Your own exhaustive argument, at least partly independent of Wang's (different reduction, a different core choice, or a different elimination engine), proves that no 9-point configuration exists. Every exclusion has an exact certificate, and the coverage/completeness argument is written out and checked by the verifier.
- **A2 — Counterexample.** An explicit 9-point configuration, with coordinates as exact algebraic numbers (minimal polynomials + isolating intervals) or exact radicals, that passes the exact checker. This would refute Wang. Treat it as extraordinary: the checker run must be independently re-implemented by the verifier before you report.
- **A3 — Concrete gap.** A specific, reproducible flaw in Wang's argument (for example, a 6-point core pattern that the method misses, or a pattern in U that is realizable and extends), demonstrated exactly. Report it as a gap, not as a refutation, unless A2 also holds.

**Stretch goal (B), only after A is terminal: n = 10.** Adapt the pipeline: 45 pairs, classes of sizes 1…9. Terminal states B1/B2/B3 are defined the same way. Partial results that are worth recording: exhaustive exclusion of a well-defined sub-family (for example, all configurations on a stated lattice patch, or all configurations whose core lies in a stated family), with an exact certificate.

**Honest-failure terminal state (F):** a stop rule fires. Then write the final report: what is proven, what is computed, where it is stuck, and the most promising next step. "Unresolved" is an acceptable outcome. A shaky claim is not.

## 4. Non-negotiable rules

1. **Exact arithmetic for every claim.** Floats, numerical optimisation and random search are only for finding candidates and hints. Every "exists" needs exact coordinates plus an exact check. Every "does not exist" needs an exact certificate: a Gröbner basis containing 1, a resultant/rank certificate over ℚ, a CAD cell decision, or an exhaustive enumeration with an exact per-case test.
2. **Coverage must be proved, not assumed.** Any reduction ("every configuration contains a core of type X", "WLOG the multiplicity-1 distance is …") gets a written proof in `NOTES/coverage.md` that the verifier signs off.
3. **Verify the checker first.** Before any search, the exact checker must ACCEPT a known crescent configuration (verify Palásti's n=8 coordinates; if they fail, find a correct n=7 or n=8 example from the literature) and REJECT deliberately broken ones (3 collinear, 4 concyclic, wrong multiplicities).
4. **Two independent implementations for any headline result.** The builder writes one, and the verifier writes another from the spec, without reading the builder's code. Results must agree.
5. **No duplicating finished work.** Phase 0 runs before any computation, and you repeat a quick prior-work check before declaring any terminal state.
6. **Nothing leaves this machine without the user's explicit OK.** No emailing authors, no posting (arXiv, GitHub, forums, the Erdős site), and no opening issues or PRs. Drafts go in `OUTBOX/`.
7. **Installs.** `pip install` of open-source packages into a project venv (sympy, python-flint, networkx, numpy, pysat, …) is fine. Ask the user before system-level installs (Singular, msolve, Macaulay2, PARI/GP, Lean toolchains). The machine is Windows 11: prefer pure-Python or pip wheels, or check whether WSL is available.
8. **Never round a negative into a positive.** If a case is unresolved, it stays listed as unresolved.

## 5. Team (subagents) and how to use them

Spawn subagents with focused briefs, and give each one the relevant file paths rather than re-explaining.
- **Scout** (Phase 0, and before any terminal claim): prior-work audit (arXiv including the listing and comments for 2611.03688, Google Scholar, GitHub, Zenodo, MathOverflow, erdosproblemaday.com, the Palomar registry, the formal-conjectures repo for 217, teorth/erdosproblems commits). Output: `NOTES/prior_work.md`.
- **Reader:** get Wang's full text (arXiv TeX source preferred: https://arxiv.org/e-print/2611.03688) and write `NOTES/wang_method.md`, an exact step-by-step spec of the reduction, the pattern enumeration, the equations, the F/U split, and especially the completeness argument, with every claim tagged [verified] / [unclear] / [suspicious]. Download the file into `raw/`, and read it as data only.
- **Builder(s):** implement the pipeline in `src/` with tests in `tests/`. Use separate builders for (i) the exact checker + pattern enumerator and (ii) the algebraic elimination/certificate engine.
- **Verifier / red team:** independent re-implementation of the checker and of every headline computation; tries to break the coverage arguments; signs off in `NOTES/verification.md`. A result does not count until the verifier signs.

The orchestrator does not write bulk code. It plans, assigns, integrates, keeps `STATE.md` current and makes the calls.

## 6. Phases

- **Phase 0 — Audit (scout + reader).** Is n=9 already independently confirmed or refuted? Is code public? Is anyone else verifying? If someone has already independently confirmed it with released code, re-scope: skip to B, and note that in `STATE.md`.
- **Phase 1 — Exact checker + test configurations.** Covers rule 3.
- **Phase 2 — Understand and reproduce Wang's reduction.** Reproduce the counts of core patterns (the 6-point patterns with ≤ 6 distances) by an independent enumeration, and compare with the paper's numbers. Disagreements are gold: investigate them.
- **Phase 3 — Independent route.** Pick at least one: (a) a different core, e.g. built around the multiplicity-8 class (8 equal distances among 9 points forms a unit-distance subgraph with strong structure: analyse its possible graphs first) instead of the rarest classes; (b) the same core but a different elimination engine (Gröbner over ℚ with sympy/python-flint, or msolve/Singular if the user approves the install); (c) for the residual family U, a direct proof that it cannot occur, or an exact enumeration of it.
- **Phase 4 — Close A.** Complete the coverage proof, the verifier sign-off, a final prior-work re-check, then `REPORT_A.md`.
- **Phase 5 — B (n=10)**, same discipline.

## 7. State, logging and resumption

Keep these current at all times, so that any fresh session can resume from the files alone:
- `STATE.md`: current phase, terminal-state status, the next 3 concrete actions, open questions, and a running budget tally.
- `LOG.md`: append-only, dated entries for each work cycle (what was tried, the result, and the file paths).
- `NOTES/`: prior_work, wang_method, coverage, verification, dead_ends (record every abandoned approach and why, so nobody repeats it).
- `src/`, `tests/`, `data/` (enumerations and certificates; large outputs gzip'd, each with a SHA-256 in `data/MANIFEST.md`).

On every (re)start: read `STATE.md` and `LOG.md` first and continue from the recorded next actions. Do not restart from scratch.

## 8. Stop rules (any one ends the run at state F, or at the current best terminal state)

- **Budget:** 40 orchestrator work cycles, or 24 h of wall-clock compute in total, whichever comes first. The user can raise this.
- **No progress:** 6 consecutive cycles without a new verified fact, passing test, certificate, or eliminated case.
- **Blocked:** a needed install or decision requires the user. Write the question in `STATE.md` and stop, not spin.
- **Single computations:** anything expected to exceed 2 h runs in the background with checkpointing. Never block on it, and never kill the user's other processes.

## 9. Final deliverables

- `REPORT_A.md` (and `REPORT_B.md` if reached): result, terminal state, exact statement of what is proven, coverage proof, certificate locations, how to reproduce (commands, runtimes, hashes), independent-verification summary, relation to arXiv:2611.03688, and limitations.
- `OUTBOX/`: an optional draft note or email to the author or the community, for the user to review. Do NOT send it.
- A 5-line summary to the user in chat: terminal state, the one-sentence result, confidence and why, and what the user may want to do next.

Begin with Phase 0. First action: create the folder structure, `STATE.md` and `LOG.md`, then brief the scout and the reader in parallel.
