# STATE: Erdős #217 attack (n=9 verification)

**Current phase:** TERMINAL. Goal A was reached in state **A1 (independent confirmation)** on 2026-11-25. Goal B (n=10) is NOT started in this session: the user changed the scope, and B runs in a separate session.

**Terminal state: A1.** No planar 9-point crescent configuration exists. Every stage of arXiv:2611.03688 was reproduced exactly with two independent implementations, every exclusion has two independent exact certificates (E1 replay by V-a, E2 Gröbner), and the verifier signed off (NOTES/verification.md).

## Deliverables
- REPORT_A.md
- NOTES/near_misses.md: (C,Q) = (0,2) over Q(√21)
- OUTBOX/draft_note_to_author.md: NOT sent
- release/: publishable package, a local git repo with one commit, no remote, not pushed. The parent session confirms the push with the user.

## Next 3 actions
1. (User) Decide whether to push release/ to GitHub; the parent session handles this.
2. (User) Decide whether to send the OUTBOX draft to the author.
3. (Separate session) Goal B, n=10, building on this pipeline.

## Open questions
- None blocking. Caveats are listed in REPORT_A.md §5: same reduction as Wang; Cor 4.5 follows from sound certificates plus a complete enumeration rather than being proved separately; 6 F patterns were reconstructed by a single engine.

## Environment
- Windows 11, Python 3.14.4, project venv `.venv/` (sympy, python-flint, numpy, networkx, pytest), Java 21.
- **STANDING RULE: max 8 concurrent worker processes total (all jobs and subagents combined). This was the user's request because of a CPU thermal limit (2026-11-24). Scripts default to ≤4 workers and never os.cpu_count().** No background jobs running.

## Budget tally (final)
- Orchestrator cycles: 21 / 40
- Compute wall-clock: about 13 h of background job wall time (E1 ~3 h, E2 ~2.5 h + retry ~2 h, V-a ~3.6 h, others < 1 h each) / 24 h
- Consecutive no-progress cycles: 0 / 6
