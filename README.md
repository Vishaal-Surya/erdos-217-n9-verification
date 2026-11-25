# Erdős Problem #217, n = 9: independent verification that no 9-point crescent configuration exists

A **crescent configuration** is a set of n points in the plane that meets three conditions:
- no 3 points are collinear and no 4 are concyclic;
- the points determine exactly n−1 distinct distances;
- for each i = 1,…,n−1, some distance occurs exactly i times.

Erdős asked for which n such configurations exist ([Problem #217](https://www.erdosproblems.com/217)). Planar examples are known for n ≤ 8; Palásti gave one for n = 8.

## What this repository claims

**We confirm independently, by exact computation, that no planar crescent configuration exists for n = 9.** This reproduces, stage by stage, the computer-assisted proof of Y. Wang, *The nonexistence of nine-point crescent configurations*, arXiv:2611.03688 (27 Oct 2026). That paper released no code or data.

- The **reduction is Wang's**. Independence here means two things: every stage was implemented twice, independently, from the written specification only (the second implementation never saw the first one's code), and the algebraic elimination was run with a second, different engine (Gröbner bases with saturation over Q, in addition to Wang's linear-algebra method).
- **Every count reported in the paper was reproduced exactly:**
  - 3 / 22 / 932 / 625,472 patterns;
  - 344,835 after the affine tests;
  - 4,637 = 4,450 (F) + 187 (U), with U = 3/9/175 by label count;
  - 13 seven-point patterns;
  - 4,119 real-root factors;
  - maximum matrix order 12;
  - 0 nine-point completions.
- **Every one of the 340,198 algebraic exclusions has two independent exact certificates**: a replayed linear-algebra certificate over Q, and a Gröbner computation that returns the unit ideal (or puts a product of nonzero factors in the ideal).
- **What is not claimed:**
  - anything for n ≥ 10;
  - a human-readable or formally verified (Lean) proof;
  - that C+Q = 1 is impossible (see the side finding below).

**Provenance and AI disclosure.** This work was carried out by Vishaal Surya using AI assistance: Anthropic's Claude models, run as an orchestrator with separate builder and verifier agents. All code, computations and write-ups were produced that way.
- The "independent" implementations are separate agents. Each was given only the written specification, and the verifiers were never shown the builders' source.
- Trust should rest on the released code, the exact certificates and the hashes, which anyone can re-run. It should not rest on who or what wrote them.
- Independent human reproduction is very welcome. Please open an issue with your results.

The full write-up is in [`REPORT_A.md`](REPORT_A.md), including the caveats (§5). The main one is that correctness rests on Wang's reduction, which we re-derived but did not replace, together with the correctness of the released code.

## Side finding: a (C, Q) = (0, 2) near-miss

This configuration has 9 points, distance multiplicities exactly 1..8 and **no three points collinear**, but it has **two** concyclic quadruples. Its squared distances lie in Q(√21). Scaled so that the multiplicity-8 squared distance is 1, they are:

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

It improves the (C, Q) = (1, 2) lattice example posted in the erdosproblems.com forum thread. Three independent exact codes confirm it; see [`NOTES/near_misses.md`](NOTES/near_misses.md) and `verify/orch/nearmiss_config2.py`.

## Layout

| Path | Contents |
|---|---|
| `src/checker.py` | exact crescent checker (number fields via python-flint) |
| `src/enum/PatternEnum.java`, `src/canon.py`, `src/affine.py` | 6-point pattern enumeration, canonical forms, exact affine tests |
| `src/algebra*.py` | Engine 1 (Wang's quartic-span method, with certificates) and Engine 2 (Gröbner/saturation); batch drivers |
| `src/extend_*.py` | exact reconstruction of finite cores (case F) and extension to 9 points |
| `src/coverage/` | combinatorial coverage search for case U (Java) and the 9-point affine post-filter |
| `tests/` | pytest suites: checker, enumeration soundness, algebra soundness, extension, coverage |
| `verify/` | the verifiers' independent implementations, written without reading `src/`: `checker_v.py`, `va/` (enumeration, certificate replay, classifier), `vb/` (coverage, rigorous extension), `orch/` (near-miss checks) |
| `data/` | pattern lists, F/U lists, certificates, verdicts, run logs; hashes in `data/MANIFEST.md` and `SHA256SUMS` |
| `NOTES/` | method spec of Wang's paper, proofs (coverage, extension), verification sign-off, prior work, near misses |
| `LOG.md`, `STATE.md`, `PROMPT.md` | research log and brief (the project was run as an orchestrated multi-agent effort) |

Files over 1 MB are stored gzipped (`*.gz`). Scripts that read the original uncompressed file name may need `gunzip -k` first. Wang's TeX source is not redistributed; get it from arXiv.

## Reproducing

Requirements:
- Python ≥ 3.12 with `pip install sympy python-flint numpy networkx pytest`;
- Java 21 (`JAVA_HOME` set, or `java` and `javac` on PATH).

Our runs used a 20-core machine capped at ≤ 8 worker processes.

```bash
python -m pytest tests/test_checker.py -q                                   # ~2 s
javac -d src/enum/build src/enum/PatternEnum.java
java -Xmx8g -cp src/enum/build PatternEnum A data 4 6                       # 6-point patterns, seconds
python src/affine.py data/patterns6_local.txt.gz data/patterns6_affine.txt.gz 4 flint
python src/algebra_full.py --workers 4 --chunk 5000                         # Engine 1, ~3 h
bash src/run_e2_all.sh                                                      # Engine 2, ~2.5 h (+ retry, see REPORT_A §7)
bash src/coverage/run_real.sh data/E1_U.txt.gz                              # case U, ~10 min
python src/coverage/affine_n.py data/coverage_E1_found.txt data/coverage_E1_affine.txt --check
# case F: see NOTES/extension.md (src/extend_run_recon.py, src/extend_run.py), ~10 min
# independent verifiers: verify/va/run_all.py (~3.6 h), verify/vb/*.py
```

Runtimes and per-stage commands are in `REPORT_A.md` §7 and in the individual NOTES. Check data integrity with `sha256sum -c SHA256SUMS`.

## License

MIT; see [`LICENSE`](LICENSE). The cited papers belong to their authors.
