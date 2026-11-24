#!/usr/bin/env bash
# Builder C: real coverage run on the E1 residual family U.  CPU cap: 2 threads total.
# Usage (from attack-217/):  bash src/coverage/run_real.sh [data/E1_U.txt.gz]
set -euo pipefail
U="${1:-data/E1_U.txt.gz}"
JDK="${JAVA_HOME:-}"; JDK="${JDK:+$JDK/bin}"; JDK="${JDK:-$(dirname "$(command -v java)")}"
J=("$JDK/java" -Xmx4g -cp src/coverage/build CoverageSearch --n 9 --U "$U" --threads 2)
"$JDK/javac" -d src/coverage/build src/coverage/CoverageSearch.java
sha256sum "$U" | tee data/coverage_E1_run.log
# 1) main search, all 7 arrangements (orbit-set membership, lookahead on)          ~10 min
"${J[@]}" --out data/coverage_E1_found.txt 2>&1 | grep -v '^  \[' | tee -a data/coverage_E1_run.log
# 1b) affine post-filter N6 on the surviving complete patterns (flint + Fraction)   ~5 s
.venv/Scripts/python src/coverage/affine_n.py data/coverage_E1_found.txt data/coverage_E1_affine.txt --check | tee -a data/coverage_E1_run.log
# 2) cross-check: canonical-form membership code path, all arrangements             slower
"${J[@]}" --mode canon --out data/coverage_E1_found_canonmode.txt 2>&1 | grep -v '^  \[' > data/coverage_E1_run_canonmode.log
# 3) comparable of Wang's "13 seven-point patterns" (no affine tests here)
"${J[@]}" --level7 data/coverage_E1_level7.txt 2>&1 | grep -v '^  \[' > data/coverage_E1_level7.log
# 4) paranoid spot-check (incremental==full local at every node, leaf re-check), no lookahead
"${J[@]}" --arr triangle,chain_d1_middle --nolook --paranoid 2>&1 | grep -v '^  \[' > data/coverage_E1_run_paranoid_spot.log
