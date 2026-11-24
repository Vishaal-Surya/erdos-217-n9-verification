"""Builder C: tests for the coverage search (src/coverage/).

Run:  .venv/Scripts/python -m pytest tests/test_coverage.py -q      (about 3-4 minutes, 1 worker)

Long runs (n=7 with U = all local 4-point patterns, planted n=9 for all 7 arrangements, stand-in U
subsets) are recorded in NOTES/coverage.md with their logs in data/coverage_*.log.
"""
import os
import re
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "src", "coverage"))
import arrangements  # noqa: E402
import canon  # noqa: E402

JDK = os.path.join(os.environ["JAVA_HOME"], "bin") if os.environ.get("JAVA_HOME") else os.path.dirname(shutil.which("java") or "java")
BUILD = os.path.join(ROOT, "src", "coverage", "build")
DATA = os.path.join(ROOT, "data")


@pytest.fixture(scope="session", autouse=True)
def compiled():
    os.makedirs(BUILD, exist_ok=True)
    subprocess.check_call([os.path.join(JDK, "javac"), "-d", BUILD, os.path.join(ROOT, "src", "coverage", "CoverageSearch.java")])
    if not os.path.exists(os.path.join(DATA, "coverage_targets_n8.txt")):
        subprocess.check_call([sys.executable, os.path.join(ROOT, "src", "coverage", "targets.py")])


def run(*args):
    out = subprocess.run([os.path.join(JDK, "java"), "-cp", BUILD, "CoverageSearch", "--threads", "1", *args],
                         cwd=ROOT, capture_output=True, text=True, check=True).stdout
    return out


def hits(out):
    m = re.search(r"TARGETS HIT: (\d+) \[(.*)\]", out)
    return set(x.strip() for x in m.group(2).split(",") if x.strip())


def total(out):
    return int(re.search(r"TOTAL leaves=(\d+)", out).group(1))


@pytest.mark.parametrize("n", [7, 8, 9])
def test_arrangements_are_the_seven_orbits(n):
    orb = arrangements.orbits(n)
    assert len(orb) == 7
    from math import comb
    p = comb(n, 2)
    assert sum(orb.values()) == p * comb(p - 1, 2)
    seen = set()
    for rep in arrangements.OURS.values():
        r0 = tuple(tuple(x - 1 - (9 - n) for x in e) for e in rep)
        c = arrangements.canon_of(r0)
        assert c in orb and c not in seen
        seen.add(c)
        # positions 1..n-3 span no low pair (the seed block)
        assert all(max(e) >= n - 3 for e in r0)


def test_n8_palasti_found_from_its_cores():
    out = run("--n", "8", "--U", "data/coverage_coresU_n8_all.txt", "--targets", "data/coverage_targets_n8.txt")
    assert hits(out) == {"palasti8", "tri8_0"}


def test_n8_negative_control_drop_one_core():
    # U minus one of Palasti's 5-point cores: Palasti must disappear (the U test really prunes).
    cores = open(os.path.join(DATA, "coverage_coresU_n8_palasti8.txt")).read().split()
    tmp = os.path.join(DATA, "coverage_tmp_U_drop1.txt")
    with open(tmp, "w") as f:
        f.write("\n".join(cores[1:]) + "\n")
    out = run("--n", "8", "--U", tmp, "--targets", "data/coverage_targets_n8.txt")
    os.remove(tmp)
    assert "palasti8" not in hits(out)


def test_n7_lattice_examples_found_from_cores():
    out = run("--n", "7", "--U", "data/coverage_coresU_n7_all.txt", "--targets", "data/coverage_targets_n7.txt",
              "--arr", "star,chain_d1_end,P3K2_d1_meeting")
    assert hits(out) == {f"tri7_{i}" for i in range(8)}


def test_lookahead_and_mode_do_not_change_results():
    base = ["--n", "8", "--U", "data/coverage_coresU_n8_all.txt", "--arr", "P3K2_d1_meeting,chain_d1_end"]
    a = run(*base)
    b = run(*base, "--nolook")
    c = run(*base, "--mode", "canon")
    strip = lambda o: [re.sub(r"nodes=\d+ |time=[\d.]+s", "", l) for l in o.splitlines() if l.startswith(("arr=", "TOTAL"))]
    assert strip(a) == strip(b) == strip(c)


def test_paranoid_planted_disjoint_n9():
    tag = "s1_n0_disjoint"
    if not os.path.exists(os.path.join(DATA, f"coverage_planted_U_{tag}.txt")):
        subprocess.check_call([sys.executable, os.path.join(ROOT, "src", "coverage", "planted.py"), "1", "0"])
    out = run("--n", "9", "--arr", "disjoint", "--U", f"data/coverage_planted_U_{tag}.txt",
              "--targets", f"data/coverage_planted_targets_{tag}.txt", "--paranoid")
    assert hits(out) == {"planted_disjoint"}


def test_targets_are_crescent_and_canonical_agree():
    # Java canonComplete vs python brute canonical (column 4 of the targets file) are compared indirectly:
    # the Java side matched the targets, here we check the python canonical is a fixed point.
    for n in (7, 8):
        for line in open(os.path.join(DATA, f"coverage_targets_n{n}.txt")):
            if line.startswith("#"):
                continue
            name, s, arr, cf = line.split()
            assert sorted(s.count(str(i)) for i in range(1, n)) == list(range(1, n))
            assert min("".join(cf[p] for p in pm) for pm in canon._perm_maps(n)) == cf
