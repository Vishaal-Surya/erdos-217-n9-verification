"""Planted known-answer tests for the n=9 coverage search (Builder C).

For each of the 7 arrangements, build a random 9-point labelled pattern P with the arrangement's low
pairs (delta1 -> class size 1, delta2 -> size 2) and high class sizes exactly 3..8, hill-climbed (label
swaps keep the class sizes) until it passes Lemma 2.1 (i)-(iii) and Lemma 2.2 (src/canon.local_ok).
P is NOT claimed realizable. Set U := canonical forms of P's all-high 6-subsets, optionally plus
`noise` random members of data/patterns6_affine.txt.gz. The search must report P (up to equivalence).

Writes data/coverage_planted_U_<arr>.txt and data/coverage_planted_targets_<arr>.txt.
Usage: python planted.py [seed=1] [noise=0]
"""
import gzip
import itertools
import os
import random
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "src", "coverage"))
import canon  # noqa: E402
import arrangements  # noqa: E402

DATA = os.path.join(ROOT, "data")
N = 9
PAIRS = list(itertools.combinations(range(N), 2))


def violations(lab):
    c = canon.matrix(lab, N)
    v = 0
    for a in range(N):
        vals = [c[a][x] for x in range(N) if x != a]
        v += sum(max(0, vals.count(t) - 3) for t in set(vals))
    for a, b in itertools.combinations(range(N), 2):
        keys = [(c[a][x], c[b][x]) for x in range(N) if x not in (a, b)]
        v += sum(max(0, keys.count(q) - 2) for q in set(keys))
        v += max(0, sum(1 for p, q in keys if p == q) - 2)
    for a, b, cc, d in itertools.combinations(range(N), 4):
        if c[a][b] == c[cc][d] and c[a][cc] == c[b][d] and c[a][d] == c[b][cc]:
            v += 1
    return v


def plant(rep, rng):
    low = {tuple(sorted((x - 1, y - 1))): (1 if i == 0 else 2) for i, (x, y) in enumerate(rep)}
    high = [p for p in PAIRS if p not in low]
    for attempt in range(200):
        sizes = [s for s in range(3, 9) for _ in range(s)]
        rng.shuffle(sizes)
        lab = {p: low[p] for p in low}
        lab.update(dict(zip(high, sizes)))
        seq = [lab[p] for p in PAIRS]
        cur = violations(seq)
        for it in range(4000):
            if cur == 0:
                return seq
            i, j = rng.sample([PAIRS.index(p) for p in high], 2)
            if seq[i] == seq[j]:
                continue
            seq[i], seq[j] = seq[j], seq[i]
            nv = violations(seq)
            if nv <= cur or rng.random() < 0.02:
                cur = nv
            else:
                seq[i], seq[j] = seq[j], seq[i]
    raise RuntimeError("could not plant")


def cores(seq):
    lab = dict(zip(PAIRS, seq))
    out = set()
    for T in itertools.combinations(range(N), 6):
        sub = [lab[(a, b)] for a, b in itertools.combinations(T, 2)]
        if any(x in (1, 2) for x in sub):
            continue
        out.add(canon.canonical(sub, 6))
    return out


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    noise = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    rng = random.Random(seed)
    extra = []
    if noise:
        with gzip.open(os.path.join(DATA, "patterns6_affine.txt.gz"), "rt") as f:
            allp = [l.strip() for l in f if l.strip()]
        extra = rng.sample(allp, noise)
    for name, rep in arrangements.OURS.items():
        seq = plant(rep, rng)
        assert canon.local_ok(seq, N)
        s = "".join(map(str, seq))
        cs = cores(seq)
        tag = f"s{seed}_n{noise}_{name}"
        with open(os.path.join(DATA, f"coverage_planted_U_{tag}.txt"), "w") as f:
            f.write("\n".join(sorted(cs | set(extra))) + "\n")
        with open(os.path.join(DATA, f"coverage_planted_targets_{tag}.txt"), "w") as f:
            f.write(f"planted_{name} {s}\n")
        print(name, s, "cores:", len(cs))


if __name__ == "__main__":
    main()
