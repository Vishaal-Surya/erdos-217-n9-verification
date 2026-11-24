"""Cross-check: builder src/checker.py vs verifier verify/checker_v.py on every saved
example (both sources) and on a perturbed variant of each (one lattice point moved)."""
import json, re, sys, itertools
sys.path.insert(0, "src"); sys.path.insert(0, "verify")
import checker as B
import checker_v as V

def spec_from_lattice(ab):
    pts = []
    for a, b in ab:
        from fractions import Fraction as Fr
        x = Fr(a) + Fr(b, 2); y = Fr(b, 2)
        pts.append([str(x), f"({y})*sqrt(3)"])
    return {"field": {"type": "sqrt", "d": 3}, "points": pts}

def both(ab):
    vb = B.check_spec(spec_from_lattice(ab)).ok
    vv = V.check([V.lattice_to_q3(a, b) for a, b in ab])[0]
    return vb, vv

sets = []
ke = json.load(open("data/known_examples.json"))
for k, r in ke["results"].items():
    for ex in r.get("examples", []):
        sets.append(("builder", int(k), [tuple(p) for p in ex["lattice_ab"]]))
for line in open("verify/lattice_n7_r3.txt"):
    m = re.search(r"lattice\(i,j\) = (\[.*?\])", line)
    if m:
        sets.append(("verifier", 7, [tuple(p) for p in eval(m.group(1))]))
dis = 0; acc = 0; rej = 0
for src, k, ab in sets:
    vb, vv = both(ab)
    if not (vb and vv): print("NOT ACCEPTED BY BOTH", src, k, ab, vb, vv); dis += 1
    else: acc += 1
    for idx, (da, db) in itertools.product(range(len(ab)), [(1, 0), (0, 1), (2, -1), (1, 1)]):
        ab2 = list(ab); ab2[idx] = (ab[idx][0] + da, ab[idx][1] + db)
        if len(set(ab2)) < len(ab2): continue
        vb2, vv2 = both(ab2)
        if vb2 != vv2: print("DISAGREE", src, k, ab2, vb2, vv2); dis += 1
        else: rej += (not vb2)
print(f"examples accepted by both: {acc}/{len(sets)}; perturbed agreements incl. {rej} joint rejects; disagreements: {dis}")
