"""Known-answer soundness test for the E2 reconstruction (src/extend_recon.py).

data/algebra_realizable.json holds exact general-position 6-point sets (coordinates in Q(sqrt3))
whose labelled pattern is computed from the points.  For every such set whose pattern is in F
(E1 verdict F, data/algebra_realizable_results.jsonl), the actual normalized label tuple must be
found EXACTLY: some reconstructed factor q has q(lam_actual) = 0 and z_l(lam_actual) = actual z_l
for every label l, the factor must be REALIZED, and lam_actual must lie in a PD root's ball.
Run: .venv/Scripts/python -m pytest tests/test_extend_recon.py -q   (or python tests/... for all)
"""
import json, os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
import extend_recon as R
import flint


class Q3:
    def __init__(s, a, b=0): s.a, s.b = Fr(a), Fr(b)
    def __add__(s, o): o = o if isinstance(o, Q3) else Q3(o); return Q3(s.a + o.a, s.b + o.b)
    __radd__ = __add__
    def __sub__(s, o): o = o if isinstance(o, Q3) else Q3(o); return Q3(s.a - o.a, s.b - o.b)
    def __mul__(s, o): o = o if isinstance(o, Q3) else Q3(o); return Q3(s.a * o.a + 3 * s.b * o.b, s.a * o.b + s.b * o.a)
    __rmul__ = __mul__
    def inv(s): n = s.a * s.a - 3 * s.b * s.b; return Q3(s.a / n, -s.b / n)
    def __eq__(s, o): o = o if isinstance(o, Q3) else Q3(o); return s.a == o.a and s.b == o.b
    def __float__(s): return float(s.a) + float(s.b) * 3 ** 0.5


def evalpoly(coeffs, x):
    v = Q3(0)
    for c in reversed(coeffs):
        v = v * x + Q3(Fr(c))
    return v


def load_cases(limit=None):
    F = set()
    for l in open(os.path.join(ROOT, "data/algebra_realizable_results.jsonl")):
        d = json.loads(l)
        if d.get("e1", {}).get("verdict") == "F":
            F.add(d["pattern"])
    J = json.load(open(os.path.join(ROOT, "data/algebra_realizable.json")))
    cases = []
    for gen, lst in J.items():
        for pat, pts in lst:
            if pat in F:
                P = [(Q3(Fr(x[0]), Fr(x[1])), Q3(Fr(y[0]), Fr(y[1]))) for x, y in pts]
                cases.append((gen, pat, P))
    return cases[:limit] if limit else cases


def check_case(pat, P):
    c, r = R.A.parse(pat)
    z = [None] * r
    for i in range(6):
        for j in range(i + 1, 6):
            dx, dy = P[i][0] - P[j][0], P[i][1] - P[j][1]
            v = dx * dx + dy * dy
            if z[c[i][j]] is None:
                z[c[i][j]] = v
            assert z[c[i][j]] == v, "pattern/points mismatch"
    o = R.reconstruct_pattern(pat)
    assert o["e2"] == "ZERODIM", (pat, o["e2"])
    a = o["a"] - 1
    zn = [x * z[a].inv() for x in z]
    lam = Q3(0)
    al = o.get("alpha", [])
    for k, l in enumerate(o["others"]):
        lam = lam + zn[l - 1] * al[k]
    hits = []
    for f in o["factors"]:
        if any(c_ != "0" for c_ in f["q"]) and evalpoly(f["q"], lam) == 0:
            if all(evalpoly(f["z"][l], lam) == zn[l] for l in range(r)):
                hits.append(f)
    assert len(hits) == 1, (pat, len(hits))
    f = hits[0]
    assert f["status"] == "REALIZED", (pat, f["status"])
    # lam_actual inside one of the PD root balls
    K = R.KField(R.flint.fmpq_poly([R.fmpq(*map(int, s.split('/'))) if '/' in s else int(s) for s in f["q"]]))
    rts = K.real_roots(128)
    lf = float(lam)
    inside = [k for k in f["pd_roots"] if abs(float(rts[k].mid()) - lf) < 1e-9 * max(1, abs(lf))]
    assert inside, (pat, "actual root not a PD root")
    return True


def test_known_F_realizations():
    cases = load_cases()
    assert len(cases) >= 500
    for gen, pat, P in cases:
        assert check_case(pat, P)


if __name__ == "__main__":
    import time
    t = time.time()
    cases = load_cases()
    n = 0
    for gen, pat, P in cases:
        check_case(pat, P); n += 1
    print(f"OK: {n} known exact F realizations recovered exactly ({time.time()-t:.1f}s)")
