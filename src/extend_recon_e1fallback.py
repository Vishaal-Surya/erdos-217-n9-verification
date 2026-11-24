"""Fallback reconstruction via Wang's Sec. 6 multiplication-matrix route (algebra.engine1 +
algebra.f_univariate, read-only) for F patterns where the E2 Groebner basis hit its limits.
Lemma 6.1: every normalized GP realization is z_j = g_j(lam) at a root of f, so taking ALL
irreducible factors of f and reconstructing each with extend_recon.reconstruct_core is complete.
python src/extend_recon_e1fallback.py PATTERN... >> data/F_recon_e1fallback.jsonl
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import algebra as A
import extend_recon as R
from extend_kfield import KField, poly_to_list
import flint


def recon_e1(s):
    t0 = time.perf_counter()
    e = A.engine1(s)
    out = {"pattern": s, "r": e["r"], "e2": "UNRESOLVED", "route": "E1", "e1": e["verdict"], "factors": []}
    if e["verdict"] != "F":
        out["status"] = "E1 not F"; return out
    U = A.f_univariate(e)
    if U.get("failed"):
        out["status"] = "UNRESOLVED (M-choice failed)"; return out
    out.update(a=e["a"], h=U["h"], c=U["c"], commute=U["commute"], f=str(U["f"]))
    c, r = A.parse(s)
    a = e["a"] - 1
    for q, mult in U["factors"][1]:
        K = KField(flint.fmpq_poly(q))
        z = [K.one if j == a else K.red(U["g"][j]) for j in range(r)]
        rec = R.reconstruct_core(c, r, K, z)
        fr = {"q": poly_to_list(K.q), "deg": K.deg, "status": rec["status"], "nreal": rec.get("nreal", K.nreal()),
              "z": [poly_to_list(x) for x in z]}
        if rec["status"] == "REALIZED":
            fr["pd_roots"] = rec["pd_roots"]
            fr["roots"] = [K.root_str(k) for k in rec["pd_roots"]]
            fr["G"] = [poly_to_list(x) for x in rec["G"]]
            fr["points"] = [[poly_to_list(p[0]), poly_to_list(p[1])] for p in rec["points"]]
        out["factors"].append(fr)
    out["t"] = time.perf_counter() - t0
    return out


if __name__ == "__main__":
    for s in sys.argv[1:]:
        print(json.dumps(recon_e1(s)), flush=True)
