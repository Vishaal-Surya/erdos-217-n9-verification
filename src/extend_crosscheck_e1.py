"""Cross-check of the E2 reconstruction against Wang's multiplication-matrix route (Engine 1).

For each F pattern: algebra.engine1 -> algebra.f_univariate (Wang Sec. 6: M_j, f, g_j), factor f,
apply Wang's retention filter (real root; z_j != 0 and z_j != z_k in K_q) -> count comparable to
Wang's "4119 retained irreducible factors"; then reconstruct each retained factor with
extend_recon.reconstruct_core and compare the realized cores with the E2 ones (data/F_recon*.jsonl)
by a numeric fingerprint (label values normalised by label 1, evaluated at every PD real root,
rounded to 25 significant digits).  The fingerprint is only a cross-check (the primary route is E2).

python src/extend_crosscheck_e1.py data/F_recon_final.jsonl data/F_recon_e1check.jsonl --workers 2
"""
import json, os, sys, time, multiprocessing as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
MAXW = 2


def fingerprints(K, z, roots):
    out = set()
    for k in roots:
        v0 = K.eval_float  # noqa
        import flint
        rts = K.real_roots(256)
        old = flint.ctx.prec
        try:
            flint.ctx.prec = K._roots_prec
            vals = []
            for x in z:
                v = flint.arb(0)
                for cf in reversed((x % K.q).coeffs()):
                    v = v * rts[k] + flint.arb(cf)
                vals.append(v)
            n0 = vals[0]
            out.add(tuple(float((v / n0).mid()) for v in vals))
        finally:
            flint.ctx.prec = old
    return out


def _round(fp):
    return tuple(float(f"{x:.12g}") for x in fp)


def work(item):
    s, e2fps = item
    import algebra as A
    import extend_recon as R
    from extend_kfield import KField
    import flint
    t0 = time.perf_counter()
    out = {"pattern": s}
    try:
        e = A.engine1(s)
        if e["verdict"] != "F":
            out["status"] = "E1 not F?"; return out
        U = A.f_univariate(e)
        if U.get("failed"):
            out["status"] = "M-choice failed (UNRESOLVED in Wang's route)"; return out
        c, r = A.parse(s)
        a = e["a"] - 1
        nf = 0; nreal = 0; nret = 0; nreal_core = 0
        fps = set()
        for q, mult in U["factors"][1]:
            nf += 1
            K = KField(flint.fmpq_poly(q))
            z = [K.one if j == a else K.red(U["g"][j]) for j in range(r)]
            if K.nreal() == 0:
                continue
            nreal += 1
            if any(x == 0 for x in z) or len({K.key(x) for x in z}) < r:
                continue
            nret += 1
            rec = R.reconstruct_core(c, r, K, z)
            if rec["status"] == "REALIZED":
                nreal_core += 1
                fps |= {_round(f) for f in fingerprints(K, z, rec["pd_roots"])}
        e2 = {_round(tuple(x)) for x in e2fps}
        out.update(status="ok", h=U["h"], c=U["c"], factors=nf, with_real_root=nreal,
                   retained=nret, realized=nreal_core, n_fp_e1=len(fps), n_fp_e2=len(e2),
                   match=(fps == e2), commute=U.get("commute"))
    except Exception as ex:
        out["status"] = "ERROR " + repr(ex)
    out["t"] = time.perf_counter() - t0
    return out


def e2_fingerprints(d):
    from extend_kfield import KField, poly_from_list
    fps = []
    for fr in d["factors"]:
        if fr["status"] != "REALIZED":
            continue
        K = KField(poly_from_list(fr["q"]))
        z = [K.red(poly_from_list(x)) for x in fr["z"]]
        fps += [list(f) for f in fingerprints(K, z, fr["pd_roots"])]
    return fps


def main():
    inp, outp = sys.argv[1], sys.argv[2]
    w = 1
    if "--workers" in sys.argv:
        w = int(sys.argv[sys.argv.index("--workers") + 1])
    w = max(1, min(w, MAXW))
    items = []
    for l in open(inp):
        d = json.loads(l)
        items.append((d["pattern"], e2_fingerprints(d)))
    done = set()
    if os.path.exists(outp):
        for l in open(outp):
            done.add(json.loads(l)["pattern"])
    todo = [it for it in items if it[0] not in done]
    print(f"{len(items)} patterns, {len(todo)} to do, workers={w}", flush=True)
    t0 = time.time()
    with open(outp, "a") as f, mp.get_context("spawn").Pool(w, maxtasksperchild=100) as pool:
        for k, res in enumerate(pool.imap_unordered(work, todo, chunksize=1)):
            f.write(json.dumps(res) + "\n"); f.flush()
            if not res.get("match", False):
                print("NOTE", json.dumps(res), flush=True)
            if k % 50 == 0:
                print(f"{k+1}/{len(todo)} {time.time()-t0:.0f}s", flush=True)
    print("DONE", time.time() - t0, flush=True)


if __name__ == "__main__":
    main()
