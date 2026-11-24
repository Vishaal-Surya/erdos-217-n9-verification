"""Driver: extend every realizable F core (from data/F_recon*.jsonl) to 9 points (exact, K_q).

python src/extend_run.py data/F_recon_final.jsonl data/F_extend.jsonl --workers 2
One JSON line per realizable (pattern, factor): stats of every stage and any 9-point set found
(coordinates as polynomial coefficient lists in u modulo q).  Resumable (skips done keys).
"""
import json, os, sys, time, multiprocessing as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
MAXW = 2  # user CPU cap for builder D


def load_core(fr):
    from extend_kfield import KField, poly_from_list
    K = KField(poly_from_list(fr["q"]))
    G = [K.red(poly_from_list(x)) for x in fr["G"]]
    S = [(K.red(poly_from_list(p[0])), K.red(poly_from_list(p[1]))) for p in fr["points"]]
    return K, G, S


def work(item):
    key, pattern, fr, n = item[:4]
    full = len(item) > 4 and item[4]
    from extend_core import Extender
    from extend_kfield import poly_to_list
    t0 = time.perf_counter()
    try:
        K, G, S = load_core(fr)
        X = Extender(K, G, S, n=n, full=full)
        found = X.run()
        out = {"key": key, "pattern": pattern, "deg": K.deg, "r": X.r, "m": X.m,
               "stats": dict(X.stats), "found": len(found)}
        if X.near:
            out["near_non_gp"] = [[[poly_to_list(p[0]), poly_to_list(p[1])] for p in P] for P in X.near]
            out["q"] = fr["q"]; out["G"] = fr["G"]
        if found:
            out["configs"] = [{"core_is_high": f["core_is_high"],
                               "points": [[poly_to_list(p[0]), poly_to_list(p[1])] for p in f["points"]]}
                              for f in found]
            out["q"] = fr["q"]; out["G"] = fr["G"]
    except Exception as ex:
        out = {"key": key, "pattern": pattern, "error": repr(ex), "status": "UNRESOLVED"}
    out["t"] = time.perf_counter() - t0
    return out


def main():
    inp, outp = sys.argv[1], sys.argv[2]
    w = 1
    if "--workers" in sys.argv:
        w = int(sys.argv[sys.argv.index("--workers") + 1])
    w = max(1, min(w, MAXW))
    n = 9
    full = "--full" in sys.argv  # cross-check: all triangles, no covering shortcuts
    items = []
    for line in open(inp):
        d = json.loads(line)
        for i, fr in enumerate(d.get("factors", [])):
            if fr["status"] == "REALIZED":
                items.append((f"{d['pattern']}#{i}", d["pattern"], fr, n, full))
    done = set()
    if os.path.exists(outp):
        for l in open(outp):
            try:
                o = json.loads(l)
                if "error" not in o:
                    done.add(o["key"])
            except Exception:
                pass
    todo = [it for it in items if it[0] not in done]
    # big fields last (slowest)
    todo.sort(key=lambda it: len(it[2]["q"]))
    print(f"{len(items)} cores, {len(done)} done, {len(todo)} to do, workers={w}", flush=True)
    t0 = time.time()
    with open(outp, "a") as f, mp.get_context("spawn").Pool(w, maxtasksperchild=100) as pool:
        for k, res in enumerate(pool.imap_unordered(work, todo, chunksize=1)):
            f.write(json.dumps(res) + "\n")
            f.flush()
            if res.get("found") or res.get("error"):
                print("NOTE", res["key"], res.get("found"), res.get("error"), flush=True)
            if k % 50 == 0:
                print(f"{k+1}/{len(todo)} {time.time()-t0:.0f}s", flush=True)
    print("DONE", time.time() - t0, flush=True)


if __name__ == "__main__":
    main()
