"""Driver: E2 reconstruction of every F pattern (resumable, explicit worker count <= 2).

python src/extend_run_recon.py data/E1_F.txt.gz data/F_recon.jsonl --workers 2 [--limits a,b,c]
Appends one JSON line per pattern; on restart, patterns already in the output are skipped.
"""
import gzip, json, os, sys, time, multiprocessing as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
MAXW = 2  # user CPU cap for builder D


def work(s, limits=(10000, 100000, 100000)):
    import extend_recon as R
    try:
        return R.reconstruct_pattern(s, limits=tuple(limits))
    except Exception as ex:  # recorded as UNRESOLVED, never as excluded
        return {"pattern": s, "e2": "ERROR", "error": repr(ex), "factors": []}


def main():
    inp, outp = sys.argv[1], sys.argv[2]
    w = 1
    if "--workers" in sys.argv:
        w = int(sys.argv[sys.argv.index("--workers") + 1])
    w = max(1, min(w, MAXW))
    limits = (10000, 100000, 100000)
    if "--limits" in sys.argv:  # e.g. --limits 100000,1000000,1000000 (redo of GB-limit cases)
        limits = tuple(int(x) for x in sys.argv[sys.argv.index("--limits") + 1].split(","))
    import functools
    op = gzip.open if inp.endswith(".gz") else open
    pats = [l.split()[0] for l in op(inp, "rt") if l.strip()]
    done = set()
    if os.path.exists(outp):
        for l in open(outp):
            try:
                done.add(json.loads(l)["pattern"])
            except Exception:
                pass
    todo = [p for p in pats if p not in done]
    print(f"{len(pats)} patterns, {len(done)} done, {len(todo)} to do, workers={w}", flush=True)
    t0 = time.time()
    with open(outp, "a") as f, mp.get_context("spawn").Pool(w, maxtasksperchild=200) as pool:
        for k, res in enumerate(pool.imap_unordered(functools.partial(work, limits=limits), todo, chunksize=1)):
            f.write(json.dumps(res) + "\n")
            if k % 50 == 0:
                f.flush()
                print(f"{k+1}/{len(todo)} {time.time()-t0:.0f}s", flush=True)
    print("DONE", time.time() - t0, flush=True)


if __name__ == "__main__":
    main()
