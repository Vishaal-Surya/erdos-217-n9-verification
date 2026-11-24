"""Stage 5: independent (certificate-free) classification of all E2-unresolved excluded patterns."""
import json, os, time
import va_alg as V
fn = 'out/classify_e2unresolved.jsonl'
done = set()
if os.path.exists(fn):
    for l in open(fn): done.add(json.loads(l)['p'])
pats = [l.strip() for l in open('../../data/E2_unresolved_excluded.txt') if l.strip()]
t = time.time()
with open(fn, 'a') as f:
    for s in pats:
        if s in done: continue
        f.write(json.dumps({'p': s, 'va': V.classify(s)}) + '\n'); f.flush()
with open('out/run_all.log', 'a') as f:
    print(time.strftime('%H:%M:%S'), 'STAGE5 done', len(pats), round(time.time() - t), file=f)
