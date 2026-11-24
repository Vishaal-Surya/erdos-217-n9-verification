"""V-a sequential driver (1 process). Stages, each resumable via its output file."""
import gzip, json, os, random, sys, time
import va_alg as V
D = '../../data/'
O = 'out/'
def log(*a):
    with open(O + 'run_all.log', 'a') as f:
        print(time.strftime('%H:%M:%S'), *a, file=f, flush=True)

# Stage 1: affine on my own local list
if not os.path.exists(O + 'va_patterns6_affine.txt'):
    t = time.time(); surv = []; rej = {}
    for s in open(O + 'va_patterns6_local.txt').read().split():
        ok, why, _ = V.affine(s)
        if ok: surv.append(s)
        else:
            k = why.split()[0].split('(')[0]; rej[k] = rej.get(k, 0) + 1
    open(O + 'va_patterns6_affine.txt', 'w').write('\n'.join(surv) + '\n')
    b = set(gzip.open(D + 'patterns6_affine.txt.gz', 'rt').read().split())
    a = set(surv)
    byr = {}
    for s in surv: r = max(s); byr[r] = byr.get(r, 0) + 1
    log('STAGE1 affine', len(a), 'byr', sorted(byr.items()), 'rej', rej, 'data', len(b), 'mine-data', len(a - b), 'data-mine', len(b - a), 'sec', round(time.time() - t))
    for s in sorted(a ^ b)[:50]: log('  DIFF', s, s in a)

aff = set(open(O + 'va_patterns6_affine.txt').read().split())
F = {}
for line in gzip.open(D + 'E1_F.txt.gz', 'rt').read().splitlines():
    if line.strip(): p = line.split(); F[p[0]] = line
U = {}
for line in gzip.open(D + 'E1_U.txt.gz', 'rt').read().splitlines():
    if line.strip(): p = line.split(); U[p[0]] = line

# Stage 2: independent classification of the 4,637 retained
fn = O + 'classify_retained.jsonl'
done = set()
if os.path.exists(fn):
    for line in open(fn): done.add(json.loads(line)['p'])
if len(done) < len(F) + len(U):
    t = time.time()
    with open(fn, 'a') as f:
        for s in sorted(F) + sorted(U):
            if s in done: continue
            c = V.classify(s)
            f.write(json.dumps({'p': s, 'builder': 'F' if s in F else 'U', 'meta': (F.get(s) or U.get(s)), 'va': c}) + '\n'); f.flush()
    log('STAGE2 done', round(time.time() - t))

# Stage 3: replay all certificates
fn = O + 'replay_all.jsonl'
start = 0
if os.path.exists(fn):
    start = sum(1 for _ in open(fn))
lines = gzip.open(D + 'E1_excluded_certs.jsonl.gz', 'rt').read().splitlines()
if start < len(lines):
    t = time.time()
    with open(fn, 'a') as f:
        for i in range(start, len(lines)):
            d = json.loads(lines[i])
            s = d['pattern']
            ok, dm, msg = V.replay(s, d['cert'])
            f.write(json.dumps({'i': i, 'p': s, 'ok': ok, 'dims': dm, 'msg': msg, 'ntrace': len(d['cert']['trace'])}) + '\n')
            if i % 5000 == 0: f.flush(); log('STAGE3', i, round(time.time() - t))
    log('STAGE3 done', round(time.time() - t))

# Stage 4: independent classification of random excluded (no certificates used)
fn = O + 'classify_excluded_sample.jsonl'
random.seed(217)
pats = [json.loads(l)['pattern'] for l in lines]
smp = random.sample(pats, 6000)
# stress add-on: every excluded pattern whose certificate needed >=2 division steps
stress = [json.loads(l)['pattern'] for l in lines if len(json.loads(l)['cert']['trace']) >= 2]
todo = list(dict.fromkeys(smp + stress)); SS = set(stress)
done = set()
if os.path.exists(fn):
    for line in open(fn): done.add(json.loads(line)['p'])
t = time.time()
with open(fn, 'a') as f:
    for s in todo:
        if s in done: continue
        c = V.classify(s)
        f.write(json.dumps({'p': s, 'va': c, 'stress': s in SS}) + '\n'); f.flush()
log('STAGE4 done', len(todo), round(time.time() - t))
