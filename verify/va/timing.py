import gzip, json, time, random
import va_alg as V
random.seed(2)
L=gzip.open('../../data/E1_excluded_certs.jsonl.gz','rt').read().splitlines()
smp=random.sample(L,200)
t=time.time(); bad=0
for line in smp:
    d=json.loads(line); ok,dm,msg=V.replay(d['pattern'],d['cert'])
    if not ok or not dm: bad+=1; print(d['pattern'],ok,dm,msg)
print('replay 200', time.time()-t, 'bad', bad)
loc=gzip.open('../../data/patterns6_local.txt.gz','rt').read().split()
smp=random.sample(loc,1000)
t=time.time(); n=sum(V.affine(s)[0] for s in smp); print('affine 1000', time.time()-t, n)
t=time.time()
for line in random.sample(L,30):
    d=json.loads(line); c=V.classify(d['pattern']); assert c[0]=='EXCL',(d,c)
print('classify 30 excl', time.time()-t)
