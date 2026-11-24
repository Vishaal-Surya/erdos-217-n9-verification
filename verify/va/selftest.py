import random, itertools, time, gzip, json
from fractions import Fraction as Fr
import va_alg as V
random.seed(1)
# 1. polynomial identities at random rational points (pattern with all 15 labels distinct is not allowed (r<=6),
# so evaluate pair-variable polys directly)
def evalp(poly, y):
    t=0
    for k,v in poly.items():
        p=Fr(v)
        for i in k: p*=y[i]
        t+=p
    return t
for trial in range(20):
    P=[(Fr(random.randint(-50,50),random.randint(1,9)),Fr(random.randint(-50,50),random.randint(1,9))) for _ in range(6)]
    y=[(P[i][0]-P[j][0])**2+(P[i][1]-P[j][1])**2 for (i,j) in V.PAIRS]
    assert all(evalp(m,y)==0 for m in V.MINORS)
    for Q,poly in V.DET4.items():
        B=[[1,P[i][0],P[i][1],P[i][0]**2+P[i][1]**2] for i in Q]
        import sympy
        dB=sympy.Matrix(B).det()
        assert evalp(poly,y)==-4*dB**2
    for (i,j,k),poly in V.TRI.items():
        u=(P[j][0]-P[i][0],P[j][1]-P[i][1]); w=(P[k][0]-P[i][0],P[k][1]-P[i][1])
        assert evalp(poly,y)==4*(u[0]*w[1]-u[1]*w[0])**2
print("identities OK; distinct minors", len({tuple(sorted(m.items())) for m in V.MINORS}))
# 2. timing on a few patterns
F=gzip.open('../../data/E1_F.txt.gz','rt').read().split('\n')[:5]
for line in F:
    s=line.split()[0]; t=time.time(); print(line, V.classify(s), round(time.time()-t,3))
