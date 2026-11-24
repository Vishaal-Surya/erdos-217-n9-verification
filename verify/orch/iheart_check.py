"""Exact (C,Q) count for the forum's IheartTesla lattice config, integer arithmetic.
Lattice (a,b) -> (2a+b, b) * (1/2, sqrt3/2) scaling: use X=2a+b, Y=b with metric X^2+3Y^2 (=4*|p|^2).
Collinear: det in (X,Y) (linear image). Concyclic: det rows (1,X,Y,X^2+3Y^2) (G-metric lift)."""
import itertools
from collections import Counter
from fractions import Fraction
P=[(4,-2),(4,-3),(2,1),(3,-1),(2,2),(-1,0),(0,0),(1,-1),(3,1)]
Q=[(2*a+b,b) for a,b in P]
def det(M):
    M=[[Fraction(x) for x in r] for r in M]; n=len(M); d=Fraction(1)
    for i in range(n):
        p=next((r for r in range(i,n) if M[r][i]!=0),None)
        if p is None: return 0
        if p!=i: M[i],M[p]=M[p],M[i]; d=-d
        d*=M[i][i]
        for r in range(i+1,n):
            f=M[r][i]/M[i][i]; M[r]=[M[r][k]-f*M[i][k] for k in range(n)]
    return d
d2=lambda p,q:(p[0]-q[0])**2+3*(p[1]-q[1])**2
sizes=sorted(Counter(d2(p,q) for p,q in itertools.combinations(Q,2)).values())
C=[t for t in itertools.combinations(range(9),3) if det([[Q[t[1]][k]-Q[t[0]][k] for k in (0,1)],[Q[t[2]][k]-Q[t[0]][k] for k in (0,1)]])==0]
Qd=[q for q in itertools.combinations(range(9),4) if det([[1,Q[i][0],Q[i][1],Q[i][0]**2+3*Q[i][1]**2] for i in q])==0]
print('sizes',sizes,'C',len(C),C,'Q',len(Qd),Qd)
