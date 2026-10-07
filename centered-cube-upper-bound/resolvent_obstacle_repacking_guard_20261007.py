from fractions import Fraction as F
from pathlib import Path
import json, hashlib
P=Path(__file__).parent

def solve(A,b):
    n=len(b); a=[list(A[i])+[b[i]] for i in range(n)]
    for j in range(n):
        k=next(k for k in range(j,n) if a[k][j]); a[j],a[k]=a[k],a[j]
        pivot=a[j][j];a[j]=[v/pivot for v in a[j]]
        for k in range(n):
            if k!=j:
                z=a[k][j];a[k]=[v-z*w for v,w in zip(a[k],a[j])]
    return [row[-1] for row in a]

def apply(A,u): return [sum((a*x for a,x in zip(row,u)),F(0)) for row in A]
rows=[]; checks=0
for n in [2,3,4]:
    N=2**n
    S=[[F(n,2) if i==j else -F(1,2) if (i^j).bit_count()==1 else F(0) for j in range(N)] for i in range(N)]
    profiles=[[F(i==0) for i in range(N)], [F(i in (0,N-1)) for i in range(N)], [F(1+(i%3),3) for i in range(N)]]
    for label,u in zip(['single spike','two separated spikes','positive full support'],profiles):
        sig=apply(S,u);kap=max(F(1),max(-x for x in sig));O=[x>0 for x in u]
        nu=[kap+sig[i] if O[i] else F(0) for i in range(N)]
        mu=[kap if O[i] else -sig[i] for i in range(N)]
        W=sum(nu,F(0))/N
        assert min(nu)>=0 and min(mu)>=0 and max(mu)<=kap
        assert sum(nu)==sum(mu);checks+=2
        for t in [F(1,4),F(1),F(4)]:
            A=[[F(i==j)+t*S[i][j] for j in range(N)] for i in range(N)]
            ru=solve(A,u);h=solve(A,sig)
            newnu=[kap+h[i] if O[i] else F(0) for i in range(N)]
            newmu=[kap if O[i] else ru[i]/t for i in range(N)]
            assert all(h[i]==(u[i]-ru[i])/t for i in range(N));checks+=1
            assert min(ru)>=0 and sum(ru)==sum(u);checks+=1
            assert min(h)>=-kap and sum(h)==0;checks+=1
            assert all(h[i]<=0 for i in range(N) if not O[i]);checks+=1
            assert sum(abs(x) for x in h)/N<=2*W;checks+=1
            assert min(newnu)>=0 and min(newmu)>=0 and max(newmu)<=kap;checks+=1
            assert all(newnu[i]-newmu[i]==h[i] for i in range(N));checks+=1
            assert sum(newnu)==sum(newmu) and sum(newnu)/N<=2*W;checks+=1
            rows.append({'n':n,'profile':label,'t':str(t),'W':str(W),'new_mass':str(sum(newnu)/N),'mass_ratio':str(sum(newnu)/(N*W))})
reg=P/'resolvent_obstacle_repacking_registration_20261007.json'
out={'status':'PASS','checks':checks,'rows':rows,'scope':'Exact finite Markov algebra diagnostic; not original-kernel numerical estimate or weak endpoint proof','registration_sha256':hashlib.sha256(reg.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(P/'resolvent_obstacle_repacking_results_20261007.json').write_text(json.dumps(out,indent=2))
print(json.dumps({'status':out['status'],'checks':checks,'cases':len(rows)}))
