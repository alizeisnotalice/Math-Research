"""Independent tiny exact positive-polynomial ratios plus numerical McD checks."""
from pathlib import Path
from fractions import Fraction as F
import itertools,math,json
import numpy as np
HERE=Path(__file__).resolve().parent
n=3;hidden=list(itertools.product([0,1],repeat=n))
z=[[[[F(1)+F((i+b+1)*(1 if h else -1),20),
        F(3,4)+F((i+2*b+1)*(1 if h else -1),30)]for b in range(2)]for h in range(2)]for i in range(n)]
def component(obs,b):
    return sum((math.prod(z[i][obs[i]][b][bits[i]]for i in range(n))
                for bits in hidden if(sum(bits)<=1 if b==0 else sum(bits)>=2)),F(0))
values=[[component(o,b)for o in hidden]for b in range(2)]
records=[];ratios=0
for b in range(2):
    cs=[];rs=[];qs=[]
    for i in range(n):
        actual=max(max(z[i][0][b][s],z[i][1][b][s])/min(z[i][0][b][s],z[i][1][b][s])for s in range(2))
        q=max(abs(z[i][0][b][s]-z[i][1][b][s])/((z[i][0][b][s]+z[i][1][b][s])/2)for s in range(2))
        assert actual<=1/(1-q)
        for j,o in enumerate(hidden):
            other=list(o);other[i]=1-other[i];index=hidden.index(tuple(other))
            assert values[b][j]/values[b][index]<=actual;ratios+=1
        cs.append(math.log(float(actual)));qs.append(q);rs.append(q*q/4)
    S=math.fsum(c*c for c in cs);q0=max(qs);r0=max(rs)
    C=4*(1+float(r0))/(1-float(q0))**2
    surrogate=C*math.fsum(math.log1p(float(r))for r in rs)
    assert S<=surrogate+1e-13
    h=np.log(np.array(list(map(float,values[b]))));center=h-float(h.mean());mgf=[]
    for t in [-10.,-1.,.5,1.,10.]:
        v=t*center;lhs=float(np.logaddexp.reduce(v)-math.log(len(v)));upper=t*t*S/8
        assert lhs<=upper+1e-12;mgf.append({'t':t,'lhs':lhs,'upper':upper})
    records.append({'b':b,'S':S,'surrogate':surrogate,'q0':str(q0),'r0':str(r0),'MGF_checks':mgf})
uppers=[sum(v)/len(v)*F(101,100)for v in values];M=sum(uppers);tails=[]
for multiplier in [F(101,100),F(11,10),F(3,2),F(2)]:
    lam=M*multiplier;actual=F(0)
    for j,o in enumerate(hidden):
        event=sum(v[j]for v in values)>lam
        union=any(v[j]>lam*u/M for v,u in zip(values,uppers))
        assert not event or union
        actual+=F(event,len(hidden))
    gap=math.log(float(multiplier));bound=min(1.,math.fsum(math.exp(-2*gap*gap/r['S'])for r in records))
    assert float(actual)<=bound+1e-12
    tails.append({'multiplier':str(multiplier),'tail_exact':str(actual),'upper':bound})
out={'n':n,'hidden_states':len(hidden),'exact_coordinate_ratio_checks':ratios,
     'component_records':records,'sum_tail_checks':tails,
     'scope':'Arbitrary independent two-point law; this is not the uniform-bin variance/3 test. No source or MC replay.'}
p=HERE/'mcd_checks.json'
if p.exists():raise RuntimeError('Refusing replacement')
p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
