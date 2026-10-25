from pathlib import Path
from itertools import product,combinations
from fractions import Fraction as F
import hashlib,json
import sympy as s
from flint import arb,ctx
ctx.prec=192
def ball(x):
 x=F(x);return arb(x.numerator)/x.denominator
def pw(x,a,b):return (x.log()*a/b).exp()
def ends(x):return [str(x.lower().fmpq()),str(x.upper().fmpq())]
out=[]
for n,tau,target in [(3,F(1,3),F(57,100000)),(4,F(0),F(19,100000))]:
 body=[(1,)+v for v in product((-1,1),repeat=n-1)];V=body[:];weights=[F(1)]*len(body)
 if tau:
  for i,j in combinations(range(n),2):
   for sg in (-1,1):
    v=[0]*n;v[i]=1;v[j]=sg;V.append(tuple(v));weights.append(tau)
 active=[i for i,v in enumerate(V) if v[0]];k=len(active);S=F(0);bases=[]
 for i in active:
  best=max((abs(int(s.Matrix([V[i],*rest]).det())),rest) for rest in combinations(body,n-1));det,rest=best;assert det>0;S+=weights[i]/det;bases.append({'direction':V[i],'completion':rest,'det':det})
 A=sum(weights[i] for i in active);cn=k*k*S;assert cn==(F(256,3) if n==3 else F(32))
 records=[]
 for eps in [F(1,1000),F(1,2000)]:
  E=pw(arb(2),-(n-1)**2,n)*ball(F(n,n+1))*pw(ball(cn)/arb.pi()**n,1,n)*pw(ball(A)/arb.pi(),n-1,n)*pw(ball(eps),n+1,n)
  root=pw(ball(target),n-1,n)-E;assert root>0;allowed=pw(root,n,n-1)
  records.append({'epsilon':str(eps),'E_interval':ends(E),'E_display':str(E),'allowed_I_H_interval':ends(allowed),'allowed_I_H_display':str(allowed)})
 out.append({'n':n,'tau':str(tau),'k':k,'A':str(A),'C_numerator_over_pi_n':str(cn),'bases':bases,'I_target':str(target),'budgets':records})
p=Path(__file__).resolve().parents[1];(p/'results/round011_face_time_budget.json').write_text(json.dumps({'status':'primary analytic lemma plus exact finite check; no I certificate yet; independent audit pending','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'cases':out},indent=2)+'\n');print(json.dumps([{k:v for k,v in z.items() if k!='bases'} for z in out],indent=2))
