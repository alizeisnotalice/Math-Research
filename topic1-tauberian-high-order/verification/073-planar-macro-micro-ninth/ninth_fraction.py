from fractions import Fraction as F
from pathlib import Path
import json, os, sys
R=11
class S:
    def __init__(self, coeff=0):
        self.c=list(coeff[:R+1])+[F(0)]*max(0,R+1-len(coeff)) if isinstance(coeff,list) else [F(coeff)]+[F(0)]*R
    def __add__(self,o):
        o=conv(o); return S([x+y for x,y in zip(self.c,o.c)])
    __radd__=__add__
    def __neg__(self): return S([-x for x in self.c])
    def __sub__(self,o): return self+-conv(o)
    def __rsub__(self,o): return conv(o)+-self
    def __mul__(self,o):
        o=conv(o); return S([sum((self.c[j]*o.c[i-j] for j in range(i+1)), F(0)) for i in range(R+1)])
    __rmul__=__mul__
    def inv(self):
        assert self.c[0]
        out=[1/self.c[0]]
        for i in range(1,R+1): out.append(-sum((self.c[j]*out[i-j] for j in range(1,i+1)),F(0))/self.c[0])
        return S(out)
    def __truediv__(self,o): return self*conv(o).inv()
    def __rtruediv__(self,o): return conv(o)*self.inv()
    def __pow__(self,n):
        assert n>=0
        out=S(1)
        for _ in range(n): out=out*self
        return out

def conv(o): return o if isinstance(o,S) else S(o)
u=S([F(0),F(1)]); a=4*u**2; v=2*u; g=1/(1+2*u)
k=F(16,3)*g**2+F(640,81)*g**3*u
lp=k**3/(32*g)-5*k**4*u/(768*g**2)+13*k**5*u**2/(8192*g**3)-33*k**6*u**3/(81920*g**4)
beta=4*u**2-16*u**3
lm=lp-beta*lp+F(5,192)*F(16,3)**4*u**3

def ratio(ll):
    d=1-k**2*u**6; rho=g-2*ll*u**6
    return d/(d-a*rho**2)
qp=ratio(lp); qm=ratio(lm)
ln=1/(1-a*g**2); p7=ln+F(4,27)*a*v**6-F(176,243)*a*v**7
res_p=qp-p7; res_m=qm-p7; inc=qm-qp
assert all(x==0 for x in res_m.c[:10])
assert res_m.c[10]/(4*2**8)==F(17638,6561)
assert res_m.c[11]/(4*2**9)==-F(23295344,2657205)
assert res_p.c[10]/(4*2**8)==F(15694,6561)
assert res_p.c[11]/(4*2**9)==-F(17871584,2657205)
assert inc.c[10]/(4*2**8)==F(8,27)
assert inc.c[11]/(4*2**9)==-F(496,243)
# Independent rational integral of (1+r)*r on [0,1].
cap=F(1,2)+F(1,3); assert cap==F(5,6)
assert 4*4*cap/F(64)==F(5,24)
data={"python":sys.version,"pid":os.getpid(),"truncation_degree":R,
      "lambda_pure":[str(x) for x in lp.c[:4]],
      "lambda_micro":[str(x) for x in lm.c[:4]],
      "pure_residual":[str(x) for x in res_p.c],
      "micro_residual":[str(x) for x in res_m.c],
      "micro_increment":[str(x) for x in inc.c],
      "normalized_micro_ninth":str(res_m.c[11]/(4*2**9)),
      "assertions_passed":9}
Path(__file__).with_suffix('.json').write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:data[k] for k in ('python','pid','lambda_pure','lambda_micro','normalized_micro_ninth','assertions_passed')},ensure_ascii=False,indent=2))
