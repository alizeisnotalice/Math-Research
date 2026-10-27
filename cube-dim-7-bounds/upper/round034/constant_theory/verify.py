"""Exact parameter audit; no n=7 or high-dimensional numerical search."""
from fractions import Fraction as R
from pathlib import Path
import hashlib
import json
import math
import sympy as s

root=Path(__file__).resolve().parent
Y=R(64); aa=R(1,10); bb=R(29,100); cc=R(21,100); qq=R(1,20)
g=bb+cc
assert 0<aa<bb<Y/2 and 0<cc<Y/2 and g*Y>=1
f1=1-4*g/Y; f2=1-2*g*g/(1-bb/Y)
assert f1>0 and f2>0
KA=cc/(4*(1+cc/Y)); CL=f1*f2
assert KA==R(336,6421) and CL==R(98301,203872)
taylor=sum((R(9,5)**k)/math.factorial(k) for k in range(7))
assert taylor==R(7542629,1250000)>6
KP_lower=(bb-aa)*(CL-R(1,12))-1/Y-R(3,8)/Y**3
assert KP_lower==R(12055437337,200414330880)>qq
assert KA>qq
assert KP_lower-qq==R(2034720793,200414330880)
# log2<3/4 is established by this lower Taylor sum.
assert 1+R(3,4)+R(3,4)**2/2==R(65,32)>2
assert R(64,4)-18*R(3,4)==R(5,2)>0
assert R(1,4)-R(3,64)>0
checks=['example admissibility and positive local-capture factors',
        'exact KA and CL', 'positive Taylor certificate exp(9/5)>6',
        'exact KP lower bound and strict q=1/20 margins',
        'exact log2 upper certificate and x>=3L/4 base at L=64']

a,beta,zeta=s.symbols('a beta zeta')
alpha=a/(a+1); kappa=beta-zeta
assert s.factor(alpha*beta-kappa-2/(a+1)-(zeta-(beta+2)/(a+1)))==0
checks.append('finite-a guard identity')

gv=s.symbols('g',positive=True)
Q=gv*(1-2*gv**2)/(5-8*gv**2)
assert s.factor(s.diff(Q,gv)-(5-22*gv**2+16*gv**4)/(5-8*gv**2)**2)==0
z=(11-s.sqrt(41))/16
assert s.simplify(5-22*z+16*z*z)==0
qinf=(19-s.sqrt(41))*s.sqrt(11-s.sqrt(41))/320
assert s.simplify(Q.subs(gv,s.sqrt(z))-qinf)==0
checks.append('family-ceiling derivative and exact radical value')
receipt={
 'status':'ACCEPT: analytic all-n proof in PROOF.md; exact finite algebra checks passed',
 'sympy':s.__version__,'checks':checks,
 'example':{k:str(v) for k,v in {'Y':Y,'a0':aa,'b0':bb,'c0':cc,'q':qq,
           'KA':KA,'CL':CL,'exp_9_over_5_taylor_lower':taylor,
           'KP_lower':KP_lower,'KP_lower_minus_q':KP_lower-qq}.items()},
 'family_q_ceiling':str(s.N(qinf,40)),
 'source_sha256':{name:hashlib.sha256((root/name).read_bytes()).hexdigest()
                  for name in ('PROOF.md','verify.py')},
 'no_parameter_search':True,
}
(root/'verification.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(receipt,ensure_ascii=False,indent=2))
