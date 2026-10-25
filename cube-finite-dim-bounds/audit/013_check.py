#!/usr/bin/env python3
"""Independent algebra for tensor remainder, quotient cancellation, clipping repair."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import hashlib,json
import sympy as s
e=s.symbols('e',positive=True);R=s.symbols('R0:4',nonnegative=True);T=s.symbols('T0:4',nonnegative=True)
higher=s.prod(R[j]+e*T[j] for j in range(4))-s.prod(R)-e*sum(T[j]*s.prod(R[k] for k in range(4) if k!=j) for j in range(4))
subset=sum(e**len(S)*s.prod(T[j] if j in S else R[j] for j in range(4)) for m in range(2,5) for S in combinations(range(4),m))
assert s.expand(higher-subset)==0
V0,D,dv,j1,rj,m1,rm=s.symbols('V0 D dv j1 rj m1 rm',real=True);V=V0+D*e+dv
lhs=(1+j1*e+rj)*V-V0*(1+m1*e+rm)
rhs=(D+V0*(j1-m1))*e+dv+j1*e*(V-V0)+rj*V-V0*rm
assert s.expand(lhs-rhs)==0
bad_eps=F(1,32);bad=max(F(0),bad_eps-F(1,64));assert bad==F(1,64)>0
assert F(1,64)<F(1,16) and F(1,64)>F(1,256)
assert s.expand(16*e**2-e+s.Rational(1,64)-16*(e-s.Rational(1,32))**2)==0
# A new fixed rational fixture for the epsilon=s^4 selection argument only.
g=F(1,10);C=F(7,2);K=F(2,3);t=F(1,2);steps=0
while C*t**4+K*t**3>g/2:t/=2;steps+=1
q=t**4;margin=g*q-C*q*q-K*t**7
assert margin>=g*q/2>0
out={'status':'independent exact identities PASS','engine':'SymPy '+s.__version__+' and stdlib Fraction','tensor_higher_coefficients_identity':True,'full_mass_Jacobian_cancellation_identity':True,'clipping_counterexample':{'cap':'1/16','epsilon':str(bad_eps),'actual_clipped_value':str(bad),'claimed_B0_invalid':True,'guarded_smaller_cap':'1/256','strict_gap_slack':str(F(1,64)-F(1,256)),'full_cap_true_B':'16','repair_square':'16*(e-1/32)^2'},'selection_fixture':{'g':str(g),'Cnet':str(C),'K':str(K),'halving_steps':steps,'fourth_root':str(t),'epsilon':str(q),'positive_margin':str(margin)},'scope':'identities and one abstract fixture only, not a real n4 lower certificate','author_imports':False,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('013_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
