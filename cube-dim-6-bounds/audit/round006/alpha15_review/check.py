#!/usr/bin/env python3
"""Independent alpha=3/2 symbolic and Arb audit; zero kernel/MC evaluations."""
import hashlib,json
from pathlib import Path
import sympy as s
from flint import arb,ctx
ctx.prec=256
HERE=Path(__file__).resolve().parent
SRC=HERE.parents[2]/'upper/round006/alternative'
a=s.Rational(3,2);d=s.Integer(6);h=d-a;q=d/h
g,c,I=s.symbols('g c I',positive=True)
cost=h/a*I*c**(-a/h)/(g-c)
assert s.simplify(s.diff(cost,c).subs(c,g/4))==0
assert s.simplify(cost.subs(c,g/4)-4**q*I/g**q)==0
exponents=[k-(d-k)/a for k in range(1,7)]
assert exponents[:2]==[-s.Rational(7,3),-s.Rational(2,3)]
aa=arb(3)/2
M=(arb(2)/3).gamma()/(aa*arb.pi())
A=(arb(5)/2).gamma()/(arb.pi()*(arb.pi()/3).sin()**(arb(5)/2))
Cold1=A*M**5*arb(64)**(arb(5)/2)
Cold2=31*A**2*M**4*arb(64)**(arb(5)/2)*(62*arb(6).sqrt())**(arb(5)/2)
Dold=3*Cold2+arb(3)/4*Cold1
C1=A*M**5*arb(32)**(arb(5)/2)
C2=31*A**2*M**4*arb(32)**(arb(5)/2)*(31*arb(6).sqrt())**(arb(5)/2)
D=arb(15)/4*C1**(arb(1)/5)*C2**(arb(4)/5)
Bcap=(6*D)**(arb(4)/3)*arb.pi()**2/28
Bcap_old=(6*Dold)**(arb(4)/3)*arb.pi()**2/28
Tmult=(C1/C2)**(arb(3)/5)
C0=M**6/160
Iglobal=arb(1)/48*(arb(4)*C0/3)**(arb(1)/3)
target=arb(1)/10**6
delta=(target/Bcap)**(arb(3)/7)
assert Dold/D>400
assert Bcap>10**9
assert delta<arb(1)/10**6
out=dict(status='PASS_CONDITIONAL_ANALYTIC_INTERFACE_AND_CAP_REFINEMENT',
 no_new_C6_bound=True,kernel_evaluations=0,MC_samples=0,
 exact=dict(a=str(a),h=str(h),q=str(q),rank_exponents=list(map(str,exponents)),
            cost='4^(4/3)*I/g^(4/3)',cap_delta_exponent='7/3',
            D_optimized='(15/4)*C1^(1/5)*C2^(4/5)',
            global_I_upper='(1/48)*(4*C0/3)^(1/3)'),
 arb256={k:v.str(55) for k,v in dict(M=M,A=A,C1=C1,C2=C2,D_old=Dold,D_optimized=D,
    D_improvement=Dold/D,cap_coefficient_old=Bcap_old,cap_coefficient_optimized=Bcap,
    cap_coefficient_improvement=Bcap_old/Bcap,time_cut_multiplier=Tmult,
    global_I_upper=Iglobal,delta_for_allocated_cap_I_1e_minus6=delta).items()},
 input_hashes={name:hashlib.sha256((SRC/name).read_bytes()).hexdigest() for name in ['report.md','sources.json']},
 script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(HERE/'check_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
