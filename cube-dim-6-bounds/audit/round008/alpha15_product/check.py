#!/usr/bin/env python3
"""F008 zero-kernel independent audit: one C007 tail case and alpha15 face products."""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import hashlib,json
import sympy as s
from flint import arb,ctx
ctx.prec=256
HERE=Path(__file__).resolve().parent
C7=HERE.parents[2]/'upper/round007/alpha15_deterministic'
def ar(x):
    x=Q(x);return arb(x.numerator)/x.denominator
def dot(x,y):return sum(a*b for a,b in zip(x,y))
V=[(1,)+v for v in product((-1,1),repeat=5)]
x=tuple(map(Q,('1','4/5','3/5','2/5','1/5','0')))
M=(arb(2)/3).gamma()/((arb(3)/2)*arb.pi())
A=(arb(5)/2).gamma()/(arb.pi()*(arb.pi()/3).sin()**(arb(5)/2))

# Independent fixed-point geometry and just ONE of C007's sixteen tail pairs.
saved=json.loads((C7/'truncation_results.json').read_text())
geom=[];sumH=arb(0);sumB=arb(0)
for i,v in enumerate(V):
    z=tuple(x[j]-v[j]*x[0] for j in range(1,6));ds=[];es=[]
    for j,w in enumerate(V):
        if i==j:continue
        vw=dot(v,w);vx=dot(v,x);wx=dot(w,x)
        ds.append(dot(x,x)-(6*vx*vx-2*vw*vx*wx+6*wx*wx)/(36-vw*vw))
        u=tuple(w[k]-v[k] for k in range(1,6))
        es.append(dot(z,z)-dot(z,u)**2/dot(u,u))
    d2,e2=min(ds),min(es);assert d2>0 and e2>0
    rec=saved['geometric_distances_exact'][i]
    assert d2==Q(rec['distance_to_pair_planes_squared']) and e2==Q(rec['quotient_distance_to_lines_squared'])
    CH=31*A**2*M**4*(31*arb(6).sqrt()/ar(d2).sqrt())**5
    CB=ar(Q(31*30,32))*A**2*M**3*(60*arb(5).sqrt()/ar(e2).sqrt())**5
    sumH+=arb(10)/3*CH**ar(Q(9,10))*(M**6/160)**ar(Q(1,10))
    sumB+=arb(10)/3*CB**ar(Q(9,10))*(M**5/160)**ar(Q(1,10))/(2*arb.pi().sqrt())
    geom.append([str(d2),str(e2)])
assert min(Q(v[0]) for v in geom)==Q(4,25)
assert min(Q(v[1]) for v in geom)==Q(6,25)
ell=arb(1)/32;R=arb(10)**18;r=64;p=arb(7)/10
small=sumH*ell**r*arb(2)**(arb(4)*r/3)*(1+arb(4)*r/3).gamma()/(1+arb(r)).gamma()
large=sumB*(1-4*p/3).gamma()/(1-p).gamma()*R**(-p-arb(1)/2)
assert small+large<arb('9.592e-11')

# Symbolic normalization, ellipsoid determinant, exponents and moment interpolation.
sv=s.Matrix([-1,1,-1,1,-1]);T=(-sv).row_join(s.eye(5));v=s.Matrix([1,*list(sv)])
assert T*v==s.zeros(5,1)
assert T*T.T==s.eye(5)+sv*sv.T
assert (T*T.T).det()==6
assert -s.Rational(7,3)/4+s.Rational(3,4)==s.Rational(1,6)
assert -s.Rational(7,3)/4-s.Rational(5,3)/4-s.Rational(2,3)*s.Rational(3,4)==-s.Rational(3,2)
variants=[]
for fourth in (Q(3),Q(8,3)):
    c=32*ar(fourth)**(-arb(1)/4)
    C5rad=(arb(10)/3).gamma()/(18*arb.pi()**3)*c**(-arb(10)/3)
    Cq=C5rad/arb(6).sqrt()
    C6=c**-4/(16*arb.pi()**3)
    variants.append(dict(fourth_moment_factor=str(fourth),c=c.str(60),C5rad=C5rad.str(60),Cquotient=Cq.str(60),C6rad=C6.str(60)))

# Use stronger c6 and determinant6 in all final face budgets.
c=32*ar(Q(8,3))**(-arb(1)/4)
Cq=(arb(10)/3).gamma()/(18*arb.pi()**3)*c**(-arb(10)/3)/arb(6).sqrt()
C6=c**-4/(16*arb.pi()**3)
assert C6<M*Cq
B0=A*arb(32)**(arb(5)/2)*Cq
D=A/M*arb(32)**(arb(5)/2)
E=M*arb(32)**(-arb(5)/3)
tc=(C6/B0)**(arb(3)/5)
t0=(E/A)**(arb(3)/5)
assert tc<t0
a0=B0**(arb(1)/4)*A**(arb(3)/4)
a1=C6**(arb(1)/4)*A**(arb(3)/4)
a2=C6**(arb(1)/4)*E**(arb(3)/4)
first=arb(6)/7*a0*tc**(arb(7)/6)
second=arb(4)/3*a1*(t0**(arb(3)/4)-tc**(arb(3)/4))
third=2*a2*t0**(-arb(1)/2)
full=first+second+third
assert full<ar(Q(1002681,10**9))
assert full**(arb(4)/3)<ar(Q(100358,10**9))
def prefix(e):
    if e<tc:return arb(6)/7*a0*e**(arb(7)/6)
    if e<t0:return first+arb(4)/3*a1*(e**(arb(3)/4)-tc**(arb(3)/4))
    return first+second+2*a2*(t0**(-arb(1)/2)-e**(-arb(1)/2))
tails=[]
for ep,tu,hu in ((Q(1,10000),1416468,1001263767),(Q(1,1000),20790903,981889332),
                 (Q(1,100),218481669,784198566),(Q(1,50),385603959,617076276)):
    val=prefix(ar(ep));assert val>0 and val<full
    assert val<ar(Q(tu,10**12)) and full-val<ar(Q(hu,10**12))
    tails.append(dict(epsilon=str(ep),T_epsilon_norm_upper=val.str(60),H_epsilon_norm_upper=(full-val).str(60),
        T_norm_strict_rational_upper=str(Q(tu,10**12)),H_norm_strict_rational_upper=str(Q(hu,10**12))))
paths=[C7/'report.md',C7/'truncation_budget.py',C7/'truncation_results.json']
C8=C7.parent.parent/'round008'/'alpha15_coercivity'
paths += [C8/'report.md',C8/'coercivity_check.py',C8/'coercivity_results.json']
out=dict(status='PASS_C007_ONE_TAIL_CASE_AND_ALPHA15_WHOLE_FACE_PRODUCTS',
 kernel_evaluations=0,MC_samples=0,scalar_scenarios=7,new_C6_bound=False,
 C007=dict(ell='1/32',R='10^18',negative_moment_order=64,positive_moment_order='7/10',
   all_32_exact_geometry_rows_match=True,sum_H=sumH.str(60),sum_B=sumB.str(60),
   small_tail=small.str(60),large_tail=large.str(60),total_tail=(small+large).str(60),
   interior_integral_certified=False),
 coercivity_variants=variants,
 normalization=dict(N=32,quotient_Jacobian='1',ellipsoid_determinant=6,
   positive_orthant_L1='f_(32t)(1)/32 <= min(A*t, M*32^(-5/3)*t^(-2/3))',
   erroneous_At_over_32_rejected=True,short_norm_time_exponent='1/6',long_norm_time_exponent='-3/2',
   coarse_At_harmonic_full_time_exponent='-1/4 (not integrable)'),
 final_constants={k:v.str(60) for k,v in dict(M=M,A=A,c6=c,Cquotient=Cq,C6rad=C6,
   B0=B0,harmonic_B=2*B0,harmonic_D=D,E=E,tc=tc,t0=t0,
   short_time_K=arb(6)/7*a0,full_G_norm_upper=full,full_I_upper=full**(arb(4)/3)).items()},
 strict_rational_bounds={'G_norm':'1002681/1000000000','I':'50179/500000000',
   'C007_deleted_parameter_tail':'1199/12500000000000'},
 finite_epsilon_budgets=tails,
 input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
 script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(HERE/'check_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
