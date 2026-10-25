#!/usr/bin/env python3
"""F009: five exact Gram identities and one scalar short-barrier certificate."""
from pathlib import Path
from itertools import product,combinations
from fractions import Fraction as F
import hashlib,json,time
import sympy as s
from flint import arb,ctx
ctx.prec=384
START=time.monotonic()
HERE=Path(__file__).resolve().parent
A=HERE.parents[2]/'upper/round009/short_potential_barrier'
def ar(q):
    q=F(q);return arb(q.numerator)/q.denominator
V=[s.Matrix([1,*u]) for u in product((-1,1),repeat=5)]
frame=sum((v*v.T for v in V),s.zeros(6,6))
assert frame==32*s.eye(6)
geometries=[]
for k in range(1,6):
    v=s.ones(6,1);w=s.Matrix([1]+[-1]*k+[1]*(5-k))
    B=v.row_join(w);G=B.T*B;D=G.det();R=G.inv()*B.T;Q=s.eye(6)-B*R
    assert D==4*k*(6-k)
    assert R*R.T==G.inv()
    assert (R*R.T)[0,0]==(R*R.T)[1,1]==s.Rational(6,D)
    assert Q*v==s.zeros(6,1) and Q*w==s.zeros(6,1)
    assert Q*frame*Q==32*Q
    eigen=sorted(G.eigenvals())
    assert eigen==sorted(set([2*k,12-2*k]))
    geometries.append({'k':k,'Gram':[[int(q) for q in row] for row in G.tolist()],
      'det':int(D),'min_Gram_eigenvalue':int(min(eigen)),
      'pseudoinverse_row_norm_squared':str((R*R.T)[0,0]),
      'pair_count':int(16*s.binomial(5,k)),'projected_tight_frame_verified':True})
assert sum(g['pair_count'] for g in geometries)==496

# Exact fourth-moment improvement and radial normalization in dimension four.
x=s.symbols('x0:6',real=True)
assert s.expand(6*sum(z**4 for z in x)-sum(z*z for z in x)**2-
                sum((x[i]**2-x[j]**2)**2 for i,j in combinations(range(6),2)))==0
S3=2*s.pi**2
a_exact=S3*s.gamma(4)/((2*s.pi)**4*32**4)
b_exact=S3*s.gamma(6)/(8*(2*s.pi)**4*(8*s.sqrt(6))**6)
assert s.simplify(b_exact/a_exact-s.Rational(5,108))==0
D,tau,X=s.symbols('D tau X',positive=True)
f=D**s.Rational(3,2)/(tau*tau*D+6*X*X)**2
assert s.simplify(s.diff(s.log(f),D)-(18*X*X-tau*tau*D)/(2*D*(tau*tau*D+6*X*X)))==0

# ONE registered L/T/delta, independent Gamma evaluation (not A's eight-L loop).
pi=arb.pi();sphere=2*pi*pi
c=8*arb(6).sqrt();crest_c=7*arb(6).sqrt()
a=sphere*arb(4).gamma()/((2*pi)**4*32**4)
b=sphere*arb(6).gamma()/(8*(2*pi)**4*c**6)
Crest=sphere*arb(4).gamma()/((2*pi)**4*crest_c**4)
Ctail=30*(2/pi)*(30*arb(6).sqrt())*Crest
kappa=(4*b/a).sqrt()
assert kappa.overlaps((arb(5)/27).sqrt())
L=arb(3)/4;T=arb(4)/125;delta=arb(1)/80000
tau=a*L/(4*Ctail)
assert tau<T
X=arb(6).sqrt()+L
assert tau*tau*36 <18*X*X
Ds=[20,32,36,32,20]
ds=[]
for determinant in Ds:
    Bsq=6*X*X/determinant
    ds.append(a/(2*pi*pi*arb(determinant).sqrt()*(tau*tau+Bsq)**2))
assert ds[0].overlaps(ds[4]) and ds[0]<ds[1] and ds[1]<ds[2]
d=ds[0];coef=d/kappa;offset=kappa/tau
assert coef>arb(1)/5000000000 and offset<40000
assert kappa*delta<tau
barrier=d*(1/(kappa*delta)-1/tau)
assert barrier>ar(F(1655621,200000000000))
simple=ar(F(1,5000000000))*(1/delta-40000)
assert simple.overlaps(ar(F(1,125000)))
margin=F(1,125000)-F(1,128000)
assert margin==F(3,16000000) and margin>0
saved=json.loads((A/'results.json').read_text())
row=next(r for r in saved['rows'] if r['L']=='3/4')
for name,value in [('tau_arb384',tau),('density_coefficient_d_arb384',d),
                  ('barrier_c_arb384',coef),('barrier_C_arb384',offset)]:
    assert value.overlaps(arb(row[name]))
assert barrier.overlaps(arb(saved['K_T_lower_at_distance_threshold_arb384']))
paths=[A/'theorem.md',A/'results.json',A/'registration.md',A/'check_barrier.py',A/'verify_barrier.py',A/'checked_results.json']
out={'status':'PASS_UNIFORM_ALL_CUBE_S2_NEIGHBORHOOD_BARRIER; COMPLEMENT_OPEN',
 'scalar_cases':1,'L':'3/4','T':'4/125','delta':'1/80000','canonical_Gram_cases':5,
 'kernel_evaluations':0,'MC_samples':0,'matlab_calls':0,'exact_pair_geometries':geometries,
 'b_over_a':'5/108','minimum_D':20,
 'constants':{key:value.str(90) for key,value in dict(a=a,b=b,Crest=Crest,Ctail=Ctail,kappa=kappa,
      tau=tau,d=d,barrier_coefficient=coef,barrier_offset=offset,
      actual_lower_budget_at_delta=barrier).items()},
 'rounded_coefficient_strict_lower':'1/5000000000','rounded_offset_strict_upper':'40000',
 'global_distance_barrier':'K_(4/125)(x) >= (1/5000000000)*(1/dist(x,S2)-40000)_+ for every x in the cube off S2; K=infinity on S2',
 'closed_tube_uniform_lower':'1/125000','short_potential_threshold':'1/128000','exact_positive_margin':str(margin),
 'whole_cube_minimum_certified':False,'new_C6_bound':False,
 'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'elapsed_seconds':time.monotonic()-START,'exit_code':0}
(HERE/'check_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
