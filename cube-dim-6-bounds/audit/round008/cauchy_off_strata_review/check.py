#!/usr/bin/env python3
"""Independent exact geometry + ONE off-S2 scalar budget; no kernel/MC."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product,combinations
from collections import Counter
import hashlib,json,time
import sympy as s
from flint import arb,ctx
ctx.prec=256
START=time.monotonic()
HERE=Path(__file__).resolve().parent
B=HERE.parents[2]/'upper/round008/off_strata_integration'
def ar(q):
    q=F(q);return arb(q.numerator)/q.denominator
def dot(a,b):return sum(u*v for u,v in zip(a,b))
def ser(a):return a.str(60)
V=[(1,)+v for v in product((-1,1),repeat=5)]
assert all(sum(v[i]*v[j] for v in V)==32*(i==j) for i in range(6) for j in range(6))

# Gram-Schmidt rather than B's expanded 3x3 determinant formula.
counts=Counter()
for v,w,z in combinations(V,3):
    d=dot(v,w);e=dot(v,z);f=dot(w,z)
    n1=F(6);n2=6-F(d*d,6)
    proj=F(6*e*e-2*d*e*f+6*f*f,36-d*d)
    n3=6-proj
    determinant=n1*n2*n3
    assert determinant.denominator==1 and n2>0 and n3>0
    counts[int(determinant)]+=1
assert counts==Counter({64:480,96:1920,128:480,160:640,192:1440})
assert sum(counts.values())==4960
S3=sum((6*n/arb(d).sqrt() for d,n in counts.items()),arb(0))

# Enumerate ALL pair planes and independently classify positive-face coverage.
classes=Counter();same_subsets=set();mixed_endpoints=set();eligible_by_h=Counter()
for v,w in combinations(V,2):
    equal=[i for i in range(1,6) if v[i]==w[i]]
    J=tuple(i for i in range(1,6) if v[i]!=w[i])
    if any(v[i]==-1 for i in equal):
        classes['distance_at_least_1_over_sqrt2']+=1
        continue
    eligible_by_h[len(J)]+=1
    if len({v[i] for i in J})==1:
        classes['same_sign_segment']+=1;same_subsets.add(J)
    else:
        classes['mixed_sign_plane']+=1;mixed_endpoints.add(J)
assert classes==Counter({'distance_at_least_1_over_sqrt2':375,'same_sign_segment':31,'mixed_sign_plane':90})
assert len(same_subsets)==31 and len(mixed_endpoints)==26
assert same_subsets=={J for k in range(1,6) for J in combinations(range(1,6),k)}
assert mixed_endpoints=={J for k in range(2,6) for J in combinations(range(1,6),k)}

# Canonical exact pair distance formula and transverse determinant, for h=1..5.
geometry=[]
ys=s.symbols('y1:6',real=True);x=s.Matrix([1,*ys])
for h in range(1,6):
    v=s.ones(6,1);w=s.Matrix([1]+[-1]*h+[1]*(5-h))
    W=v.row_join(w);P=s.eye(6)-W*(W.T*W).inv()*W.T
    d2=(x.T*P*x)[0]
    j=ys[:h];u=[1-a for a in ys[h:]]
    proposed=sum(a*a for a in u)-sum(u)**2/s.Integer(6-h)+sum(a*a for a in j)-sum(j)**2/s.Integer(h)
    assert s.expand(d2-proposed)==0
    dim=5-h;Q=s.eye(dim)-s.ones(dim,dim)/s.Integer(6-h)
    determinant=Q.det() if dim else s.Integer(1)
    assert determinant==s.Rational(1,6-h)
    geometry.append({'h':h,'fixed_block_dimension':dim,'fixed_block_metric_determinant':str(determinant),
                     'axis_length_squared':h,'tube_cross_section_factor_squared':h*(6-h)})

# Fourier radial integral and exact antiderivative checks.
t=s.symbols('t',positive=True)
primitive=1/(3*t**3)-1024/t+32768*(s.pi/2-s.atan(32*t))
assert s.simplify(s.diff(primitive,t)+t**-4/(1+1024*t*t))==0
assert s.Rational(5,6)+1==s.Rational(11,6)
assert s.Rational(-4,6)==s.Rational(-2,3)
assert s.Rational(5,6)-1>-1 and s.Rational(5,6)-2<-1
pi=arb.pi();c=32/arb(3).sqrt()
Q3=1/(pi*pi*c**3)
B3=(32/pi)**3*S3*Q3
assert ar(F(46141420,10**6))<B3<ar(F(46141421,10**6))
Ctube=pi*pi/2*sum((arb(int(s.binomial(5,h)))*arb(h*(6-h)).sqrt() for h in range(1,6)),arb(0))
Cmixed=26*arb(2)**ar(F(5,2))
assert Ctube<arb('423.617312') and Cmixed<arb('147.078211')

# Exactly one layer and one short-time epsilon, registered in advance.
delta=arb(1)/16;eps=arb(1)/1000;target=arb(1)/11250
volume=Ctube*delta**4+Cmixed*delta**5
assert volume<1
Cf=1/(243*pi**6)
A=Cf**ar(F(1,6))*pi**(-ar(F(5,6)))
K=arb(6)/7*A
whole=K*eps**ar(F(7,6))
off=arb(6)/11*B3**ar(F(1,6))*pi**(-ar(F(5,6)))*delta**(-ar(F(2,3)))*eps**ar(F(11,6))
split=(Cf*pi*volume)**ar(F(1,5))
tail=Cf*(1/(3*split**3)-1024/split+32768*(pi/2-(32*split).atan()))
cap=K*split**ar(F(7,6))+volume**ar(F(5,6))*tail
cap_loose=K*split**ar(F(7,6))+Cf/3*volume**ar(F(5,6))*split**-3
# Lower bound on the minimum INTEGRATED ENVELOPE, not on the actual cap norm.
# The common time denominator keeps the crossover at `split`.
optimal_envelope_lower=K*split**ar(F(7,6))/(1+1024*split**2)+volume**ar(F(5,6))*tail
assert off<whole and off<ar(F(7990283,10**12))
assert volume<ar(F(6604152,10**9))
assert cap>11*target and cap<ar(F(1051196,10**9))
assert cap<cap_loose
assert optimal_envelope_lower>5*target
saved=json.loads((B/'results.json').read_text())
for key,value in [('B3',B3),('ordered_jacobian_sum3',S3),('Ctube_delta4',Ctube),('Cmixed_delta5',Cmixed)]:
    assert value.overlaps(arb(saved[key]['ball']))
layer=next(x for x in saved['layers'] if x['delta']=='1/16')
row=next(x for x in saved['rows'] if x['delta']=='1/16' and x['epsilon']=='1/1000')
assert volume.overlaps(arb(layer['tube_volume_upper']['ball']))
assert cap.overlaps(arb(layer['cap_norm_budget_exact_time_denominator']['ball']))
assert off.overlaps(arb(row['short_time_off_norm_upper']['ball']))
paths=[B/'registration.md',B/'check_budget.py',B/'results.json']
if (B/'report.md').exists():paths.append(B/'report.md')
out={'status':'PASS_GEOMETRY_AND_ONE_ALLOCATED_BUDGET; NO_C6_BOUND',
 'kernel_evaluations':0,'MC_samples':0,'matlab_calls':0,'geometric_layers':1,'scalar_short_time_pairs':1,
 'exact_gramdet3_unordered_counts':dict(sorted(counts.items())),
 'exact_positive_face_pair_classes':dict(classes),'same_sign_segments':len(same_subsets),'mixed_endpoint_boxes':len(mixed_endpoints),
 'canonical_tube_geometry':geometry,'first_three_time_denominators':[1024,992,960],
 'constants':{k:ser(v) for k,v in dict(S3=S3,B3=B3,Ctube=Ctube,Cmixed=Cmixed).items()},
 'one_budget':{'delta':'1/16','epsilon':'1/1000',**{k:ser(v) for k,v in dict(tube_volume_upper=volume,
   off_short_norm_upper=off,global_short_norm_upper=whole,off_over_global=off/whole,
   cap_split=split,cap_allocated_norm_upper=cap,cap_allocated_upper_over_target=cap/target,cap_loose=cap_loose,
   optimum_two_envelope_budget_lower=optimal_envelope_lower,
   optimum_two_envelope_budget_lower_over_target=optimal_envelope_lower/target).items()},
   'cap_allocated_upper_exceeds_target':True,'actual_cap_norm_lower_bound_claimed':False},
 'strict_rational_upper_bounds':{'tube_volume':'825519/125000000','off_short_norm':'7990283/1000000000000',
    'cap_allocated_norm_budget':'262799/250000000'},
 'derivative_conclusion':'Pointwise O(1) heat times t^-m is not integrable at zero for m>=1. The L6/5 envelope is t^(5/6), so m=1 DOES give local W^(1,6/5); m>=2 is not covered at zero. No nonexistence claim about actual derivatives.',
 'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'elapsed_seconds':time.monotonic()-START,'exit_code':0}
(HERE/'check_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
