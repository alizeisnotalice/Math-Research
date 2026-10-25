#!/usr/bin/env python3
"""One exact translation and two same-point comparisons. No kernel evaluation."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,time
from flint import arb,ctx,fmpz
ctx.prec=384
START=time.monotonic()
HERE=Path(__file__).resolve().parent
B=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round010/short_minimum_diagnostic')
OWN=HERE.parents[1]
P1=B/'point1.json';P2=B/'point2.json';R1=B/'point1_results.json';R2=B/'batch2_results.json'
p1=json.loads(P1.read_text());p2=json.loads(P2.read_text())
b1=json.loads(R1.read_text());b2=json.loads(R2.read_text())
x1=[F(n,p1['denominator']) for n in p1['candidate_numerators']]
x2=[F(n,p2['denominator']) for n in p2['candidate_numerators']]
assert b1['point_numerators']==p1['candidate_numerators'] and b1['denominator']==p1['denominator']
assert b2['point_numerators']==p2['candidate_numerators'] and b2['denominator']==p2['denominator']
assert x1[0]==x2[0]==1 and x1[1]==x2[1]
for x in (x1,x2):assert all(x[i]>x[i+1]>0 for i in range(5))
# Six distinct positive coordinates imply both points are off S2 (indeed S3).
h=[b-a for a,b in zip(x1,x2)];r=max(map(abs,h));imax=[i+1 for i,a in enumerate(h) if abs(a)==r]
assert imax==[3]
def ar(q):
    q=F(q);return arb(q.numerator)/q.denominator
G2=arb(fmpz(b2['G']['G_times_pi5_numerator']))/arb(fmpz(b2['G']['G_times_pi5_denominator']))/arb.pi()**5
assert G2.overlaps(arb(b2['G']['G_arb384']['ball']))
T32=F(4,125);T64=F(8,125)
assert F(b2['T'])==T32
H32=arb(b2['H_arb384'])
hrow=next(row for row in b1['H_rows'] if F(row['T'])==T64)
H64_at_x1=arb(hrow['H_arb384'])
u=ar(r/T64)/12
phi=(u+(1+u*u).sqrt())**12
assert phi>1
H64_lower_formula=H64_at_x1/phi
K32=G2-H32
K64_upper_formula=G2-H64_lower_formula
theta=arb(1)/4096
n32=K32/ar(T32);n64_upper=K64_upper_formula/ar(T64)
assert K32>0 and K64_upper_formula>0
assert K32.overlaps(arb(b2['K_T']['ball']))
assert n32<theta and n64_upper<theta
assert n64_upper<n32
assert n32<ar(F(243309,10**9))<theta
uniform_cost_lower=1/(64*n32)
assert uniform_cost_lower>ar(F(642187,10000))
gap_exact=F(1,4096)-F(243309,10**9)
assert gap_exact>0
# Check source hashes for the actual immutable point files used in B's records.
for results,pointfile in ((b1,P1),(b2,P2)):
    records=[v for v in results['inputs'] if v['path']==str(pointfile)]
    assert len(records)==1 and records[0]['sha256']==hashlib.sha256(pointfile.read_bytes()).hexdigest()
paths=[P1,P2,R1,R2,B/'checked_results.json',B/'check_results.py',B/'report.md',
       OWN/'round007/sparse_shift_bound.md',OWN/'round009/short_potential/report.md']
out={'status':'PASS_FINITE_POSITIVE_MIXTURE_TRANSFER_AND_COMMON_POINT_EXCLUSION',
 'kernel_evaluations':0,'MC_samples':0,'matlab_calls':0,'fixed_geometry_shifts':1,'target_comparisons':2,
 'point1_exact':[str(x) for x in x1],'point2_exact':[str(x) for x in x2],
 'shift_exact':[str(a) for a in h],'maxnorm_shift_exact':str(r),'maxnorm_coordinate_1_based':3,
 'point2_strict_chamber_and_off_S2':True,
 'T_support':['4/125','8/125'],
 'scalar_formulas':{key:value.str(85) for key,value in dict(G2=G2,Phi=phi,
  H64_at_x1=H64_at_x1,H64_x2_lower_formula=H64_lower_formula,
  K32_at_x2=K32,K64_x2_upper_formula=K64_upper_formula,
  normalized_K32_at_x2=n32,normalized_K64_x2_upper_formula=n64_upper,
  threshold_minus_normalized_K32=theta-n32,
  threshold_minus_normalized_K64_upper_formula=theta-n64_upper,
  uniform_model_cost_lower=uniform_cost_lower).items()},
 'H64_x2_was_not_evaluated':True,'K64_x2_upper_formula_is_not_a_two_sided_value_interval':True,
 'common_point_normalized_mixture_strict_rational_upper':'243309/1000000000',
 'gap_below_1_over_4096_at_least_exact':str(gap_exact),
 'all_support_mixtures_model_cost_strict_lower':'642187/10000',
 'probability_weights_quantifier':'every nonnegative pair with sum1, including both endpoints',
 'new_T_or_other_supports_excluded':False,'new_C6_bound':False,
 'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'elapsed_seconds':time.monotonic()-START,'exit_code':0}
(HERE/'check_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
