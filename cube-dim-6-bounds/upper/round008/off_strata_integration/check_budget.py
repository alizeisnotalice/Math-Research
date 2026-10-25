#!/usr/bin/env python3
"""B008 bounded scalar audit. Zero body32 kernel evaluations and no MATLAB."""
import hashlib, itertools, json, math, os, platform, sys, time
from collections import Counter
from fractions import Fraction
from pathlib import Path
import flint
from flint import arb, ctx
ctx.prec = 384
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
started=time.monotonic()

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ball(x): return {'ball':str(x),'lower':str(x.lower()),'upper':str(x.upper())}
def q(a,b=1): return arb(a)/b

inputs=[
 ROOT/'AGENTS.md', ROOT/'cube-dim-6-bounds/STATE.md',
 ROOT/'cube-dim-6-bounds/upper/round006/integral/report.md',
 ROOT/'cube-dim-6-bounds/upper/round007/mixed_time_error/report.md',
 ROOT/'cube-dim-6-bounds/audit/round006/report.md',
 ROOT/'cube-dim-6-bounds/audit/round007/sparse_shift_bound.md',
 ROOT/'cube-dim-6-bounds/upper/round007/g_time_moments/report.md',
 ROOT/'cube-dim-6-bounds/upper/theory/001_n6_body_frame_and_singular_strata.md',
 ROOT/'cube-dim-6-bounds/upper/round005/cover/010_product_inequality_whole_face.md',
 HERE/'registration.md',HERE/'check_budget.py']
manifest={'inputs':[{'path':str(p),'sha256':sha(p)} for p in inputs]}
(HERE/'input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

V=[(1,)+s for s in itertools.product((-1,1),repeat=5)]
assert all(sum(v[i]*v[j] for v in V)==32*(i==j) for i in range(6) for j in range(6))
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def gramdet3(a,b,c):
 x,y,z=dot(a,b),dot(a,c),dot(b,c)
 return 216+2*x*y*z-6*(x*x+y*y+z*z)
counts=Counter(gramdet3(*tr) for tr in itertools.combinations(V,3))
assert sum(counts.values())==4960 and min(counts)==64
assert set(counts)=={64,96,128,160,192}
# Every ordered deletion triple is counted once; no replacement by unordered sum.
s1=arb(32)/arb(6).sqrt()
s2=sum((arb(32*math.comb(5,h))/(2*arb(h*(6-h)).sqrt()) for h in range(1,6)),arb(0))
s3=sum((arb(6*n)/arb(g).sqrt() for g,n in sorted(counts.items())),arb(0))
pi=arb.pi(); c=32/arb(3).sqrt()
Q5=2/(pi**3*c**5); Q4=3/(4*pi**2*c**4); Q3=1/(pi**2*c**3)
B1=(32/pi)*s1*Q5
B2=(32/pi)**2*s2*Q4
B3=(32/pi)**3*s3*Q3
Cf=1/(243*pi**6); A=Cf**q(1,6)*pi**(-q(5,6)); K=q(6,7)*A
Bfull=A*arb(32)**(-q(7,6))*pi/(2*(7*pi/12).sin())
target=q(1,11250)
D6=Fraction(73354040202609,16384000000)
D12=Fraction(2354965600826387008780486926466707,8589934592000000000000)
d6=q(D6.numerator,D6.denominator); d12=q(D12.numerator,D12.denominator)
# Positive-face S2 tube volume, valid 0<delta<1/sqrt(2).
Ctube=(pi**2/2)*sum((arb(math.comb(5,h))*arb(h*(6-h)).sqrt() for h in range(1,6)),arb(0))
Cmixed=26*arb(2)**q(5,2)
# Since d1>=d2>=delta and |x|>=1, p<=B3/delta^4.
# Interpolate with positive-orthant marginal <=t/pi; integrate in time.
# Uoff=(6/11)*B3^(1/6)*pi^(-5/6)*delta^(-2/3)*epsilon^(11/6).
rows=[]; layers=[]
for den in [16,8,4,2]:
 delta=q(1,den)
 vol_raw=Ctube*delta**4+Cmixed*delta**5
 vol=vol_raw if vol_raw<1 else arb(1)
 # Analytic optimizer of K eps^(7/6)+(Cf/3)V^(5/6)eps^-3.
 epsstar=(Cf*pi*vol)**q(1,5)
 short=K*epsstar**q(7,6)
 long_loose=Cf/3*vol**q(5,6)*epsstar**(-3)
 # exact integrated t^-4/(1+1024 t^2), at the same explicit epsstar.
 tailint=Cf*(1/(3*epsstar**3)-1024/epsstar+32768*(pi/2-(32*epsstar).atan()))
 cap=short+vol**q(5,6)*tailint
 cap_loose=short+long_loose
 layers.append({'delta':f'1/{den}','tube_volume_raw_upper':ball(vol_raw),'tube_volume_upper':ball(vol),
 'epsilon_analytic_split':ball(epsstar),'cap_norm_budget_exact_time_denominator':ball(cap),
 'cap_norm_budget_loose_closed_form':ball(cap_loose),'cap_budget_over_target':ball(cap/target),
 'cap_budget_certifies_target':bool(cap<target),
 'pointwise_off_constant_Cdelta':ball(B3/delta**4),
 'modulus_crossover_T':ball((Cf*delta**4/B3)**q(1,4)),
 'modulus_uniform_radius_necessary_for_target_tenth':ball((4*target/10)/(5*(B3/delta**4))),
 'uniform_cell_count_lower_for_modulus_budget':ball((1-vol)/(2*((4*target/10)/(5*(B3/delta**4))))**5)})
 for ed in [1000,100,32]:
  eps=q(1,ed)
  whole=K*eps**q(7,6)
  off=q(6,11)*B3**q(1,6)*pi**(-q(5,6))*delta**(-q(2,3))*eps**q(11,6)
  # First and mixed-order moments using the near-zero off-stratum scalar envelope.
  # This is an upper budget for J_m after epsilon (extend t^(5/6) envelope to infinity).
  C_off=B3**q(1,6)*pi**(-q(5,6))*delta**(-q(2,3))
  m6=C_off*eps**(q(11,6)-6)/(6-q(11,6))
  m12=C_off*eps**(q(11,6)-12)/(12-q(11,6))
  rows.append({'delta':f'1/{den}','epsilon':f'1/{ed}',
  'short_time_off_norm_upper':ball(off),'short_time_whole_face_upper':ball(whole),
  'allocated_off_over_whole':ball(off/whole),'off_upper_beats_whole':bool(off<whole),
  'allocated_M6_upper':ball(m6),'allocated_M12_upper':ball(m12),
  'allocated_one_axis_Gauss3_relative_budget_h_2epsilon':ball(d6/46080*(2*eps)**6*m6/target),
  'allocated_two_axis_mixed_Gauss3_relative_budget_h_2epsilon':ball(d12/(46080**2)*(2*eps)**12*m12/target)})
assert len(rows)==12 and len(layers)==4
# Exact derivative coefficient identities, m=1,6,12, no kernel calls.
derivatives={str(m):str(Fraction(math.factorial(63+m),math.factorial(63)*32**m)) for m in [1,6,12]}
# Exact Gaussian mixed constants are imported as mathematical constants, not rerun kernels.
D6=Fraction(73354040202609,16384000000)
D12=Fraction(2354965600826387008780486926466707,8589934592000000000000)
assert D12>D6*D6
out={'status':'PASS_EXPLICIT_BOUND_AND_FINITE_BUDGET_DIAGNOSTIC_NO_C6_IMPROVEMENT',
'pid':os.getpid(),'python':sys.version,'python_executable':sys.executable,'platform':platform.platform(),
'python_flint':flint.__version__,'arb_precision_bits':ctx.prec,
'new_kernel_points':0,'geometric_distance_layers':4,'scalar_combinations':12,'matlab_calls':0,
'gramdet3_unordered_counts':dict(sorted(counts.items())),
'quotient_coercivity':'32/sqrt(3) times Euclidean norm, from second/fourth moments and Holder',
'quotient_Q5':ball(Q5),'quotient_Q4':ball(Q4),'quotient_Q3':ball(Q3),
'ordered_jacobian_sum1':ball(s1),'ordered_jacobian_sum2':ball(s2),'ordered_jacobian_sum3':ball(s3),
'B1':ball(B1),'B2':ball(B2),'B3':ball(B3),'Cf':ball(Cf),'A':ball(A),'K':ball(K),
'Ctube_delta4':ball(Ctube),'Cmixed_delta5':ball(Cmixed),
'baseline_full_G_norm_budget':ball(Bfull),'conditional_target':str(Fraction(1,11250)),
'baseline_budget_over_target':ball(Bfull/target),
'derivative_coefficients_Cm':derivatives,'layers':layers,'rows':rows,
'proof_status':'Bounds proved in report.md; scalar arithmetic certified by Arb384. These allocated budgets do not lower-bound the real norm or derivatives.',
'input_manifest_sha256':sha(HERE/'input_manifest.json'),'source_sha256':sha(__file__),
'elapsed_seconds':time.monotonic()-started,'exit_code':0}
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'pid':out['pid'],'elapsed_seconds':out['elapsed_seconds'],
'B1':str(B1),'B2':str(B2),'B3':str(B3),'Ctube':str(Ctube),
'layers':[{'delta':x['delta'],'volume':x['tube_volume_upper']['ball'],'cap_over_target':x['cap_budget_over_target']['ball']} for x in layers],
'rows':[{'delta':x['delta'],'epsilon':x['epsilon'],'off_over_whole':x['allocated_off_over_whole']['ball']} for x in rows]},indent=2))
