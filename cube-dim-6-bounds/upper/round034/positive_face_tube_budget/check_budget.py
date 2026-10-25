#!/usr/bin/env python3
"""F034: ten fixed symbolic/scalar configurations, zero kernels/grids."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,math,os,sys,time
import sympy as sp
import flint
from flint import arb,ctx
ctx.prec=512
HERE=Path(__file__).resolve().parent
start=time.time()
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
inp=json.loads((HERE/'input_manifest.json').read_text())
assert start-inp['started_unix']<1200
assert not (HERE/'results.json').exists()
for p,h in inp['input_sha256'].items():assert sha(p)==h,p
oldroot=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round008/off_strata_integration')
old=json.loads((oldroot/'results.json').read_text())
def q(a,b=1):return arb(a)/b
def save(x):
 return {'ball':str(x),'lower_exact':str(x.lower().fmpq()),'upper_exact':str(x.upper().fmpq())}
def check_overlap(x,s):assert x.overlaps(arb(s))
def cap(V,e):
 R=1/(3*e**3)-1024/e+32768*(pi/2-(32*e).atan())
 return K*e**q(7,6)+V**q(5,6)*Cf*R

gram=[]
for k in range(1,6):
 m=6-k
 # Columns of the projection of (0,-u) into the A zero-sum space.
 P=sp.eye(k)-sp.ones(k,k)/k
 J=-P[:,1:k]
 M=J.T*J
 det=M.det() if k>1 else sp.Integer(1)
 assert det==sp.Rational(1,k)
 assert M==sp.eye(k-1)-sp.ones(k-1,k-1)/k
 gram.append({'k':k,'m':m,'Gram':str(M),'Gram_determinant':str(det),
              'face_to_residual_inverse_Jacobian_squared':k,
              'mean_B_Jacobian_squared':m,'cone_fraction':f'1/{k}',
              'tube_factor_squared':str(sp.Rational(m,k)),
              'partition_multiplicity':math.comb(5,k-1)})
new_alg=sum(sp.binomial(5,k-1)*sp.sqrt(sp.Rational(6-k,k)) for k in range(1,6))
old_alg=sum(sp.binomial(5,k-1)*sp.sqrt(k*(6-k)) for k in range(1,6))
assert sp.simplify(new_alg-(10+10*sp.sqrt(2)+2*sp.sqrt(5)))==0
assert sp.simplify(old_alg-3*new_alg)==0
pi=arb.pi();v4=pi**2/2;v5=8*pi**2/15
Cnew=v4*(10+10*arb(2).sqrt()+2*arb(5).sqrt())
Cold=3*Cnew;C5=26*arb(2)**q(5,2)
C1=4*arb(6).sqrt()*pi**2/45
assert (arb(6).sqrt()/6*v5-C1).contains(0)
check_overlap(Cold,old['Ctube_delta4']['ball']);check_overlap(C5,old['Cmixed_delta5']['ball'])
Cf=1/(243*pi**6);A=Cf**q(1,6)*pi**(-q(5,6));K=q(6,7)*A
check_overlap(Cf,old['Cf']['ball']);check_overlap(K,old['K']['ball'])
target=q(1,11250)
rows=[]
for den,oldrow in zip((16,8,4,2),old['layers']):
 assert oldrow['delta']==f'1/{den}'
 d=q(1,den);voldraw=Cold*d**4+C5*d**5;vnewraw=Cnew*d**4
 vold=voldraw if voldraw<1 else arb(1)
 vnew=vnewraw if vnewraw<1 else arb(1)
 eold=(Cf*pi*vold)**q(1,5)
 enew=(Cf*pi*vnew)**q(1,5)
 Uold=cap(vold,eold);Usame=cap(vnew,eold);Unew=cap(vnew,enew)
 check_overlap(vold,oldrow['tube_volume_upper']['ball'])
 check_overlap(Uold,oldrow['cap_norm_budget_exact_time_denominator']['ball'])
 Rnew=1/(3*enew**3)-1024/enew+32768*(pi/2-(32*enew).atan())
 # Exact lower bound on integral of the pointwise min of the TWO allocated envelopes.
 Lmethod=K*enew**q(7,6)/(1+1024*enew**2)+vnew**q(5,6)*Cf*Rnew
 assert Lmethod>target
 if den!=2:assert Unew<Uold and Usame<Uold
 else:assert (Unew-Uold).contains(0)
 rows.append({'delta':f'1/{den}','old_volume':save(vold),'new_volume':save(vnew),
              'old_raw':save(voldraw),'new_raw':save(vnewraw),
              'old_over_new_volume_ratio':save(vold/vnew),
              'old_epsilon':save(eold),'new_same_selection_rule_epsilon':save(enew),
              'old_allocated_budget':save(Uold),'new_at_old_epsilon':save(Usame),
              'new_allocated_budget':save(Unew),'new_budget_over_target':save(Unew/target),
              'best_two_envelope_method_lower_bound':save(Lmethod),
              'method_lower_over_target':save(Lmethod/target),
              'method_still_fails_target':True})
# The smallest registered tube already forces the same two-envelope budget above target.
assert rows[0]['method_still_fails_target']
# A rational safe lower bound, without interpreting this as a true G-norm lower bound.
Lrat=Fraction(rows[0]['method_lower_over_target']['lower_exact'])
Ls=Fraction(Lrat.numerator*10**6//Lrat.denominator,10**6)
out={'status':'PASS_NEW_TUBE_UPPER_CONSTANT_BUT_OLD_CAP_METHOD_STILL_FAILS',
     'pid':os.getpid(),'python':sys.version,'python_flint':flint.__version__,'precision_bits':ctx.prec,
     'exact_geometry_scalar_configurations':10,'Gram_classes':gram,
     'Cnew_exact':'(pi^2/2)*(10+10sqrt(2)+2sqrt(5))','Cold_equals_3_Cnew':True,
     'Cnew':save(Cnew),'Cold':save(Cold),'old_C5':save(C5),'S1_cap_coefficient':save(C1),
     'nearest_unsigned_reduction_preexisting':'A012 theorem section2',
     'new_contribution':'maximum-coordinate cone fraction1/k plus removal of redundant signed tubes',
     'rows':rows,'smallest_layer_method_ratio_safe_lower':str(Ls),
     'target_conditional_on_unproved_global_g':'1/11250',
     'new_kernel_calls':0,'grid_runs':0,'optimization_runs':0,'matlab_calls':0,
     'new_I_certificate':False,'new_C6_bound':False,
     'source_sha256':sha(Path(__file__)),'registration_sha256':sha(HERE/'registration.md'),
     'input_manifest_sha256':sha(HERE/'input_manifest.json'),
     'elapsed_seconds':time.time()-start,'elapsed_since_registration':time.time()-inp['started_unix']}
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'pid':out['pid'],'elapsed_seconds':out['elapsed_seconds'],
 'Cnew':str(Cnew),'S1cap':str(C1),'method_safe_lower':str(Ls),
 'rows':[{'delta':r['delta'],'new_volume':r['new_volume']['ball'],
          'new_budget':r['new_allocated_budget']['ball'],'target_ratio':r['new_budget_over_target']['ball'],
          'method_lower_ratio':r['method_lower_over_target']['ball']} for r in rows]},indent=2))
