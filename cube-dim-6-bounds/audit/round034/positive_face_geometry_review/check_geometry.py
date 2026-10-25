#!/usr/bin/env python3
"""A034: nine fixed independent geometry/analytic checks; no kernel or budget rerun."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,math,os,sys,time,datetime
import sympy as s
import flint
from flint import arb,ctx
HERE=Path(__file__).resolve().parent
start=time.time()
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
inp=json.loads((HERE/'input_manifest.json').read_text())
for p,h in inp['source_sha256'].items():
 assert sha(p)==h,p
assert not (HERE/'results.json').exists()
ctx.prec=384
checks=[]
for k in range(1,6):
 m=6-k
 P=s.eye(k)-s.ones(k)/k
 # Independently form columns from e1=sum(u)/k, ei=sum(u)/k-ui.
 J=s.Matrix(k,k-1,lambda i,j:s.Rational(1,k)-(1 if i==j+1 else 0))
 M=J.T*J
 assert M==s.eye(k-1)-s.ones(k-1)/k
 detA=M.det() if k>1 else s.Integer(1)
 assert detA==s.Rational(1,k)
 inv=s.Matrix(k-1,k,lambda i,j:(1 if j==0 else 0)-(1 if j==i+1 else 0))
 assert inv*J==s.eye(k-1)
 assert s.ones(1,k)*J==s.zeros(1,k-1)
 # Rational zero-sum B-basis; divide its metric determinant out to obtain
 # orthogonal dz. This checks the mean Jacobian without assuming an ON basis.
 Z=s.Matrix(m,m-1,lambda i,j:(1 if i==j else 0)-(1 if i==m-1 else 0))
 D=s.ones(m,1).row_join(Z)
 zero_det=(Z.T*Z).det() if m>1 else s.Integer(1)
 mean_jac_sq=s.cancel((D.T*D).det()/zero_det)
 assert mean_jac_sq==m
 if k>1:
  Q=s.eye(k);Q.row_swap(0,k-1)
  assert Q.T*Q==s.eye(k) and Q*P==P*Q
 checks.append({'configuration':k,'k':k,'m':m,'face_Gram':str(M),
 'det_face_Gram':str(detA),'inverse_face_Jacobian_squared':k,
 'B_mean_orthogonal_Jacobian_squared':str(mean_jac_sq),
 'cone_fraction':f'1/{k}','unsigned_partition_count':math.comb(5,k-1),'pass':True})
new=sum(s.binomial(5,k-1)*s.sqrt(s.Rational(6-k,k)) for k in range(1,6))
old=sum(s.binomial(5,k-1)*s.sqrt(k*(6-k)) for k in range(1,6))
assert s.simplify(new-(10+10*s.sqrt(2)+2*s.sqrt(5)))==0
assert s.simplify(old-3*new)==0
assert sum(math.comb(5,k-1) for k in range(1,6))==31
checks.append({'configuration':6,'new_factor':str(s.simplify(new)),
 'old_equals_3_new':True,'unsigned_planes':31,'pass':True})
P6=s.eye(6)-s.ones(6)/6
J6=-P6[:,1:6]
assert (J6.T*J6).det()==s.Rational(1,6)
w=s.Matrix([1,-1,0,0,0,0]);assert (w.T*w)[0]==2
K1=s.sqrt(6)/6*(8*s.pi**2/15)
assert s.simplify(K1-4*s.sqrt(6)*s.pi**2/45)==0
x=s.Matrix([1,0,0,0,0,0]);v=s.ones(6,1);vp=s.Matrix([1,-1,1,1,1,1])
def line_d2(a):return (x.T*x)[0]-(x.T*a)[0]**2/(a.T*a)[0]
assert line_d2(v)==line_d2(vp)==s.Rational(5,6)
checks.append({'configuration':7,'S1cap_coefficient':str(K1),
 'valid_delta':'0<delta<=1/sqrt(2)','support_linear_functional_norm_squared':2,
 'nonunique_nearest_line_point':[1,0,0,0,0,0],'two_line_distances_squared':'5/6','pass':True})
e,t=s.symbols('e t',positive=True)
R=1/(3*e**3)-1024/e+32768*(s.pi/2-s.atan(32*e))
assert s.simplify(s.diff(R,e)+e**(-4)/(1+1024*e**2))==0
assert s.limit(R,e,s.oo)==0
# Symbolic crossing relation, with positive Cf,V. No old scalar budget is evaluated.
Cf,V=s.symbols('Cf V',positive=True)
A=Cf**s.Rational(1,6)*s.pi**(-s.Rational(5,6))
eps=(Cf*s.pi*V)**s.Rational(1,5)
assert s.simplify((A*eps**s.Rational(1,6))/(Cf*V**s.Rational(5,6)*eps**(-4)))==1
assert s.integrate(t**s.Rational(1,6),(t,0,e))==s.Rational(6,7)*e**s.Rational(7,6)
checks.append({'configuration':8,'R_derivative_identity':True,'R_limit_infinity':0,
 'crossing_identity':True,'short_interval_denominator_direction':'>=',
 'Umin_monotone_in_V':'pointwise monotonicity, not monotonicity of L',
 'old_budget_reevaluations':0,'pass':True})
source=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-literature/cube-dim-6-bounds/upper/round034/positive_face_tube_budget/results.json')
d=json.loads(source.read_text());row=d['rows'][0];assert row['delta']=='1/16'
a=row['best_two_envelope_method_lower_bound'];b=row['method_lower_over_target']
lo,hi=Fraction(a['lower_exact']),Fraction(a['upper_exact'])
rlo,rhi=Fraction(b['lower_exact']),Fraction(b['upper_exact'])
safe=Fraction(1059393,200000)
assert 0<lo<=hi and 0<rlo<=rhi
assert 11250*lo>safe and rlo>safe>1
assert max(11250*lo,rlo)<=min(11250*hi,rhi)
checks.append({'configuration':9,'source_saved_delta':'1/16',
 'source_saved_ratio_lower_exact':str(rlo),'safe_lower_exact':str(safe),
 'strict_comparison_margin_exact':str(rlo-safe),
 'scaled_saved_L_lower_exceeds_safe':True,'saved_intervals_consistent':True,
 'independent_atan_or_budget_recomputation':False,'pass':True})
def save(a):return {'ball':str(a),'lower_exact':str(a.lower().fmpq()),'upper_exact':str(a.upper().fmpq())}
cn=arb.pi()**2/2*(10+10*arb(2).sqrt()+2*arb(5).sqrt())
k1=4*arb(6).sqrt()*arb.pi()**2/45
out={'status':'PASS_INDEPENDENT_GEOMETRY_AND_ANALYTIC_DIRECTION_REVIEW',
 'pid':os.getpid(),'started_utc':datetime.datetime.fromtimestamp(start,datetime.timezone.utc).isoformat(),
 'python':sys.version,'sympy':s.__version__,'python_flint':flint.__version__,'precision_bits':ctx.prec,
 'fixed_configurations':len(checks),'checks':checks,'Cnew':save(cn),'K1':save(k1),
 'source_hashes_verified':True,'source_sha256':sha(Path(__file__)),
 'registration_sha256':sha(HERE/'registration.md'),'input_manifest_sha256':sha(HERE/'input_manifest.json'),
 'kernel_calls':0,'native_calls':0,'matlab_calls':0,'MC_calls':0,'grid_calls':0,
 'old_4960_enumerations':0,'old_delta_budget_reevaluations':0,
 'new_I_certificate':False,'new_C6_bound':False,'elapsed_seconds':time.time()-start}
assert out['fixed_configurations']==9
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'pid':out['pid'],'configurations':len(checks),'elapsed_seconds':out['elapsed_seconds'],'Cnew':str(cn),'K1':str(k1),'saved_ratio_safe_lower':str(safe)},indent=2))
