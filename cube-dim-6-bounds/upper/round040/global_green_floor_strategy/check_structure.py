#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import json,hashlib,os,sys,time
import sympy as s
import flint
from flint import arb,ctx
ctx.prec=384
r=Path(__file__).resolve().parent;start=time.time()
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
inp=json.loads((r/'input_manifest.json').read_text())
for p,h in inp['source_sha256'].items():assert sha(p)==h,p
assert not (r/'results.json').exists()
V=[(1,)+q for q in product((1,-1),repeat=5)]
v0=(1,)*6
flip=lambda j:tuple(-1 if i==j else 1 for i in range(6))
cols=[v0,flip(4),flip(5),flip(1),flip(2),flip(3)]
A=s.Matrix.hstack(*map(s.Matrix,cols));Ai=A.inv()
assert abs(A.det())==32 and Ai*A==s.eye(6)
checks=[{'configuration':1,'basis_columns':cols,'determinant':str(A.det()),'inverse':str(Ai),'pass':True}]
rest=[v for v in V if v not in cols];assert len(rest)==26
R=s.Matrix.hstack(*map(s.Matrix,rest));M=Ai*R
bounds=[sum(abs(z) for z in M.row(i)) for i in range(6)]
assert bounds==[49,15,15,15,15,15]
checks.append({'configuration':2,'body_directions':32,'remaining_directions':26,'exact_row_absolute_sums':list(map(str,bounds)),'pass':True})
b,c=s.symbols('b c',real=True)
x=s.Matrix([1,1,1,1,b,c]);a=Ai*x
assert a==s.Matrix([(b+c)/2,(1-b)/2,(1-c)/2,0,0,0])
A0=s.Rational(1,100)**2+s.Rational(149,100)**2
A1=s.Rational(1,100)**2+s.Rational(13,20)**2
assert A0==s.Rational(11101,5000) and A1==s.Rational(2113,5000)
C=arb(5000)**3/(arb(2)**31*arb.pi()**6*arb(226)**3*11101*arb(2113)**2)
assert C>arb(1)/10**20
checks.append({'configuration':3,'S3_parameters':str(a),'time_max':'1/100','A0':str(A0),'A1':str(A1),'density_lower_exact':'5000^3/(2^31*pi^6*226^3*11101*2113^2)','density_lower_ball':str(C),'density_lower_endpoint':str(C.lower().fmpq()),'density_upper_endpoint':str(C.upper().fmpq()),'strictly_above':'1/10^20','pass':True})
xq=[F(1)]*4+[F(1,2),F(0)]
def var(q):
 u=sum(q)/len(q);return sum((z-u)**2 for z in q)
d2=[var(xq[:k])+var(xq[k:]) for k in range(1,6)]
assert min(d2)==F(1,8) and min(d2)>F(1,125000)**2
checks.append({'configuration':4,'point':list(map(str,xq)),'five_squared_distances':list(map(str,d2)),'nearest_squared_distance':'1/8','outside_paid_tube':True,'pass':True})
source=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-curvature-astra/cube-dim-6-bounds/upper/round038/green_floor_diagnostic/results.json')
d=json.loads(source.read_text());Gup=F(d['points'][2]['check']['G_interval'][1])
tube=(F(125000)-40000)/5000000000
assert tube==F(17,1000000)>Gup>F(1,65000)
checks.append({'configuration':5,'delta':'1/125000','A009_paid_G_floor':str(tube),'saved_point3_G_upper_exact':str(Gup),'strictly_exceeds_saved_G_upper_and_target':True,'pass':True})
w=[F(1)]*5+[F(0)]
for j,p in enumerate(d['points'],start=6):
 q=[F(n,p['denominator']) for n in p['numerators']]
 h=[wi-qi for wi,qi in zip(w,q)]
 gs=[[F(z) for z in a] for a in p['check']['gradient_intervals']]
 lower=sum((a if t>=0 else b)*t for (a,b),t in zip(gs,h))
 upper=sum((b if t>=0 else a)*t for (a,b),t in zip(gs,h))
 assert lower<=upper<0
 checks.append({'configuration':j,'A038_point':p['index'],'feasible_vertex':list(map(str,w)),'direction':list(map(str,h)),'directional_derivative_lower_exact':str(lower),'directional_derivative_upper_exact':str(upper),'strictly_negative':True,'proved_neighborhood_radius':None,'pass':True})
out={'status':'PASS_NEW_S3_MOMENT_OBSTRUCTION_AND_PAID_COMPACT_REDUCTION_GLOBAL_FLOOR_UNRESOLVED','pid':os.getpid(),'fixed_configurations':len(checks),'checks':checks,'source_sha256':sha(Path(__file__)),'input_manifest_sha256':sha(r/'input_manifest.json'),'registration_sha256':sha(r/'registration.md'),'python':sys.version,'sympy':s.__version__,'python_flint':flint.__version__,'precision_bits':ctx.prec,'elapsed_seconds':time.time()-start,'kernel_G_H_p_calls':0,'native_calls':0,'matlab_calls':0,'MC_grid_optimizer_calls':0,'new_global_g_certificate':False,'new_C6_bound':False}
assert len(checks)==8
(r/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'pid':out['pid'],'configurations':len(checks),'elapsed_seconds':out['elapsed_seconds'],'S3_density_lower':str(C),'paid_tube_floor':str(tube),'saved_gradient_dots_strictly_negative':3},indent=2))
