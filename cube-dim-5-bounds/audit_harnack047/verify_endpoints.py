#!/usr/bin/env python3
"""Read-only audit: exact algebra and saved endpoints, never evaluate Green."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import argparse
import hashlib
import json
import os
import signal
import sys
import time
import sympy as sp

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument('--output-dir', type=Path, default=HERE)
args = parser.parse_args()
OUT = args.output_dir.resolve()
OUT.mkdir(parents=True, exist_ok=True)
start = time.monotonic()

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def checkpoint(stage):
    d = {'pid': os.getpid(), 'stage': stage, 'elapsed_seconds': time.monotonic()-start,
         'deadline_seconds': 120,
         'recoverable_command': f'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 {sys.executable} -B {__file__} --output-dir {OUT}'}
    (OUT/'checkpoint.json').write_text(json.dumps(d,indent=2)+'\n')
def timeout(_s,_f):
    checkpoint('timeout_no_certificate')
    raise TimeoutError('registered 120-second deadline')
signal.signal(signal.SIGALRM,timeout)
signal.alarm(120)
checks = {}
def identity(name,lhs,rhs):
    assert sp.simplify(lhs-rhs)==0, name
    checks[name]={'lhs':str(lhs),'rhs':str(rhs),'residual':'0'}

checkpoint('hashes_and_saved_state')
reg=json.loads((HERE/'registration.json').read_text())
src=Path(reg['source_root'])
for name,digest in reg['sources'].items(): assert sha(src/name)==digest, name
d=json.loads((src/'results.json').read_text())
oldreg=json.loads((src/'registration.json').read_text())
assert sha(src/'registration.json')==d['source_sha256']['registration.json']
saved_checkpoint=Path(d['checkpoint']['path'])
assert sha(saved_checkpoint)==d['checkpoint']['sha256']
state=json.loads(saved_checkpoint.read_text())
assert state['next_simplex']==3504 and state['status']=='complete'
assert state['fingerprint']==d['checkpoint']['fingerprint']
fingerprint=hashlib.sha256(json.dumps({'source':reg['sources']['certify.py'],
    'reg':sha(src/'registration.json'),'inputs':oldreg['read_only_inputs']},sort_keys=True).encode()).hexdigest()
assert fingerprint==state['fingerprint']
assert sum(map(len,state['stencil']))==704 and len(state['stencil'])==5
checks['saved_state']={'processed_simplices':3504,'stencil_quantities':5,'stencil_entries':704,
                       'fingerprint_verified':True,'source_producer_executed':False}

checkpoint('scalar_and_matrix_algebra')
r=sp.symbols('r',nonnegative=True)
f=2*(r-1)/(1+r)**2
identity('scalar_lower_gap',f+2,2*r*(r+3)/(1+r)**2)
identity('scalar_upper_gap',sp.Rational(1,4)-f,(r-3)**2/(4*(1+r)**2))
t,w=sp.symbols('t w',real=True)
identity('score_bound_square',(t*t+w*w)**2-4*t*t*w*w,(t*t-w*w)**2)
identity('Hessian_upper_constant',sp.Rational(1,4)+5,sp.Rational(21,4))
identity('paraboloid_lower_constant',sp.Rational(-2,2),-1)
identity('paraboloid_upper_constant',sp.Rational(21,4)/2,sp.Rational(21,8))
N=sp.Matrix.hstack(*(sp.Matrix(v) for v in d['N_columns']))
assert N.det()==-48
Ninv=N.inv()
assert [[str(x) for x in row] for row in Ninv.tolist()]==d['N_inverse']
Q=(Ninv.T*Ninv)[1:5,1:5]
assert Q==sp.eye(4)/4-sp.ones(4)/36
Qi=4*sp.eye(4)+sp.Rational(4,5)*sp.ones(4)
assert Q*Qi==sp.eye(4)
assert Q.eigenvals()=={sp.Rational(1,4):3,sp.Rational(5,36):1}
checks['fixed_basis_and_Q']={'det_N':'-48','Q':'I_4/4-J_4/36','Q_inverse':'4I_4+(4/5)J_4',
    'eigenvalues':{'1/4':3,'5/36':1},'positive_definite':True}
v,k,l,A,B,P0,P1,delta,den=sp.symbols('v k l A B P0 P1 delta den')
identity('initial_node_physical_derivative',sp.diff(sp.log(A*A+B*B*(k+l*v)**2)/2,v).subs(v,0),
    B*B*k*l/(A*A+B*B*k*k))
identity('first_jet_quotient_derivative',sp.diff((P0+P1*v)/(den+delta*v),v).subs(v,0),
    (P1-delta*(P0/den))/den)
D=sp.symbols('D',positive=True)
identity('cone_time_scale',D**5/(144*D),D**4/144)
identity('scaled_epsilon',sp.Rational(144*104729,64),sp.Rational(942561,4))

checkpoint('saved_intervals_and_cube')
eps=F(d['epsilon']); target=F(d['paraboloid_regions']['target'])
Hlo,Hhi=map(F,d['H_center_directed_interval'])
assert 0<target<Hlo<Hhi
gi=[tuple(map(F,row)) for row in d['tail_loggradient_directed_intervals']]
assert len(gi)==4 and all(lo<=hi for lo,hi in gi)
assert len(d['four_tail_gradient_H_balls'])==4
region=d['paraboloid_regions']; Delta=F(region['Delta_lower'])
assert Delta>0
h=F(region['safe_centered_cube_halfwidth']); h0=F(d['fixed_closed_box']['tail_halfwidth'])
def quad(sign):return sum(x*x for x in sign)/F(4)-sum(sign)**2/F(36)
def dot_lower(sign):return sum((lo if s==1 else -hi) for s,(lo,hi) in zip(sign,gi))
signs=list(product((-1,1),repeat=4))
assert max(map(quad,signs))==1
saved_rows={tuple(row['sign']):row for row in region['maximal_h_root_rows']}
assert set(saved_rows)==set(signs)
saved_residuals=[]; safe_interval_residuals=[]; fixed_residuals=[]
root_lows=[];root_ups=[]
for sign in signs:
    row=saved_rows[sign]; c=quad(sign); b=F(row['log_gradient_dot_lower'])
    assert c==F(row['Q_sign']) and c>0
    lo,hi=F(row['lower']),F(row['upper'])
    poly=lambda x:Delta+b*x-c*x*x/(eps*eps)
    assert 0<lo<=hi and poly(lo)>=0 and poly(hi)<=0
    root_lows.append(lo);root_ups.append(hi)
    saved_residuals.append(poly(h))
    safe_interval_residuals.append(Delta+dot_lower(sign)*h-c*h*h/(eps*eps))
    fixed_residuals.append(Delta+dot_lower(sign)*h0-c*h0*h0/(eps*eps))
assert min(saved_residuals)==F(region['minimum_cube_log_residual'])>0
assert min(safe_interval_residuals)>0 and min(fixed_residuals)>0
assert [min(root_lows),min(root_ups)]==list(map(F,region['maximal_safe_h_for_conservative_corner_test_interval']))
assert h==F((min(root_lows).numerator*(1<<24))//min(root_lows).denominator,1<<24)
point=list(map(F,d['point'])); assert point[0]==1
margin=min(min(x,1-x) for x in point[1:])
assert margin==F(614,104729)
assert F(region['safe_halfwidth_clipped_to_unit_face'])==min(h,margin)
fixed_lower=F(d['fixed_closed_box']['H_uniform_lower_bound_dyadic'])
assert fixed_lower>F(11,100000)>target
checks['cube_and_fixed_box']={'sign_corners':16,'Q_sign_max':'1','safe_halfwidth':str(h),
    'face_halfwidth':str(min(h,margin)),'saved_exact_min_residual':str(min(saved_residuals)),
    'min_residual_using_independent_saved_gradient_rectangle':str(min(safe_interval_residuals)),
    'fixed_box_primary_margin_from_paraboloid_positive':True,
    'saved_fixed_box_lower_exceeds_both_targets':True,'all_saved_root_brackets_verified':True}

checkpoint('ellipsoid')
ellipse=region['centered_ellipsoid']; Bup=F(ellipse['dual_loggradient_norm_upper']);rho=F(ellipse['rho'])
gradient_vertices=list(product(*gi))
dual_sq=lambda g:4*sum(x*x for x in g)+F(4,5)*sum(g)**2
max_dual=max(map(dual_sq,gradient_vertices))
assert Bup>0 and Bup*Bup>=max_dual
ell_res=Delta-Bup*rho-rho*rho/(eps*eps)
assert ell_res==F(ellipse['log_residual_lower'])>0
assert all(F(ellipse['Q'][i][j])==F(str(Q[i,j])) for i in range(4) for j in range(4))
checks['ellipsoid']={'gradient_rectangle_vertices':16,'B_upper_covers_entire_gradient_rectangle':True,
    'B_squared_minus_rectangle_max':str(Bup*Bup-max_dual),'rho':str(rho),
    'exact_log_residual':str(ell_res),'residual_positive':True}

signal.alarm(0)
results={'status':'pass_043_theory_source_review_and_saved_endpoint_checks',
    'source_sha256':reg['sources'],'checks':checks,'transcendental_calls':0,
    'Green_or_H_calls':0,'source_producer_executed':False,'new_global_g_claim':False,
    'numerical_scope':'four tail gradients on x1=1 only; no fifth gradient',
    'dependency':'Saved Delta and H/gradient enclosures are accepted after source review; no new log/exp evaluation.'}
(OUT/'results.json').write_text(json.dumps(results,indent=2)+'\n')
checkpoint('complete')
receipt={'pid':os.getpid(),'elapsed_seconds':time.monotonic()-start,'exit_code':0,
         'python':sys.version.split()[0],'sympy':sp.__version__,'record_groups':len(checks),'threads':1}
(OUT/'run_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
