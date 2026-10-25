#!/usr/bin/env python3
"""One exact independent structural replay; no kernel or numerical triple."""
from pathlib import Path
from itertools import product
import hashlib,json,time
import sympy as s
START=time.monotonic()
HERE=Path(__file__).resolve().parent
A=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-curvature-astra/cube-dim-6-bounds/upper/round010/short_chamber_curvature')
V=[s.Matrix([1,*a]) for a in product((-1,1),repeat=5)]
AA=sum((v*v.T for v in V),s.zeros(6,6));assert AA==32*s.eye(6)
v=s.ones(6,1);w=s.Matrix([1,-1,-1,-1,-1,-1])
columns=[v,w]
for j in range(1,5):
    col=v.copy();col[j]=-1;columns.append(col)
B=s.Matrix.hstack(*columns);det=B.det();inv=B.inv()
c=s.Matrix([1]+[s.Rational(1,2)]*5)
z=s.Matrix([0,51,22,2,-24,-51])/200
assert det==-32 and inv*c==s.Matrix([s.Rational(3,4),s.Rational(1,4),0,0,0,0])
rest=[u for u in V if all(u!=col for col in columns)]
assert len(rest)==26
R=inv*s.Matrix.hstack(*rest)
rows=[sum(abs(R[i,j]) for j in range(26)) for i in range(6)]
assert rows==[49,15,15,15,15,15]
C=max(rows);t0=min(s.Integer(1),1/C)
coefficient=1/(s.Integer(2)**26*25*abs(det)*(1+C*C)**4)
assert t0==s.Rational(1,49)
assert coefficient==s.Rational(1,1787153625614625485619200)

# All path parameters 0<sigma<=1, not just one point.
sigma=s.symbols('sigma',real=True);path=c+sigma*z
def positive_on_unit_open_left(expr):
    expr=s.expand(expr)
    assert s.degree(expr,sigma)<=1
    assert expr.subs(sigma,0)>=0 and expr.subs(sigma,1)>0
    return str(expr)
path_gaps=[positive_on_unit_open_left(path[i]-path[i+1]) for i in range(5)]
positive_on_unit_open_left(path[5])
assert sum(z[1:])==0
partials=[sum(z[1:i+1]) for i in range(1,5)]
assert partials==[s.Rational(a,200) for a in (51,73,75,51)]
gg=s.symbols('g2:7')
assert s.expand(sum(z[i+1]*gg[i] for i in range(5))-
 sum(partials[i]*(gg[i]-gg[i+1]) for i in range(4)))==0
wall_checks=[]
for i in range(1,5):
    wall=list(path);wall[i]=wall[i+1]=(path[i]+path[i+1])/2
    gaps=[positive_on_unit_open_left(wall[j]-wall[j+1]) for j in range(5) if j!=i]
    positive_on_unit_open_left(wall[-1])
    assert s.simplify(wall[i]-wall[i+1])==0
    wall_checks.append({'pair':[i+1,i+2],'all_other_gaps_positive_for_every_0_lt_sigma_le_1':True,
                        'nonzero_gap_formulas':gaps})
patterns=list(product((-1,1),repeat=3))
opposite_pairs={p if p[0]==1 else tuple(-a for a in p) for p in patterns}
assert len(opposite_pairs)==4

# Derivative bounds are for the real density, not differentiated heat envelopes.
u,t=s.symbols('u t',real=True,positive=True)
f=t/(s.pi*(t*t+u*u))
assert s.simplify(s.diff(f,u)/f+2*u/(t*t+u*u))==0
assert s.simplify(s.diff(f,u,2)/f-(6*u*u-2*t*t)/(t*t+u*u)**2)==0
rho=s.symbols('rho',nonnegative=True)
assert s.expand(2*(1+rho)**2-(6*rho-2)-(2*(rho-s.Rational(1,2))**2+s.Rational(7,2)))==0
assert s.expand(2*(1+rho)**2+(6*rho-2)-(2*rho*rho+10*rho))==0
assert s.expand((1+rho)**2-4*rho-(rho-1)**2)==0

saved=json.loads((A/'structural_certificate.json').read_text())
assert saved['basis_determinant']==str(det) and saved['C']==str(C) and saved['t0']==str(t0)
assert saved['positive_p_t_c_lower_coefficient_times_pi6']==str(coefficient)
paths=[A/name for name in ('theorem.md','report.md','verify_structure.py','structural_certificate.json','registration.md','input_manifest.json')]
out={'status':'PASS_ANALYTIC_NONCONVEXITY_STRUCTURE; EXISTENCE_ONLY',
 'structural_replays':1,'kernel_evaluations':0,'MC_samples':0,'matlab_calls':0,
 'basis_determinant':str(det),'B_inverse_c':[str(a) for a in inv*c],
 'row_l1_bounds':[str(a) for a in rows],'C':str(C),'t0':str(t0),
 'density_t_minus2_coefficient_times_pi6':str(coefficient),
 'path_all_parameters_verified':True,'path_gap_formulas':path_gaps,
 'positive_partial_sums':[str(a) for a in partials],'summation_by_parts_verified':True,
 'average_walls':wall_checks,'S3_absolute_value_count_upper':4,
 'first_Cauchy_relative_derivative_bound':'1/t','second_Cauchy_relative_derivative_bound':'2/t^2',
 'analytic_theorem_scope':'For every fixed T>0, there exists an interior strict ordered face point and one adjacent free-coordinate difference direction with negative second derivative of K_T.',
 'explicit_negative_point_or_triple_certified':False,'whole_cube_minimum_certified':False,'new_C6_bound':False,
 'a_e_to_pointwise_clarification':'Use the open eta/4 neighborhood of a compact E at distance eta from S3 before continuous-density extension; the same delta=eta/(64sqrt6) works.',
 'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'elapsed_seconds':time.monotonic()-START,'exit_code':0}
(HERE/'check_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
