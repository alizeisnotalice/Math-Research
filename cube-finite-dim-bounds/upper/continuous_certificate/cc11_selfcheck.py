from pathlib import Path
import runpy,json,time
from flint import arb,ctx
root=Path(__file__).resolve().parent;d=runpy.run_path(str(root/'cc11_diagnostic.py'));r=runpy.run_path(str(root/'cc11_regularized.py'));ctx.prec=160
point_checks=[]
for y,z in [('901/1000','403/1000'),('701/1000','203/1000'),('301/1000','103/1000'),('803/1000','337/1000')]:
 Y=d['F'](y);Z=d['F'](z);v,status=d['evaluate']((Y,Y,Z,Z));assert v is not None
 e=d['ns']['green_expression']((1,Y,Z),d['F'](1,3));expected=d['aq'](e['constant'])
 for a,c in e['log_coefficients'].items():expected+=d['aq'](c)*d['aq'](a).log()
 expected/=8*arb.pi()**3
 assert v.overlaps(expected),(y,z,v,expected)
 point_checks.append({'y':y,'z':z,'grouped':str(v),'exact_point':str(expected),'overlap':True})
box_checks=[]
for h in ['1/4096','1/8192','1/16384']:
 v,c=r['integral_box']('1','2272387/5000000',h,tol='1e-10')
 box_checks.append({'center':['1','2272387/5000000'],'halfwidth':h,'epsilon':'1/50','lower':str(v.real.lower()),'upper':str(v.real.upper()),'above_target':bool(v.real>r['q']('13/2500')),'integrand_calls':c})
assert box_checks[-1]['above_target']
ratios=[]
for cone,ks,w in r['ns']['_prepared'](r['F'](1,3)):
 for ray,k in zip(cone,ks):
  norm=abs(ray[1])+abs(ray[2])
  if norm:ratios.append(k/norm)
assert min(ratios)==3
rho=arb(3);n=20;h=arb(1)/32;s=h*(rho-1/rho)/2;delta=arb(1)/50-s/3
MC=15/(1088*arb.pi()**3*delta**2)
cheb={'patches':136,'nodes_per_patch':441,'total_nodes':59976,'h':'1/32','rho':3,'degree':20,'delta':str(delta),'M':str(MC),'uniform_interpolation_error':str((2*n+2)*4*MC*rho**(-n)/(rho-1)),'node_error_amplification':1681}
M=arb(15)/(1088*arb.pi()**3*(arb(7)/1600)**2);q=arb(1)/4
bounds={str(N):str(M*q**(N+1)*((N+2)-(N+1)*q)/(1-q)**2) for N in (12,14,16,18)}
(root/'selfcheck.json').write_text(json.dumps({'precision_bits':ctx.prec,'exact_log_point_overlap':point_checks,'regularized_local_boxes':box_checks,'cauchy_total_degree_remainder':{'epsilon':'1/50','polydisc_radius':'1/32','box_halfwidth':'1/128','M':str(M),'bounds':bounds},'improved_tube_minimum_ray_ratio':str(min(ratios)),'recommended_chebyshev_batch':cheb,'global_certificate':False},indent=2)+'\n')
ddchecks=[]
for Y,Z in [('1','2272387/5000000'),('1','0'),('1','1'),('1/2','1/2'),('0','0')]:
 v=r['point_divided_difference'](Y,Z);w,c=r['integral_box'](Y,Z,tol='1e-10');assert v.overlaps(w.real)
 ddchecks.append({'y':Y,'z':Z,'DD':str(v),'quadrature_overlap':True})
start=time.perf_counter()
for _ in range(100):r['point_divided_difference'](r['F'](901,1000),r['F'](403,1000))
elapsed=time.perf_counter()-start
(root/'dd_benchmark.json').write_text(json.dumps({'scope':'benchmark repeats one fixed generic point, no new minimum scan','precision_bits':ctx.prec,'checks':ddchecks,'calls':100,'wall_seconds':elapsed,'estimated_59976_node_seconds':elapsed*599.76,'warning':'extrapolation excludes nonrational Chebyshev node enclosure, polynomial validation, and confluent fallback cost'},indent=2)+'\n')
print(json.dumps({'point_checks_passed':len(point_checks),'local_box_checks':box_checks,'remainder_bounds':bounds},indent=2))
