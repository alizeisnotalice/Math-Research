"""One fixed-point necessary condition probe; no strict interval or global claim."""
from pathlib import Path
import importlib.util,json,math,hashlib,time
from scipy.integrate import quad
p=Path(__file__).resolve().parents[1];module=p/'upper/continuous_bridge/cb9_frame_green.py';sha=hashlib.sha256(module.read_bytes()).hexdigest();assert sha=='d23501eeaf417a736730a6e61398e27b9364c231d35ec6ee8bf573e357058bf5'
spec=importlib.util.spec_from_file_location('cb9',module);cb=importlib.util.module_from_spec(spec);spec.loader.exec_module(cb)
ge=(4*math.log(2)-1)/(8*math.pi**3);c=ge/3;lb=32/(9*math.pi);rows=[];start=time.time()
for tol in (1e-6,1e-8):
 def f(z):
  if z==0:return 0.
  s=z/(1-z)
  plus=cb.green((1+s,s,s));minus=cb.green((1-s,-s,-s))
  return (max(c-plus,0)+max(c-minus,0))/(z*z)
 total=0.;err=0.
 for a,b in zip([0,1/3,.5,.75],[1/3,.5,.75,1]):
  v,e=quad(f,a,b,epsabs=tol,epsrel=tol,limit=150);total+=v;err+=e
 mu=4/math.pi*total
 rows.append({'tolerance':tol,'mu_density_estimate':mu,'outer_quad_error_estimate_only':4/math.pi*err,'epsilon_max_estimate':mu/lb})
old=json.loads((p/'results/round009_cost_probe.json').read_text());cost=old['rows'][-1]['cost_estimate'];massb=32*math.pi/105;need=8*(ge-c)*(cost-8)/massb
out={'status':'numerical fixed-point necessary-condition probe only','module_sha256':sha,'x':[1,0,0],'g_edge':ge,'c':c,'Lb_axis_exact':'32/(9*pi)','rows':rows,'epsilon_needed_using_unproved_previous_cost':need,'scope':'no exact mu interval; no certified g/J; no new bound or excluded family','elapsed_seconds':time.time()-start}
(p/'results/round009_bump_necessary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
