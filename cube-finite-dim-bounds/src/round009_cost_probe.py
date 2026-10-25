"""One finite n3 cost probe; numerical evidence only, no global minimum proof."""
from pathlib import Path
import importlib.util,hashlib,json,sys,time
import numpy as np
from scipy.optimize import minimize
modpath=Path(sys.argv[1]).resolve(); outpath=Path(sys.argv[2]).resolve()
spec=importlib.util.spec_from_file_location('cb9_frame_green',modpath);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
start=time.time(); records=[]
def f(z):
 y,z=map(float,z)
 return float(m.green_face(y,z))
# One deterministic closed triangular-face discovery, not a certified minimum.
for y in np.linspace(0,1,9):
 for z in np.linspace(0,y,9):
  if y==1 and z==1:continue
  records.append((float(y),float(z),f((y,z))))
best=min(records,key=lambda z:z[2]);opt=minimize(f,best[:2],method='SLSQP',bounds=[(0,1),(0,1)],constraints=[{'type':'ineq','fun':lambda a:a[0]-a[1]},{'type':'ineq','fun':lambda a:1.999999-a[0]-a[1]}],options={'maxiter':120,'ftol':1e-13})
g=min(best[2],float(opt.fun));report={'status':'numerical only: sampled minimum is an upper bound on true g, so these costs are not certified upper bounds','source':str(modpath),'source_sha256':hashlib.sha256(modpath.read_bytes()).hexdigest(),'grid_best':best,'optimizer':{'x':opt.x.tolist(),'fun':float(opt.fun),'success':bool(opt.success),'message':str(opt.message),'nfev':int(opt.nfev)},'g_used':g,'axis_exact_numeric':3/(8*np.pi**3),'rows':[]}
def save():
 report['elapsed_seconds']=time.time()-start;outpath.write_text(json.dumps(report,indent=2)+'\n')
save()
for order in (32,64):
 nodes,weights=np.polynomial.legendre.leggauss(order);nodes=(nodes+1)/2;weights=weights/2
 total=0.
 for y,wy in zip(nodes,weights):
  vals=np.array([f((y,y*v)) for v in nodes]);assert np.all(vals>0) and np.all(np.isfinite(vals))
  total+=float(2*y*wy*np.dot(weights,(vals/g)**1.5))
 report['rows'].append({'order':order,'face_average_estimate':total,'cost_estimate':3*np.sqrt(3)*total,'target':8,'below_target_numerically':bool(3*np.sqrt(3)*total<8)});save();print(json.dumps(report['rows'][-1]),flush=True)
assert hashlib.sha256(modpath.read_bytes()).hexdigest()==report['source_sha256'],'source changed during run; invalidate/review result'
