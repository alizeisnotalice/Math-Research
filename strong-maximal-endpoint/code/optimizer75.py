import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import json,sys,math,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from scipy.optimize import brentq
MU=brentq(lambda x:x-1-math.exp(-x),1,2,xtol=1e-15)
def run(n):
 rho=math.sqrt(math.log(n)/n);cr=MU-math.sqrt(MU)*rho/2;rr=math.exp(-cr);ell0=1+rr;l0=-rr;H=2*rho;delta=rr/(1+rr)*rho/8;eta=B=.05;d=B*ell0*H/(ell0-2*H);assert rho<.1 and H<ell0/4
 def centered(x,c):
  ell=x*(1+math.exp(-c));a=x-ell
  assert abs(a-l0)<H and abs(x-1)<H
  val=eta*(l0+H-a)+B*(H*H-(a-l0)**2)/2+(eta-d)*(ell0-2*H)+(eta+B*ell0)*(x-1+H)+B*((x-1)**2-H*H)/2
  return val/ell
 rows=[]
 for h in [-1.,-.5,0.,.5,1.]:
  x=1+h*delta;root=brentq(lambda c:centered(x,c)-(eta+B*(-x*math.exp(-c)-l0)),cr-rho,cr+rho,xtol=1e-13)
  expected=-math.log(ell0/x-1);samples=[]
  for t in [-1.,-.5,0.,.5,1.]:
   c=cr+t*rho;ell=x*(1+math.exp(-c));value=centered(x,c);formula=eta+B*(x-1)-B*(ell-ell0)**2/(2*ell);samples.append(dict(c=c,centered_average=value,formula_error=value-formula))
  X=x*(1+math.exp(-root))/(1+rr);assert abs(root-expected)<1e-10 and abs(X-1)<1e-10
  rows.append(dict(x=x,x_offset=h,c_optimizer=root,c_formula=expected,reanchored_receiver=X,centered_maximum=centered(x,root),samples=samples))
 return dict(N=n,c_reference=cr,rho=rho,H=H,delta=delta,source_tensor_log_cap_upper=n*math.log1p((eta+B*(ell0+H))/n),receiver_cube_log_min_average=n*math.log1p((eta-B*delta)/n),maximum_formula_error=max(abs(s['formula_error']) for r in rows for s in r['samples']),rows=rows)
if __name__=='__main__':
 out=Path(sys.argv[1]);out.mkdir(exist_ok=True,parents=True);start=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(run,[1024,4096,16384,65536,262144,1048576]))
 (out/'optimizer75.json').write_text(json.dumps(dict(workers=12,seconds=time.monotonic()-start,rows=rows,scope='actual compact tensor input; unique maximizers only in the prescribed critical action window; no local efficiency qualification'),indent=2)+'\n');print('PASS 6 tensor sources and 30 maximizing policy roots')
