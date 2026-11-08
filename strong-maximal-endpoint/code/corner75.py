import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import math,json,sys,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from scipy.optimize import brentq
MU=brentq(lambda x:x-1-math.exp(-x),1,2,xtol=1e-15)
ZETA=math.log1p(math.log(2))/math.log(2)
def run(n):
 L=math.log(n);rho=math.sqrt(L/n);cr=MU-math.sqrt(MU)*rho/2;ca=cr+L**.25/math.sqrt(n);ra=math.exp(-ca);rr=math.exp(-cr);p=1+ZETA*(n-1);bound=n*p*math.log1p(1/(p-1));eps=1/n**2
 return dict(N=n,p=p,log_Hardy_bound=bound,log_bound_minus_N=bound-n,limiting_offset=1/(2*ZETA),c_reference=cr,c_actual=ca,log_reanchor_Jacobian=-n*math.log1p(rr),log_weight_at_receiver_1=n*(math.log1p(ra)-math.log1p(rr)),log_weighted_density_center=n*math.log1p(ra),log_weighted_density_lower=n*math.log1p(ra/(1+eps)),log_weighted_density_upper=n*math.log1p(ra/(1-eps)))
if __name__=='__main__':
 out=Path(sys.argv[1]);out.mkdir(exist_ok=True,parents=True);ns=[2,4,8,16,32,64,128,256,1024,4096,16384,65536,262144,1048576];start=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(run,ns))
 (out/'corner75.json').write_text(json.dumps(dict(mu=MU,zeta=ZETA,workers=12,seconds=time.monotonic()-start,rows=rows,scope='fixed-corner global Hardy upper and qualified-patch bare geometry; not the general endpoint upper'),indent=2)+'\n');print('PASS',len(rows),'fixed-corner and affine geometry constants')
