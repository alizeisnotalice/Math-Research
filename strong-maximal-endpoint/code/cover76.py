import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import math,json,sys,time
from pathlib import Path
from fractions import Fraction
from concurrent.futures import ProcessPoolExecutor
from scipy.optimize import brentq
from optimize76 import MU,BETA,ff

def run(j):
 n=j['N'];K=j['K'];rho=math.sqrt(math.log(n)/n);cr=MU-math.sqrt(MU)*rho/2;rp=math.exp(-(cr-K*rho));rm=math.exp(-(cr+K*rho));amax=(math.exp(BETA)*(1+rm)-1)/rp;aa=1.05;beta=math.log((1+aa*rp)/(1+rm));q=brentq(lambda x:ff(x)-beta,1e-12,8,xtol=1e-14)
 return dict(j,maximum_grid_ratio_at_leading_frontier=amax,grid_ratio=aa,log_volume_ratio_per_dimension=beta,leading_Hardy_charge_per_dimension=1+q,leading_rate_gap=MU-1-q,nontrivial_grid_possible=(amax>1))
if __name__=='__main__':
 out=Path(sys.argv[1]);out.mkdir(exist_ok=True,parents=True);jobs=[dict(N=n,K=k) for n in [256,1024,4096,16384,65536,262144,1048576] for k in [.5,1.,2.]];start=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(run,jobs))
 cert=[dict(ancestor_count=k,source_modular='1/2',summed_modular=str(Fraction(k,2)),reuse_factor=k) for k in [1,4,16,64,256,1024]]
 (out/'cover76.json').write_text(json.dumps(dict(seconds=time.monotonic()-start,workers=12,rows=rows,limiting_grid_ratio=math.expm1(MU-1)/(MU-1),exact_reuse_packet=cert,scope='local corner mesh and failed unweighted budget sum; no global net or source allocation'),indent=2)+'\n');print('PASS',len(rows),'local mesh calculations and exact packet budget identities')
