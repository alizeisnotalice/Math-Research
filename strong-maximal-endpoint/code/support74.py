import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import math,json,sys,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from profile74 import MU,law

def run(j):
 n=j['N'];L=math.log(n);rho=math.sqrt(L/n);D=1.5 if j['D_type']=='constant' else L**(-.25);cr=MU-math.sqrt(MU)/2*rho;ca=cr+D*rho;fa,_,_,ma=law(n,ca);fr,_,_,mr=law(n,cr);lv=n*(math.log1p(math.exp(-ca))-math.log1p(math.exp(-cr)));p=1/(1+math.exp(MU));rows=[]
 for t in [-7.,-3.,-1.,0.,3.,7.]:
  s=n*MU+t*math.sqrt(n*L)
  if s<=0:continue
  ld=fa(s)-fr(s);rows.append(dict(S=s,radial_node=t,log_density_ratio=ld,log_H_ratio=lv+ld))
 return dict(j,c_reference=cr,c_actual=ca,Delta_c=ca-cr,log_volume_ratio=lv,leading_log_volume_ratio=-p*n*(ca-cr),normalized_log_volume_ratio=lv/(n*(ca-cr)),max_abs_log_density_ratio=max(abs(r['log_density_ratio']) for r in rows),rows=rows)
if __name__=='__main__':
 out=Path(sys.argv[1]);out.mkdir(exist_ok=True,parents=True);jobs=[dict(N=n,D_type=t) for n in [64,256,1024,4096,16384,65536,262144,1048576] for t in ['constant','vanishing']];start=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(run,jobs))
 (out/'support74.json').write_text(json.dumps(dict(seconds=time.monotonic()-start,workers=12,rows=rows),indent=2)+'\n');print('PASS',len(rows),'nested support comparisons')
