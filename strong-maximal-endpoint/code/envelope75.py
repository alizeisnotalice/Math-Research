import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import sys,json,math,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from profile74 import MU,law

def run(j):
 n=j['N'];K=j['K'];L=math.log(n);rho=math.sqrt(L/n);cr=MU-math.sqrt(MU)*rho/2;cm=cr-K*rho;cp=cr+K*rho;fm,_,_,_=law(n,cm);fp,_,_,_=law(n,cp);zm=n*math.log1p(math.exp(-cm));zp=n*math.log1p(math.exp(-cp));rows=[]
 def lm(s):return fm(s)+s+zm
 for t in [-7.,-5.,-3.,-1.,0.,1.,3.,5.,7.]:
  s=n*MU+t*math.sqrt(n*L);hp=fp(s)+s+zp;hm=lm(s);ss=s-n*(cp-cm);margin=None if ss<=0 else hp-lm(ss)
  assert hm>=hp-1e-9 and (margin is None or margin>=-1e-9)
  rows.append(dict(S=s,radial_node=t,log_envelope_ratio=hm-hp,shifted_S=ss,log_shift_comparison_margin=margin))
 return dict(j,c_lower=cm,c_upper=cp,maximum_sample_log_envelope=max(r['log_envelope_ratio'] for r in rows),normalized_max_sample_log_envelope=max(r['log_envelope_ratio'] for r in rows)/math.sqrt(n*L),shifted_lower_positive=rows[0]['shifted_S']>0,rows=rows)
if __name__=='__main__':
 out=Path(sys.argv[1]);out.mkdir(exist_ok=True,parents=True);jobs=[dict(N=n,K=k) for n in [256,1024,4096,16384,65536,262144,1048576] for k in [.5,1.,2.]];start=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(run,jobs))
 (out/'envelope75.json').write_text(json.dumps(dict(seconds=time.monotonic()-start,workers=12,rows=rows,scope='finite radial samples, not certified global supremum; uniform envelope follows from analytic normalizer comparison'),indent=2)+'\n');print('PASS',len(rows),'coarse vector-envelope samples')
