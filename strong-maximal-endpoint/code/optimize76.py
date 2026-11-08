import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import math,json,sys,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from scipy.optimize import brentq,minimize_scalar
MU=brentq(lambda x:x-1-math.exp(-x),1,2,xtol=1e-15)
ZETA=math.log1p(math.log(2))/math.log(2)
ALPHA=(9-MU)/24
BETA=MU-1-math.log(MU)
def ff(q):
 if q>=.125:return q-math.log1p(q)
 term=q*q;total=term/2
 for k in range(3,100):
  term*=-q;v=term/k;total+=v
  if abs(v)<abs(total)*1e-17:break
 return total

def hm(q):
 if q>=.125:return (1+1/q)*math.log1p(q)-1
 term=q;total=term/2
 for k in range(2,100):
  term*=-q;v=term/(k*(k+1));total+=v
  if abs(v)<abs(total)*1e-17:break
 return total

def run(j):
 n=j['N'];kind=j['kind'];L=math.log(n);qmin=1/(ZETA*(n-1));crit=n*ff(qmin)
 tau={'zero':0.,'half_transition':crit/2,'transition':crit,'double_transition':2*crit,'inverseN2':1/n**2,'one':1.,'logN':L,'sqrtNL':math.sqrt(n*L),'linear01':n/100,'leading_frontier':n*BETA}[kind]
 if tau<=crit:q=qmin;branch='cap';corr=n*hm(q)+tau*(1+1/q)
 else:q=brentq(lambda u:ff(u)-tau/n,qmin,max(2.,2*tau/n+2),xtol=1e-15,rtol=2e-15);branch='interior';corr=n*q
 p=1+1/q
 def obj(z):
  a=math.exp(z);return n*hm(a)+tau*(1+1/a)
 qhi=max(8.,4*q,4*qmin);res=minimize_scalar(obj,bounds=(math.log(qmin),math.log(qhi)),method='bounded',options={'xatol':2e-12,'maxiter':500});gap=res.fun-corr
 assert abs(gap)<2e-6*max(1,abs(corr)) and p<=1+ZETA*(n-1)+1e-7
 z=brentq(lambda v:.5*math.erfc(v/math.sqrt(2))-math.sqrt(MU/n)*v,0,math.sqrt(L)+2,xtol=1e-14);d=math.sqrt(MU*n)*z;target=n*MU-d+ALPHA*L;minimum=n+n*hm(qmin);frontier=None;fb=None
 if target>=minimum:
  qt=target/n-1
  if qt>=qmin:frontier=n*ff(qt);fb='interior'
  else:frontier=(target-minimum)/(1+1/qmin);fb='cap'
 expansion=n*BETA-(MU-1)/MU*d+(((MU-1)/MU)*ALPHA+1/(2*MU))*L
 return dict(j,tau=tau,q=q,p=p,branch=branch,transition_tau=crit,log_bound=n+corr,correction_over_N=corr,q_equation_residual=ff(q)-tau/n if branch=='interior' else None,direct_minimizer_log_bound=n+res.fun,direct_minimizer_gap=gap,direct_minimizer_success=bool(res.success),reference_target=target,reference_tau_frontier=frontier,frontier_branch=fb,frontier_expansion=expansion,frontier_expansion_error_over_logN=None if frontier is None else (frontier-expansion)/L,scope='optimized Hardy upper; reference target is asymptotic model, not exact finite Q_N')
if __name__=='__main__':
 out=Path(sys.argv[1]);out.mkdir(exist_ok=True,parents=True);jobs=[dict(N=n,kind=k) for n in [2,4,8,16,32,64,128,256,1024,4096,16384,65536,262144,1048576] for k in ['zero','half_transition','transition','double_transition','inverseN2','one','logN','sqrtNL','linear01','leading_frontier']];start=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(run,jobs))
 (out/'optimize76.json').write_text(json.dumps(dict(mu=MU,zeta=ZETA,beta_frontier=BETA,workers=12,seconds=time.monotonic()-start,rows=rows),indent=2)+'\n');print('PASS',len(rows),'corner-threshold optimizations')
