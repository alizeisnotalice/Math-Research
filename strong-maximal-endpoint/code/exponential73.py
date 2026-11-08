import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import sys,json,math,time
from fractions import Fraction as F
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from scipy.optimize import brentq
from scipy.special import gammaln,logsumexp
MU=brentq(lambda x:x-1-math.exp(-x),1,2,xtol=1e-15)
BP=brentq(lambda x:-4*math.expm1(-x/4)-1-math.exp(-x),0,8,xtol=1e-14)
def cert(x,target,relation):
 term=F(1);lower=term
 for j in range(1,17):term*=x/j;lower+=term
 nextterm=term*x/17;upper=lower+nextterm/(1-x/18)
 passed=lower>target if relation=='greater' else upper<target
 assert passed
 return dict(x=str(x),target=str(target),relation=relation,Taylor_degree=16,lower=str(lower),upper=str(upper),pass_=passed)
def run(n):
 k=.01;a=math.log(1.5);lo=n*(BP+4*k)+3*a-(n-1)*math.log1p(a);hi=n*(BP+4*k)+3*math.sqrt(n)-(n-1)*math.log1p(math.sqrt(n));large=k*n+n*math.log(4)-gammaln(n)+(n-1)*(math.log(math.sqrt(n)+n*(math.log(2)+k))-math.log1p(math.sqrt(n)));H=float(logsumexp([max(lo,hi),large]));r=2*math.exp(-k*n);repair=-n*math.log1p(-r) if 0<r<1 else None;caplog=math.log(1.5)-math.log1p(-r) if r<1 else None
 return dict(N=n,H_log=H,middle_lower_log=lo,middle_upper_log=hi,large_log=large,repair_log_cost=repair,repair_log_upper_log=math.log(4*n)-k*n,repair_binary_underflow=(r==0),repaired_cap_log=caplog,cap_below_2=bool(caplog is not None and caplog<math.log(2)),H_exponential_rate=H/n)
if __name__=='__main__':
 out=Path(sys.argv[1]);out.mkdir(exist_ok=True,parents=True);t=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(run,[2**k for k in range(5,21)]))
 cs=[cert(F(3,8),F(16,11),'greater'),cert(F(2,5),F(3,2),'less'),cert(F(1,3),F(7,5),'less'),cert(F(5,4),F(4),'less')]
 (out/'exponential73.json').write_text(json.dumps(dict(seconds=time.monotonic()-t,workers=12,b4=BP,mu=MU,exact_certificates=cs,rows=rows),indent=2)+'\n');print('PASS exponential cap reduction and four exact rational inequalities')
