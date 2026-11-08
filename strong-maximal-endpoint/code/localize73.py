import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import sys,json,math,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from scipy.optimize import brentq
from scipy.special import gammaln,logsumexp
MU=brentq(lambda x:x-1-math.exp(-x),1,2,xtol=1e-15)
def run(j):
 n=j['N'];M=j['M'];L=math.log(n);rho=math.sqrt(L/n);b=math.sqrt(MU/(2*M));p=b/rho;A=math.sqrt(MU)+math.sqrt(2*MU*M)+1
 assert p>1 and A*rho<math.sqrt(n)
 bp=brentq(lambda v:-p*math.expm1(-v/p)-1-math.exp(-v),0,64,xtol=1e-14)
 def mid(u):return n*bp+M*p*L+(p-1)*u-(n-1)*math.log1p(u)
 low=mid(A*rho);upper=mid(math.sqrt(n));large=M*L+n*math.log(4)-gammaln(n)+(n-1)*(math.log(math.sqrt(n)+n*math.log(2)+M*L)-math.log1p(math.sqrt(n)))
 H=float(logsumexp([max(low,upper),large]));repair=-n*math.log1p(-2*math.exp(-M*L));caplog=A*rho-math.log1p(-2*math.exp(-M*L))
 return dict(j,p=p,A=A,bp=bp,bp_root_residual=-p*math.expm1(-bp/p)-1-math.exp(-bp),mid_lower_log=low,mid_upper_log=upper,large_log=large,H_log=H,normalized_H=(H-n*MU)/math.sqrt(n*L),predicted_normalized_H=-math.sqrt(MU)-1,repair_log_cost=repair,repair_analytic_upper=4*math.exp((1-M)*L),repaired_cap_log=caplog,cap_below_2=caplog<math.log(2))
if __name__=='__main__':
 out=Path(sys.argv[1]);out.mkdir(exist_ok=True,parents=True);jobs=[dict(N=2**k,M=m) for k in [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] for m in [2,3,4]];start=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(run,jobs))
 (out/'localize73.json').write_text(json.dumps(dict(mu=MU,seconds=time.monotonic()-start,workers=12,rows=rows),indent=2)+'\n');print('PASS',len(rows),'localization bounds')
