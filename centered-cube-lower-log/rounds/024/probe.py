import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['VECLIB_MAXIMUM_THREADS']='1'
from pathlib import Path
import numpy as np,json,time
from fractions import Fraction as Q
from math import isqrt
P=Path(__file__).resolve().parent
src=(P.parent/'017/probe.py').read_text();exec(src.split('rows=[]')[0],globals());P=Path(__file__).resolve().parent
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--grid-size',type=int,choices=[262144,1048576],default=262144);parser.add_argument('--replications',type=int,default=8);parser.add_argument('--seed',type=int,default=2026102907);parser.add_argument('--output',default='results.json');args=parser.parse_args()
start=time.perf_counter();G=args.grid_size;a=.5/K;n=900000000;x=(np.arange(G)+.5)/G;pp=prim(x);fp=np.fft.fft(pp);auto=np.fft.ifft(abs(fp)**2).real/G
kk=np.arange(int(np.ceil((10/9)*G/2)),int(np.floor((11/9)*G/2))+1);L=2*kk/G;c=np.sqrt(2*b/3);sins=np.sin(np.pi*L[:,None]*m);cosx=np.cos(2*np.pi*x[:,None]*m)
psin=pp@np.sin(2*np.pi*x[:,None]*m)/G
variance=(2*(np.mean(pp*pp)-auto[(2*kk)%G])+4*(sins*sins)@(c*psin)+.5*(sins*sins)@(c*c))/(L*L)
# Explicit covariance control at sparse scales.
ids=np.linspace(0,len(kk)-1,17).astype(int);varerr=0.;rng=np.random.default_rng(args.seed);records=[]
for rep in range(args.replications):
 winner=int(rng.integers(len(kk)));k=int(kk[winner]);ll=L[winner]
 X=(np.roll(pp,-k)-np.roll(pp,k)+cosx@(c*np.sin(np.pi*m*ll)))/ll
 q=(1+a*X)/G;q/=q.sum();counts=rng.multinomial(n,q);corr=np.fft.ifft(np.conj(np.fft.fft(counts))*fp).real
 lin=(corr[kk%G]-corr[(-kk)%G]+sins@(c*(counts@cosx)))/L
 approx=a*lin-n*a*a*variance/2;errors=[];linearerr=0
 for j in ids:
  k=int(kk[j]);ll=L[j];X=(np.roll(pp,-k)-np.roll(pp,k)+cosx@(c*np.sin(np.pi*m*ll)))/ll
  direct=float(counts@np.log1p(a*X));errors.append(direct-approx[j]);varerr=max(varerr,abs(float(np.mean(X*X))-variance[j]));linearerr=max(linearerr,abs(float(counts@X)-lin[j]))
 records.append(dict(rep=rep,winner=winner,max_approx=float(approx.max()),max_sparse_log_error=float(np.max(np.abs(errors))),sparse_log_errors=errors,max_linear_error=linearerr))
A=(154/27)/(10/9);qmax=float(variance.max());t=110.;dQ=np.sqrt(2*n*A*A*qmax*t)+2*A*A*t/3;err=a*a*dQ/2+a**3*A*(n*qmax+dQ)/(3*(1-a*A))
# Conservative q<=2 for this uniform grid, with exact rational verification.
counts4=[255,256,256,256];uu=[Q(2*j+1,8) for j in range(4)]
EP2=sum(c4*(u4*u4-u4+Q(1,3)) for c4,u4 in zip(counts4,uu))/K
Sbound=sum((4*Q(729,29)*(min(Q(10,27),Q(1,3*j))**4+3*min(Q(1,9),Q(1,3*j))**4) for j in range(1,17)),Q(0))
q_bound=(8*(EP2+Q(2*K,G))+2*Sbound/3)*Q(81,100);assert q_bound<2
aa=Q(1,2*K);AA=Q(77,15);ell=110;varQ=2*n*AA*AA*2*ell
root=Q(isqrt((varQ.numerator*10**24)//varQ.denominator)+1,10**12)
dQ_bound=root+2*AA*AA*ell/3
error_bound=aa*aa*dQ_bound/2+aa**3*AA*(n*2+dQ_bound)/(3*(1-aa*AA))
safe_error=Q((error_bound.numerator*1000+error_bound.denominator-1)//error_bound.denominator,1000)
out=dict(certified_q_upper=str(q_bound),conservative_q=2,certified_log_error_upper=str(safe_error),seconds=time.perf_counter()-start,G=G,r=len(kk),n=n,a=a,variance_min=float(variance.min()),variance_max=qmax,variance_formula_max_error=varerr,A=A,log_union_budget=t,failure_bound=float(2*len(kk)*np.exp(-t)),uniform_log_error_bound=float(err),records=records)
(P/args.output).write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k!='records'},indent=2));print([(r['max_sparse_log_error'],r['max_linear_error']) for r in records])
