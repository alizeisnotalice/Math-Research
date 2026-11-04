import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import math,json,time,sys,argparse,random
from fractions import Fraction as Q
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scipy.optimize import differential_evolution, minimize_scalar
P=Path(__file__).resolve().parent

def merge(I):
 out=[]
 for a,b in sorted(I):
  if a>=b:continue
  if out and a<=out[-1][1]:out[-1]=(out[-1][0],max(b,out[-1][1]))
  else:out.append((a,b))
 return out

def blocks(atoms):
 z=sorted(atoms); out=[]
 for j,(xj,_) in enumerate(z):
  mass=0
  for xk,m in z[j:]:
   mass+=m; a=xk-mass/2;b=xj+mass/2
   if a<b:out.append((a,b))
 return merge(out)

def subsets(atoms):
 n=len(atoms); out=[]
 for mask in range(1,1<<n):
  S=[atoms[i] for i in range(n) if mask>>i&1]; m=sum(z[1] for z in S); a=max(x for x,_ in S)-m/2; b=min(x for x,_ in S)+m/2
  if a<b:out.append((a,b))
 return merge(out)

def powdiff(a,b,tau):
 if a==0:return b**tau
 return a**tau*math.expm1(tau*math.log1p((b-a)/a))
def massweighted(intervals,tau):
 s=0.
 for a,b in intervals:
  a=float(a);b=float(b)
  if b<=0:s+=powdiff(-b,-a,tau)
  elif a>=0:s+=powdiff(a,b,tau)
  else:s+=(-a)**tau+b**tau
 return s/tau

def objective_atoms(atoms,tau):
 H=blocks(atoms); den=sum(m*abs(x)**(tau-1) for x,m in atoms)
 return massweighted(H,tau)/den

def one(r,tau):
 return (((1+r)**tau- (1-r)**tau) if r<=1 else ((1+r)**tau+(r-1)**tau))/(2*tau*r)

def rationalize(atoms,bits=30):
 q=1<<bits; z={}
 for x,m in atoms:
  xx=Q(round(float(x)*q),q);mm=Q(round(float(m)*q),q)
  if xx==0 or mm<=0:raise ValueError('invalid rationalized atom')
  z[xx]=z.get(xx,Q(0))+mm
 return sorted(z.items())

def task(inp):
 tau,n,rep,it,pop=inp; seed=190000+round(100*tau)*1000+n*10+rep; file=P/'checkpoints'/f'tau{round(tau*100):03d}_n{n:02d}_r{rep}.json'
 if file.exists():return json.loads(file.read_text())
 start=time.time();rng=np.random.default_rng(seed);d=2*n-1; bounds=[(-1.,1.)]*(n-1)+[(-12.,5.)]*n
 def unpack(v):return [(float(x),2.**float(y)) for x,y in zip(list(v[:n-1])+[1.],v[n-1:])]
 calls=[0];rejected=[0]
 def obj(v):
  calls[0]+=1; atoms=unpack(v)
  if min(abs(x) for x,m in atoms)<1e-7:rejected[0]+=1;return 1e6
  val=objective_atoms(atoms,tau)
  if not math.isfinite(val):rejected[0]+=1;return 1e6
  return -val
 # Single atom is exact only in a separate reference: smallest added masses are 2^-12.
 best=minimize_scalar(lambda r:-one(r,tau),bounds=(1.000000001,15.999999),method='bounded')
 rr=float(best.x); sz=pop*d; init=np.array([[rng.uniform(a,b) for a,b in bounds] for _ in range(sz)])
 for i in range(min(12,sz)):
  init[i,:n-1]=rng.uniform(-1,1,n-1);init[i,n-1:2*n-2]=-12+rng.uniform(0,.2,n-1);init[i,-1]=math.log2(2*rr)
 result=differential_evolution(obj,bounds,init=init,seed=seed,maxiter=it,popsize=pop,tol=0,atol=0,polish=False)
 atoms=rationalize(unpack(result.x)); h1=blocks(atoms);h2=subsets(atoms); assert h1==h2
 rational_value=objective_atoms([(float(x),float(m)) for x,m in atoms],tau)
 out={'tau':tau,'n':n,'rep':rep,'seed':seed,'status':'completed','evaluations':calls[0],'invalid_evaluations':rejected[0],'float_value':float(-result.fun),'rational_value_float':rational_value,'atoms':[[str(x),str(m)] for x,m in atoms],'halo':[[str(a),str(b)] for a,b in h1],'single_reference_R':rr,'single_reference':one(rr,tau),'seconds':time.time()-start,'message':str(result.message)}
 tmp=file.with_suffix('.tmp');tmp.write_text(json.dumps(out,indent=2,allow_nan=False));tmp.replace(file);return out

def controls():
 rng=random.Random(19019)
 for n in range(2,9):
  for rep in range(30):
   A=[(Q(rng.choice([i for i in range(-100,101) if i]),100),Q(rng.randint(1,200),100)) for _ in range(n)]
   assert blocks(A)==subsets(A)
 for tau in [.25,.5,.75,.9,1.]:
  for R in [.25,.75,1.,1.25,2.,4.]:
   assert math.isclose(objective_atoms([(1.,2*R)],tau),one(R,tau),rel_tol=1e-13,abs_tol=1e-13)
 print('CONTROLS passed 210 rational geometry and 30 single atom tests',flush=True)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=15);ap.add_argument('--maxiter',type=int,default=600);ap.add_argument('--popsize',type=int,default=12);a=ap.parse_args();controls();(P/'checkpoints').mkdir(exist_ok=True)
 jobs=[(tau,n,rep,a.maxiter,a.popsize) for tau in [.25,.5,.75,.9,1.] for n in [2,3,4,6,8] for rep in range(3)]
 print(f'START parent_pid={os.getpid()} workers={min(a.workers,15)} tasks={len(jobs)} maxiter={a.maxiter} popsize={a.popsize}',flush=True)
 with ProcessPoolExecutor(max_workers=min(a.workers,15)) as ex:
  futs={ex.submit(task,x):x for x in jobs}
  for f in as_completed(futs):
   try:
    z=f.result();print(f"DONE tau={z['tau']} n={z['n']} rep={z['rep']} calls={z['evaluations']} ratio={z['rational_value_float']:.12f} single={z['single_reference']:.12f} seconds={z['seconds']:.1f}",flush=True)
   except Exception as e: print('FAILED',futs[f],repr(e),flush=True)
if __name__=='__main__':main()
