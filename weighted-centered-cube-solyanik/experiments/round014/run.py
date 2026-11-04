import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
 os.environ[key]='1'
import math, json, time, argparse
from fractions import Fraction as Q
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from scipy.optimize import differential_evolution
import halo13

P=Path(__file__).resolve().parent

def exact_halo(E,alpha):
 ends=[v for a,b in E for v in (a,b)]; ell=sum(b-a for a,b in E); lo=ends[0]-ell/(2*alpha); hi=ends[-1]+ell/(2*alpha); pieces=list(E)
 for e in ends:
  cuts=sorted({lo,hi,*((e+f)/2 for f in ends if lo<(e+f)/2<hi)})
  def g(x):
   l=min(e,2*x-e); r=max(e,2*x-e)
   return sum(max(Q(0),min(r,b)-max(l,a)) for a,b in E)-2*alpha*abs(x-e)
  for a,b in zip(cuts,cuts[1:]):
   ga,gb=g(a),g(b)
   if ga>0 and gb>0:pieces.append((a,b))
   elif ga>0:pieces.append((a,a+(b-a)*ga/(ga-gb)))
   elif gb>0:pieces.append((a+(b-a)*(-ga)/(gb-ga),b))
 out=[]
 for a,b in sorted(pieces):
  if a>=b:continue
  if out and a<=out[-1][1]:out[-1]=(out[-1][0],max(out[-1][1],b))
  else:out.append((a,b))
 return out

def difference(H,E):
 out=[]
 for a,b in H:
  cur=a
  for c,d in E:
   if d<=cur:continue
   if c>=b:break
   if c>cur:out.append((cur,min(c,b)))
   cur=max(cur,d)
   if cur>=b:break
  if cur<b:out.append((cur,b))
 return out

def posdiff(a,b,tau):
 if a==0:return b**tau
 return a**tau * math.expm1(tau*math.log1p((b-a)/a))
def weight(E,tau):
 s=0.
 for a,b in E:
  a,b=float(a),float(b)
  if b<=0:s+=posdiff(-b,-a,tau)
  elif a>=0:s+=posdiff(a,b,tau)
  else:s+=(-a)**tau+b**tau
 return s

def score(E,t,tau,exact=False):
 H=exact_halo(E,1-Q(1,round(1/t))) if exact else halo13.halo(E,1-t)
 return weight(difference(H,E),tau)/weight(E,tau)/(t**tau),H

def task(x):
 tau,k,rep,it,pop=x; seed=140000+int(100*tau)*1000+k*10+rep; file=P/'checkpoints'/f'tau{int(tau*100)}_k{k}_seed{rep}.json'
 if file.exists():return json.loads(file.read_text())
 t=2**(-k); rng=np.random.default_rng(seed); baseline=1/((1+t)**tau-t**tau); start=time.time(); results=[]
 for n in (3,4):
  bounds=[(.15,1.85),(.000001,3.5),(.000001,3.5)]+[(1e-10,3.5),(.000000001,.999999999)]*(n-2)
  dim=len(bounds); sz=pop*dim
  init=np.array([[rng.uniform(a,b) for a,b in bounds] for j in range(sz)])
  root=[1.,1/(1+t),1/(1+t)]+[v for j in range(n-2) for v in (1e-10,.999999999)]
  init[0]=root
  for j in range(1,min(24,sz)):
   # Mixed tiny holes and close perturbations of the symmetric benchmark.
   row=np.array(root); row[:3]*=np.exp(rng.normal(0,.03,3));
   for h in range(n-2):row[3+2*h]=10**rng.uniform(-9,-1);row[4+2*h]=rng.uniform(.05,.999999999)
   init[j]=[min(max(z,a),b) for z,(a,b) in zip(row,bounds)]
  calls=[0]
  def obj(v):
   calls[0]+=1; E=halo13.mk(v,t,n)
   if E is None:return 1000.
   return -score(E,t,tau)[0]
  opt=differential_evolution(obj,bounds,init=init,seed=seed+n,maxiter=it,popsize=pop,tol=0,atol=0,polish=False)
  E=halo13.mk(opt.x,t,n); EF=[(Q(float(a)).limit_denominator(10**12),Q(float(b)).limit_denominator(10**12)) for a,b in E]
  valid=all(a<b for a,b in EF) and all(EF[j][1]<EF[j+1][0] for j in range(len(EF)-1))
  if valid:
   val,H=score(EF,t,tau,True)
   record={'segments':n,'evaluations':calls[0],'float_score':-opt.fun,'fraction_score':val,'intervals':[[str(a),str(b)] for a,b in EF],'halo':[[str(a),str(b)] for a,b in H],'parameters':opt.x.tolist()}
  else:record={'segments':n,'evaluations':calls[0],'float_score':-opt.fun,'rationalization_failed':True}
  results.append(record)
 record={'tau':tau,'k':k,'rep':rep,'seed':seed,'status':'completed','baseline':baseline,'results':results,'seconds':time.time()-start}
 tmp=file.with_suffix('.tmp');tmp.write_text(json.dumps(record,indent=2));tmp.replace(file)
 return record

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=15);ap.add_argument('--maxiter',type=int,default=1500);ap.add_argument('--popsize',type=int,default=24);a=ap.parse_args()
 jobs=[(tau,k,rep,a.maxiter,a.popsize) for tau in [.25,.5,.75] for k in [8,12,16] for rep in range(4)]
 print(f'RUN workers={min(a.workers,15)} jobs={len(jobs)} maxiter={a.maxiter} popsize={a.popsize} pid={os.getpid()}',flush=True)
 with ProcessPoolExecutor(max_workers=min(a.workers,15)) as ex:
  fs={ex.submit(task,j):j for j in jobs}
  for f in as_completed(fs):
   try:
    r=f.result();print(f"DONE tau={r['tau']} k={r['k']} rep={r['rep']} seed={r['seed']} values={[z.get('fraction_score') for z in r['results']]} baseline={r['baseline']:.12g} seconds={r['seconds']:.1f}",flush=True)
   except Exception as e:print('FAILED',fs[f],repr(e),flush=True)
if __name__=='__main__':main()
