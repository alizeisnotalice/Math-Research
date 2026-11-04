import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import sys,math,time,json,argparse
from pathlib import Path
from fractions import Fraction as Q
from concurrent.futures import ProcessPoolExecutor,as_completed
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent))
import run as base
import numpy as np
from scipy.optimize import differential_evolution

def score(A,tau):
 H=base.subsets(A);return base.massweighted(H,tau)/sum(float(m)*abs(float(x))**(tau-1) for x,m in A)
def small_chain(n):
 if n%2==0:
  return [(Q(-1)+Q(2*j,n-1),Q(4,3*(n-1))) for j in range(n)]
 return [(Q(j,n),Q(2,3*n)) for j in range(1,n+1)]
def two_with_tiny(n):
 xs=[Q(-1)+Q(2*j,n-1) for j in range(1,n-1)];xs=[Q(1,3) if x==0 else x for x in xs]
 return [(Q(-1),Q(4,3))]+[(x,Q(1,4096)) for x in xs]+[(Q(1),Q(4,3))]
def to_v(A):
 A=sorted(A); non=[(x,m) for x,m in A if x!=1]; anchor=[m for x,m in A if x==1];assert len(anchor)==1
 return [float(x) for x,m in non]+[math.log2(float(m)) for x,m in non]+[math.log2(float(anchor[0]))]
def unpack(v,n):return [(float(x),2**float(m)) for x,m in zip(list(v[:n-1])+[1.],v[n-1:])]
def path_prev(tau,n):
 records=json.loads((P.parent/'phase1_seeds.json').read_text())
 matches=[r for r in records if r['tau']==tau and r['n']==n]
 if len(matches)!=1:raise ValueError(f'Expected one archived seed for tau={tau}, n={n}; found {len(matches)}')
 return matches[0]

def expand_previous(tau,n):
 prevA=[tuple(map(Q,x)) for x in path_prev(tau,n)['atoms']]
 # Preserve coalesced phase-1 measures in n slots by splitting mass at the same position.
 while len(prevA)<n:
  idx=next(i for i,(x,m) in enumerate(prevA) if x!=1 and m>Q(1,2048));x,m=prevA[idx];prevA[idx]=(x,m-Q(1,4096));prevA.append((x,Q(1,4096)))
 return prevA
def encode(A):return [[str(x),str(m)] for x,m in A]

def task(job):
 tau,n,it,pop=job; name=f'tau{round(100*tau):03d}_n{n:02d}'; file=P/'checkpoints'/f'{name}.json'
 if file.exists():return json.loads(file.read_text())
 started=time.time();pid=os.getpid();seed=290000+round(tau*100)*100+n; rng=np.random.default_rng(seed)
 prevA=expand_previous(tau,n)
 seeds=[('chain',small_chain(n)),('two_plus_tiny',two_with_tiny(n)),('phase1_best',prevA)]
 for name0,A in seeds:
  assert len(A)==n and (name0=='phase1_best' or len(set(x for x,m in A))==n) and max(x for x,m in A)==1 and all(x!=0 and m>0 for x,m in A),(name0,A)
 d=2*n-1;bounds=[(-1.,1.)]*(n-1)+[(-12.,5.)]*n; init=np.array([[rng.uniform(a,b) for a,b in bounds] for _ in range(pop*d)])
 for j,(_,A) in enumerate(seeds):init[j]=np.clip(to_v(A),[x[0] for x in bounds],[x[1] for x in bounds])
 calls=[0];invalid=[0]
 def obj(v):
  calls[0]+=1;A=unpack(v,n)
  if min(abs(x) for x,m in A)<1e-7:invalid[0]+=1;return 1e6
  v=base.objective_atoms(A,tau)
  if not math.isfinite(v):invalid[0]+=1;return 1e6
  return -v
 seed_float=[-obj(init[j]) for j in range(3)]; opt=differential_evolution(obj,bounds,init=init,seed=seed,maxiter=it,popsize=pop,tol=0,atol=0,polish=False)
 deval=float(-opt.fun);assert math.isfinite(deval) and deval+1e-10>=max(seed_float),(name,deval,seed_float)
 deA=base.rationalize(unpack(opt.x,n)); assert base.blocks(deA)==base.subsets(deA)
 candidates=[{'origin':'seed_'+name0,'atoms':encode(A),'score':score(A,tau)} for name0,A in seeds]
 candidates.append({'origin':'de_rationalized','atoms':encode(deA),'score':score(deA,tau)})
 chosen=max(candidates,key=lambda x:x['score']); A=[tuple(map(Q,x)) for x in chosen['atoms']];H=base.subsets(A)
 assert H==base.blocks(A) and chosen['score']+1e-14>=max(z['score'] for z in candidates[:3])
 out={'tau':tau,'n':n,'seed':seed,'status':'completed','worker_pid':pid,'evaluations':calls[0],'invalid_evaluations':invalid[0],'seed_float_scores':seed_float,'de_float_score':deval,'candidates':candidates,'best_origin':chosen['origin'],'best_score':chosen['score'],'halo':[[str(a),str(b)] for a,b in H],'seconds':time.time()-started,'message':str(opt.message)}
 tmp=file.with_suffix('.tmp');tmp.write_text(json.dumps(out,indent=2,allow_nan=False));tmp.replace(file);return out

def main():
 a=argparse.ArgumentParser();a.add_argument('--workers',type=int,default=15);a.add_argument('--maxiter',type=int,default=2000);a.add_argument('--popsize',type=int,default=24);a.add_argument('--check-seeds',action='store_true',help='validate the 25 archived seed geometries without running optimization');p=a.parse_args()
 if p.check_seeds:
  for tau in [.25,.5,.75,.9,1.]:
   for n in [2,3,4,6,8]:
    A=expand_previous(tau,n);assert len(A)==n and all(x!=0 and m>0 for x,m in A) and sum(x==1 for x,m in A)==1
    assert base.blocks(A)==base.subsets(A)
    assert len(to_v(A))==2*n-1
  print('Validated 25 archived seeds and exact halo comparisons');return
 (P/'checkpoints').mkdir(exist_ok=True)
 jobs=[(tau,n,p.maxiter,p.popsize) for tau in [.25,.5,.75,.9,1.] for n in [2,3,4,6,8]]
 print(f'START parent_pid={os.getpid()} workers={min(p.workers,15)} tasks={len(jobs)} maxiter={p.maxiter} popsize={p.popsize}',flush=True)
 with ProcessPoolExecutor(max_workers=min(p.workers,15)) as ex:
  fs={ex.submit(task,j):j for j in jobs}
  for f in as_completed(fs):
   try:
    x=f.result();print(f"DONE tau={x['tau']} n={x['n']} seed={x['seed']} worker={x['worker_pid']} calls={x['evaluations']} initial={max(x['seed_float_scores']):.12g} optimized={x['de_float_score']:.12g} final={x['best_score']:.12g} origin={x['best_origin']} seconds={x['seconds']:.1f}",flush=True)
   except Exception as e:print('FAILED',fs[f],repr(e),flush=True)
if __name__=='__main__':main()
