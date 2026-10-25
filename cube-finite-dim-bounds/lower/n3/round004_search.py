"""Unit lattice plus genuinely small free satellite: finite joint-union search."""
from pathlib import Path
import json,time,itertools
import numpy as np
from scipy.stats import qmc
from scipy.optimize import minimize
from witness import union
ROOT=Path(__file__).resolve().parent
BASE=ROOT.parents[1]/'baseline/original_certificates.json'
cert=json.loads(BASE.read_text())[0];u=2**(1/3);v=3**(1/3)
abase=np.array([1+u/2,(2+u*u*v)/4,(1+u)/2])

def decode(z):
 z=np.asarray(z);return abase*np.exp(z[:3]),np.array([[0.,0.,0.],z[3:6]]),np.array([1.,np.exp(z[6])])

def make_boxes(z,keys):
 p,ph,w=decode(z);kl=np.array([k[:6] for k in keys]).reshape(-1,2,3);kh=np.array([k[6:] for k in keys]).reshape(-1,2,3)
 cnt=np.prod(np.maximum(0,kh-kl+1),axis=2);m=cnt@w;r=m**(1/3)/2
 high=(kh+ph[None,:,:])*p;low=(kl+ph[None,:,:])*p;high[cnt==0]=-np.inf;low[cnt==0]=np.inf
 lo=np.maximum(np.max(high,axis=1)-r[:,None],0);hi=np.minimum(np.min(low,axis=1)+r[:,None],p)
 good=np.all(hi>lo+1e-13,axis=1)
 return np.stack([lo[good],hi[good]],axis=2),[k for k,g in zip(keys,good) if g]

def value(z,keys):
 z=np.asarray(z)
 if max(abs(z[:3]))>.08 or min(z[3:6])<0 or max(z[3:6])>1 or not np.log(2**-12)<=z[6]<=np.log(2**-3):return 0.
 b,kk=make_boxes(z,keys);return union(b)[0]/(1+np.exp(z[6])) if len(b) else 0.

def primary_keys():
 keys=[]
 for lengths in cert['witness_counts']:
  for lo in itertools.product(range(-3,4),repeat=3):
   hi=[lo[j]+lengths[j]-1 for j in range(3)];k=tuple(lo)+(0,0,0)+tuple(hi)+(-1,-1,-1);keys.append(k)
 # Include all old witnesses potentially active under allowed period perturbations.
 out=set()
 for shift in [-.08,0,.08]:
  z=np.zeros(7);z[:3]=shift;z[3:6]=.5;z[6]=np.log(2**-8);b,kk=make_boxes(z,keys);out.update(kk)
 # Singleton satellite boxes for all possible unit-cell-intersecting translations.
 for k in itertools.product(range(-1,2),repeat=3):out.add((0,0,0)+tuple(k)+(-1,-1,-1)+tuple(k))
 return sorted(out)

def capture(z,t):
 p,ph,w=decode(z);q=np.mod(t[:,None,:]-ph[None,:,:],1);rr=np.abs((q[:,:,:,None]-np.arange(-3,4))*p[None,None,:,None]).reshape(len(t),-1);keys=set()
 for j in range(rr.shape[1]):
  r=rr[:,j];lo=np.ceil(t[:,None,:]-ph[None,:,:]-r[:,None,None]/p-1e-11).astype(int);hi=np.floor(t[:,None,:]-ph[None,:,:]+r[:,None,None]/p+1e-11).astype(int);cnt=np.prod(np.maximum(0,hi-lo+1),axis=2);mass=cnt@w
  good=(mass>8*r**3+1e-10)&(r<=2.8)
  for l,h in zip(lo[good],hi[good]):
   for b in range(2):
    if np.any(h[b]<l[b]):l[b]=[0,0,0];h[b]=[-1,-1,-1]
   keys.add(tuple(map(int,np.concatenate([l.ravel(),h.ravel()]))))
 return sorted(keys)

def main():
 start=time.monotonic();base=primary_keys();pool=set(base);t=qmc.Sobol(3,scramble=True,seed=20261025402).random_base2(12)
 phases=qmc.Sobol(3,scramble=True,seed=20261025401).random_base2(5).tolist()+[[.5,.5,.5],[.5,0,0],[0,.5,0],[0,0,.5],[.5,.5,0],[.5,0,.5],[0,.5,.5],[.25,.25,.25],[.75,.75,.75]]
 cases=[]
 for i,ph in enumerate(phases):
  for eps in [2**-8,2**-6,2**-4]:
   z=np.array([0.,0.,0.]+ph+[np.log(eps)]);kk=sorted(set(base)|set(capture(z,t)));pool.update(kk);val=value(z,kk);cases.append({'phase_id':i,'phase':ph,'epsilon':eps,'z':z.tolist(),'joint_union':val,'keys':kk})
 print('scan',len(cases),'pool',len(pool),'best',max(c['joint_union'] for c in cases),flush=True)
 # Any apparent excess is saved before optimizer continuation for exact prioritization.
 from fractions import Fraction as F
 target=float(F(cert['lower_bound_interval'][0]));excess=[c for c in cases if c['joint_union']>target+1e-12]
 if excess:
  c=max(excess,key=lambda x:x['joint_union']);save({'status':'apparent excess; stop and exact-certify immediately','initial_cases':[{k:v for k,v in q.items() if k!='keys'} for q in cases],'optimization_evaluations':0,'elapsed_seconds':time.monotonic()-start},c['z'],c['keys']);return
 pool=sorted(pool);starts=[]
 for c in sorted(cases,key=lambda c:c['joint_union'],reverse=True):
  if all(np.linalg.norm(np.array(c['phase'])-np.array(s['phase']))>.1 for s in starts):starts.append(c)
  if len(starts)==4:break
 runs=[];count=0
 for c in starts:
  def obj(z):
   nonlocal count;count+=1;return -value(z,pool)
  res=minimize(obj,c['z'],method='Nelder-Mead',options={'maxfev':400,'xatol':1e-8,'fatol':1e-10});z=res.x;before=value(z,pool);rec=capture(z,qmc.Sobol(3,scramble=True,seed=20261025403).random_base2(14));newkeys=sorted(set(pool)|set(rec));after=value(z,newkeys)
  row={'phase_id':c['phase_id'],'start':c['z'],'z':z.tolist(),'before_recapture':before,'after_recapture':after,'extra_keys':len(set(rec)-set(pool)),'nfev':res.nfev,'stop':str(res.message),'keys':newkeys};runs.append(row);print({k:v for k,v in row.items() if k!='keys'},flush=True)
  if after>target+1e-12:save({'status':'apparent excess; stop and exact-certify immediately','initial_cases':[{k:v for k,v in q.items() if k!='keys'} for q in cases],'runs':[{k:v for k,v in q.items() if k!='keys'} for q in runs],'optimization_evaluations':count,'elapsed_seconds':time.monotonic()-start},z,newkeys);return
 best=max(runs,key=lambda r:r['after_recapture']);keys=best['keys']
 def obj2(z):
  nonlocal count;count+=1;return -value(z,keys)
 res=minimize(obj2,best['z'],method='Nelder-Mead',options={'maxfev':400,'xatol':1e-8,'fatol':1e-10});z=res.x;before=value(z,keys);rec=capture(z,qmc.Sobol(3,scramble=True,seed=20261025404).random_base2(15));finalkeys=sorted(set(keys)|set(rec));after=value(z,finalkeys)
 out={'status':'numerical finite joint unions only','initial_phase_seed':20261025401,'initial_capture_seed':20261025402,'recapture_seeds':[20261025403,20261025404],'initial_cases':[{k:v for k,v in q.items() if k!='keys'} for q in cases],'primary_key_count':len(base),'pooled_key_count':len(pool),'runs':[{k:v for k,v in q.items() if k!='keys'} for q in runs],'optimization_evaluations':count,'elapsed_seconds':time.monotonic()-start,'last_before_recapture':before,'last_after_recapture':after,'last_extra_keys':len(set(rec)-set(keys)),'last_nfev':res.nfev,'last_stop':str(res.message)}
 save(out,z,finalkeys)

def save(out,z,keys):
 out.setdefault('initial_phase_seed',20261025401);out.setdefault('initial_capture_seed',20261025402);out.setdefault('independent_recapture_seed',20261025403);out.setdefault('initial_receivers_per_case',4096);out.setdefault('recapture_receivers_per_survivor',16384)
 p,ph,w=decode(z);bb,active=make_boxes(z,keys);out['final']={'z':np.asarray(z).tolist(),'periods':p.tolist(),'phases':ph.tolist(),'weights':w.tolist(),'joint_union':value(z,keys),'active_key_count':len(active)}
 # Recoverable compact receipt is written before key serialization.
 (ROOT/'round004_results.json').write_text(json.dumps(out,indent=2)+'\n');(ROOT/'round004_keys.json').write_text(json.dumps([list(map(int,k)) for k in active],separators=(',',':'))+'\n');print('FINAL',out['final'],flush=True)
if __name__=='__main__':main()
