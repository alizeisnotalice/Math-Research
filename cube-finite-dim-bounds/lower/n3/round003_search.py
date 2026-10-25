"""Three-layer staircase capture, joint union optimization, and recapture."""
from pathlib import Path
import json,time
import numpy as np
from scipy.stats import qmc
from scipy.optimize import minimize
from witness import union
ROOT=Path(__file__).resolve().parent
u=2**(1/3);v=3**(1/3)
a=np.array([1+u/2,(2+u*u*v)/4,(1+u)/2]);pbase=a*np.array([3,1,1])

def decode(z):
 z=np.asarray(z);p=pbase*np.exp(z[:3]);ph=np.array([[0.,0.,0.],[1/3,z[3],z[4]],[2/3,z[5],z[6]]]);w=np.array([1+z[7],1+z[8],1-z[7]-z[8]])
 return p,ph,w

def boxes(z,keys):
 p,ph,w=decode(z);kl=np.array([k[:9] for k in keys]).reshape(-1,3,3);kh=np.array([k[9:] for k in keys]).reshape(-1,3,3)
 cnt=np.prod(np.maximum(0,kh-kl+1),axis=2);mass=cnt@w;r=mass**(1/3)/2
 high=(kh+ph[None,:,:])*p;low=(kl+ph[None,:,:])*p;high[cnt==0]=-np.inf;low[cnt==0]=np.inf
 lo=np.maximum(np.max(high,axis=1)-r[:,None],0);hi=np.minimum(np.min(low,axis=1)+r[:,None],p)
 ok=np.all(hi>lo+1e-13,axis=1)
 return np.stack([lo[ok],hi[ok]],axis=2),[k for k,kok in zip(keys,ok) if kok]

def value(z,keys):
 z=np.asarray(z)
 if np.max(np.abs(z[:3]))>.07 or np.max(np.abs(z[3:7]))>.5 or np.max(np.abs(z[7:]))>.3 or abs(sum(z[7:]))>.4:return 0.
 b,kk=boxes(z,keys)
 return union(b)[0]/3 if len(b) else 0.

def capture(z,t):
 p,ph,w=decode(z);q=np.mod(t[:,None,:]-ph[None,:,:],1)
 rr=np.abs((q[:,:,:,None]-np.arange(-3,4))*p[None,None,:,None]).reshape(len(t),-1);keys=set()
 for j in range(rr.shape[1]):
  r=rr[:,j];kl=np.ceil(t[:,None,:]-ph[None,:,:]-r[:,None,None]/p-1e-11).astype(int);kh=np.floor(t[:,None,:]-ph[None,:,:]+r[:,None,None]/p+1e-11).astype(int)
  cnt=np.prod(np.maximum(0,kh-kl+1),axis=2);mass=cnt@w;good=(mass>8*r**3+1e-10)&(r<=2.8)
  for lo,hi in zip(kl[good],kh[good]):
   for b in range(3):
    if np.any(hi[b]<lo[b]):lo[b]=[0,0,0];hi[b]=[-1,-1,-1]
   keys.add(tuple(map(int,np.concatenate([lo.ravel(),hi.ravel()]))))
 return sorted(keys)

def main():
 start=time.monotonic();t=qmc.Sobol(3,scramble=True,seed=20261025303).random_base2(12);cases=[];pool=set()
 zs=[('baseline',np.zeros(9))]
 for mode in ['monotone','opposite','rotated']:
  for f in [.01,.04,.1,.2]:
   for mass in [[0,0],[.1,-.05],[-.1,.1]]:
    z=np.zeros(9)
    z[3:7]={'monotone':[f,0,2*f,0],'opposite':[f,0,-f,0],'rotated':[f,0,0,f]}[mode]
    z[7:]=mass;zs.append((f'{mode}_{f}_{mass}',z))
 for name,z in zs:
  keys=capture(z,t);val=value(z,keys);pool.update(keys);cases.append({'name':name,'z':z.tolist(),'joint_union':val,'key_count':len(keys)})
 print('initial cases',len(cases),'pool',len(pool),'best',max(c['joint_union'] for c in cases),flush=True)
 pool=sorted(pool);runs=[];count=0
 for c in sorted(cases[1:],key=lambda c:c['joint_union'],reverse=True)[:4]:
  def obj(z):
   nonlocal count;count+=1;return -value(z,pool)
  res=minimize(obj,c['z'],method='Nelder-Mead',options={'maxfev':350,'xatol':1e-8,'fatol':1e-10})
  z=res.x;old=value(z,pool);newkeys=capture(z,qmc.Sobol(3,scramble=True,seed=20261025304).random_base2(14));expanded=sorted(set(pool)|set(newkeys));new=value(z,expanded)
  row={'start':c['name'],'z':z.tolist(),'before_recapture':old,'after_recapture':new,'extra_keys':len(set(newkeys)-set(pool)),'nfev':res.nfev,'stop':str(res.message),'keys':[list(k) for k in expanded]};runs.append(row)
  print({k:v for k,v in row.items() if k!='keys'},flush=True)
 best=max(runs,key=lambda r:r['after_recapture']);keys=best['keys']
 def obj2(z):
  nonlocal count;count+=1;return -value(z,keys)
 res=minimize(obj2,best['z'],method='Nelder-Mead',options={'maxfev':400,'xatol':1e-8,'fatol':1e-10})
 z=res.x;before=value(z,keys);rec=capture(z,qmc.Sobol(3,scramble=True,seed=20261025305).random_base2(14));finalkeys=sorted(set(map(tuple,keys))|set(rec));after=value(z,finalkeys)
 p,ph,w=decode(z);bb,active=boxes(z,finalkeys)
 out={'status':'numerical finite-witness unions only','initial_seed':20261025303,'recapture_seeds':[20261025304,20261025305],'initial_cases':cases,'initial_pool_size':len(pool),'runs':[{k:v for k,v in r.items() if k!='keys'} for r in runs],'optimization_evaluations':count,'elapsed_seconds':time.monotonic()-start,'final':{'z':z.tolist(),'periods':p.tolist(),'phases':ph.tolist(),'weights':w.tolist(),'before_recapture':before,'joint_union':after,'extra_keys':len(set(rec)-set(map(tuple,keys))),'active_key_count':len(active),'nfev_final':res.nfev,'stop':str(res.message)}}
 (ROOT/'round003_keys.json').write_text(json.dumps([list(k) for k in active],separators=(',',':'))+'\n');(ROOT/'round003_results.json').write_text(json.dumps(out,indent=2)+'\n')
 print(out['final'],flush=True)
if __name__=='__main__':main()
