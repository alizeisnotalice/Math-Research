"""Seven-parameter full-union refinement with updated-capture key families."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,time
import numpy as np
from scipy.stats import qmc
from scipy.optimize import minimize
from witness import union
ROOT=Path(__file__).resolve().parent
MASTER=Path('/Users/zhengzhihao/Desktop/T/Math-Research-finite-dim/cube-finite-dim-bounds')
SEED=json.loads((ROOT/'round004_certificate_input.json').read_text());P0=np.array(list(map(float,SEED['periods'])));PH0=np.array(list(map(float,SEED['phases'][1])));EPS0=float(F(SEED['weights'][1]));OLD=sorted(set(tuple(k) for k in SEED['witness_keys']))
TARGET=json.loads((MASTER/'audit/006-satellite_n3_results.json').read_text())['ratio_interval'][1]

def decode(z):
 z=np.asarray(z);return P0*np.exp(z[:3]),np.array([[0.,0.,0.],z[3:6]]),np.array([1.,np.exp(z[6])])

def make_boxes(z,keys):
 p,ph,w=decode(z);lo=np.array([k[:6] for k in keys]).reshape(-1,2,3);hi=np.array([k[6:] for k in keys]).reshape(-1,2,3);cnt=np.prod(np.maximum(0,hi-lo+1),axis=2);r=(cnt@w)**(1/3)/2
 low=(lo+ph[None,:,:])*p;high=(hi+ph[None,:,:])*p;low[cnt==0]=np.inf;high[cnt==0]=-np.inf
 low=np.maximum(np.max(high,axis=1)-r[:,None],0);high=np.minimum(np.min((lo+ph[None,:,:])*p+np.where(cnt[:,:,None]==0,np.inf,0),axis=1)+r[:,None],p)
 good=np.all(high>low+1e-13,axis=1);return np.stack([low[good],high[good]],axis=2),[k for k,g in zip(keys,good) if g]

def value(z,keys):
 z=np.asarray(z)
 if max(abs(z[:3]))>.04 or min(z[3:6])<0 or max(z[3:6])>1 or not np.log(2**-14)<=z[6]<=np.log(2**-2):return 0.
 b,kk=make_boxes(z,keys);return float(union(b)[0]/(1+np.exp(z[6]))) if len(b) else 0.

def capture(z,t):
 p,ph,w=decode(z);q=np.mod(t[:,None,:]-ph[None,:,:],1);rr=np.abs((q[:,:,:,None]-np.arange(-4,5))*p[None,None,:,None]).reshape(len(t),-1);keys=set()
 for j in range(rr.shape[1]):
  r=rr[:,j];lo=np.ceil(t[:,None,:]-ph[None,:,:]-r[:,None,None]/p-1e-11).astype(int);hi=np.floor(t[:,None,:]-ph[None,:,:]+r[:,None,None]/p+1e-11).astype(int);cnt=np.prod(np.maximum(0,hi-lo+1),axis=2);mass=cnt@w;good=(mass>8*r**3+1e-10)&(r<=3.2)
  for l,h in zip(lo[good],hi[good]):
   for b in range(2):
    if np.any(h[b]<l[b]):l[b]=[0,0,0];h[b]=[-1,-1,-1]
   keys.add(tuple(map(int,np.concatenate([l.ravel(),h.ravel()]))))
 return sorted(keys)

class Discovery(Exception):pass

def main():
 start=time.monotonic();z0=np.r_[np.zeros(3),PH0,np.log(EPS0)];rec=capture(z0,qmc.Sobol(3,scramble=True,seed=20261025501).random_base2(15));pool=sorted(set(OLD)|set(rec));seed_value=value(z0,pool);print('seed',seed_value,'old',len(OLD),'new',len(set(rec)-set(OLD)),'pool',len(pool),flush=True)
 count=0;history=[];runs=[];found=None
 starts=[z0,z0+np.array([0,0,0,.01,-.01,.005,np.log(2)]),z0+np.array([0,0,0,-.01,.01,-.005,-np.log(2)])]
 budgets=[1000,600,600]
 for stage,(zstart,budget) in enumerate(zip(starts,budgets)):
  best=[-1.,None];stagecount=0
  def obj(z):
   nonlocal count,stagecount,found
   val=value(z,pool);count+=1;stagecount+=1
   if val>best[0]:best[:]=[val,np.asarray(z).copy()]
   if val>float(F(TARGET))+1e-8:
    found=np.asarray(z).copy();raise Discovery()
   return -val
  simplex=np.tile(zstart,(8,1));steps=[.00005,.00005,.00005,.005,.005,.005,.1]
  for j in range(7):simplex[j+1,j]+=steps[j]
  try:
   res=minimize(obj,zstart,method='Nelder-Mead',options={'maxfev':budget,'initial_simplex':simplex,'xatol':1e-10,'fatol':1e-12});z=res.x;stop=str(res.message)
  except Discovery:z=found;stop='material full-union excess: interrupted optimization immediately'
  before=value(z,pool);extra=capture(z,qmc.Sobol(3,scramble=True,seed=20261025502).random_base2(15));expanded=sorted(set(pool)|set(extra));after=value(z,expanded);row={'stage':stage,'start':zstart.tolist(),'z':z.tolist(),'nfev':stagecount,'stop':stop,'before_recapture':before,'after_recapture':after,'extra_keys':len(set(extra)-set(pool))};runs.append(row);print(row,flush=True)
  if found is not None or after>float(F(TARGET))+1e-8:
   pool=expanded;break
  pool=expanded
 else:
  z=np.array(max(runs,key=lambda r:r['after_recapture'])['z']);extra=capture(z,qmc.Sobol(3,scramble=True,seed=20261025503).random_base2(16));pool=sorted(set(pool)|set(extra))
 p,ph,w=decode(z);b,active=make_boxes(z,pool)
 # Preserve original 170 as keys even when a moved box is empty.
 saved=sorted(set(active)|set(OLD));assert set(OLD)<=set(saved)
 out={'status':'material improvement; exact certification required' if value(z,pool)>float(F(TARGET))+1e-8 else 'no material finite-union gain','target_exact_upper':TARGET,'seed_input_sha256':hashlib.sha256((ROOT/'round004_certificate_input.json').read_bytes()).hexdigest(),'seed_value':seed_value,'initial_capture_seed':20261025501,'initial_extra_keys':len(set(rec)-set(OLD)),'initial_pool_size':len(set(OLD)|set(rec)),'recapture_seeds':[20261025502,20261025503],'optimization_evaluations':count,'runs':runs,'elapsed_seconds':time.monotonic()-start,'final':{'z':z.tolist(),'periods':p.tolist(),'phases':ph.tolist(),'weights':w.tolist(),'joint_union':value(z,pool),'active_key_count':len(active),'saved_key_count':len(saved),'original_170_preserved':set(OLD)<=set(saved)}}
 (ROOT/'round005_results.json').write_text(json.dumps(out,indent=2)+'\n')
 # 15 decimals retain the discovered margin; certificate never depends on float geometry.
 inp={'periods':[format(x,'.15f') for x in p],'phases':[['0','0','0'],[format(x,'.15f') for x in ph[1]]],'weights':['1',format(w[1],'.15f')],'witness_keys':[list(k) for k in saved]}
 (ROOT/'round005_certificate_input.json').write_text(json.dumps(inp,indent=2)+'\n');print('FINAL',out['final'],flush=True)
if __name__=='__main__':main()
