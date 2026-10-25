"""Extract finite joint cube witnesses and compute their union on the full cell."""
import json, time
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from search import ROOT,decode,eval_receivers

def extract(z,pow=13):
 t=qmc.Sobol(3,scramble=True,seed=202610255).random_base2(pow)
 p,phase,w=decode(z);ph=np.array([np.zeros(3),phase]);q=np.mod(t[:,None,:]-ph[None,:,:],1)
 rr=np.abs((q[:,:,:,None]-np.arange(-3,4))*p[None,None,:,None]).reshape(len(t),-1)
 keys=set()
 for j in range(rr.shape[1]):
  r=rr[:,j]
  kl=np.ceil(t[:,None,:]-ph[None,:,:]-r[:,None,None]/p-1e-11).astype(int)
  kh=np.floor(t[:,None,:]-ph[None,:,:]+r[:,None,None]/p+1e-11).astype(int)
  cnt=np.prod(np.maximum(0,kh-kl+1),axis=2); mass=cnt@w
  good=(mass>8*r**3+1e-10)&(r<=2.8)
  for lo,hi in zip(kl[good],kh[good]):
   # Canonical empty orbit marker.
   for b in range(2):
    if np.any(hi[b]<lo[b]):lo[b]=[0,0,0];hi[b]=[-1,-1,-1]
   keys.add(tuple(np.concatenate([lo.ravel(),hi.ravel()])))
 boxes=[];kept=[]
 for key in sorted(keys):
  lo=np.array(key[:6]).reshape(2,3);hi=np.array(key[6:]).reshape(2,3)
  cnt=np.prod(np.maximum(0,hi-lo+1),axis=1);mass=cnt@w;r=mass**(1/3)/2
  has=cnt>0
  lower=np.max((hi[has]+ph[has])*p,axis=0)-r
  upper=np.min((lo[has]+ph[has])*p,axis=0)+r
  lower=np.maximum(lower,0);upper=np.minimum(upper,p)
  if np.all(upper>lower):boxes.append(np.stack([lower,upper],axis=1));kept.append(list(map(int,key)))
 return np.array(boxes),kept

def union(boxes):
 # Exact coordinate decomposition modulo the floating-point endpoints.
 cuts=[np.unique(boxes[:,j,:].ravel()) for j in range(3)]
 sh=[len(c) for c in cuts];delta=np.zeros(sh,dtype=np.int32)
 for b in boxes:
  ix=[[int(np.searchsorted(c,b[j,k])) for k in range(2)] for j,c in enumerate(cuts)]
  for e in range(8):
   bits=[(e>>j)&1 for j in range(3)];idx=tuple(ix[j][bits[j]] for j in range(3));delta[idx]+=(-1)**sum(bits)
 for j in range(3):delta=np.cumsum(delta,axis=j)
 occ=delta[:-1,:-1,:-1]>0
 return float(np.einsum('ijk,i,j,k->',occ,np.diff(cuts[0]),np.diff(cuts[1]),np.diff(cuts[2]))),[len(c) for c in cuts]

def main():
 inp=json.loads((ROOT/'results.json').read_text());out=[]
 for c in inp['candidates']:
  start=time.monotonic();boxes,keys=extract(c['z']);vol,shape=union(boxes)
  row={k:c[k] for k in ['z','p','phase','weights','validation_value']};row.update({'joint_union_value':vol/2,'witness_count':len(keys),'cuts':shape,'elapsed_seconds':time.monotonic()-start,'witness_keys':keys})
  print({k:v for k,v in row.items() if k!='witness_keys'},flush=True);out.append(row)
 family=sorted({tuple(k) for c in out for k in c.pop('witness_keys')})
 (ROOT/'witness_family.json').write_text(json.dumps([list(k) for k in family],separators=(',',':'))+'\n')
 (ROOT/'witness_results.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
