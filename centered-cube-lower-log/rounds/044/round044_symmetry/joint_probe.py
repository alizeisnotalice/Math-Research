"""Whole multi-table score/gate versus actual prior, before a proposed bound."""
from pathlib import Path
import json,math,time
import numpy as np
import probe as p
HERE=Path(__file__).resolve().parent

def main():
 target=HERE/'JOINT_SCAN.json';assert not target.exists();pf=HERE/'JOINT_PLAN.json';assert not pf.exists()
 config=dict(cases=[dict(n=1,grid=512),dict(n=2,grid=128),dict(n=3,grid=64)],lengths=[1.2,1.5,1.8,2.4],batch=4096,selector='count k -> table k mod4')
 p.dump(pf,dict(config=config,script_SHA256=p.sha(__file__),source_SHA256=p.sha(p.__file__)))
 start=time.time();out=[]
 for case in config['cases']:
  n=case['n'];grid=case['grid'];total=grid**n;s=np.arange(n+1)%4;pi=np.bincount(s,weights=[math.comb(n,k)/2**n for k in range(n+1)],minlength=4)
  c=p.C.copy();c['lengths']=config['lengths'];sums=np.zeros((4,len(c['lengths'])));maxcs=-np.inf
  for at in range(0,total,config['batch']):
   idx=np.arange(at,min(at+config['batch'],total));xx=np.empty((len(idx),n));v=idx.copy()
   for i in range(n):xx[:,i]=(v%grid+.5)/grid;v//=grid
   feat,fixed,_=p.features(xx,c);F=p.old.value(feat,s);G=np.einsum('b,bnj->nj',pi,fixed)
   separated=np.zeros_like(F)
   for b in range(4):
    if pi[b]:
     comp=feat[b][...,s==b].sum(-1);separated+=comp*comp/(pi[b]*fixed[b])
   maxcs=max(maxcs,float(np.max(F*F/G-separated)))
   sums+=np.array([(F*F/G).sum(0),separated.sum(0),abs(F-G).sum(0),F.sum(0)])
  out.append(dict(n=n,grid=grid,total=total,pi=pi.tolist(),whole_F2_over_G=(sums[0]/total).tolist(),CS_separated_integral=(sums[1]/total).tolist(),L1_difference=(sums[2]/total).tolist(),F_mass_quadrature=(sums[3]/total).tolist(),max_pointwise_CS_violation=maxcs))
 result=dict(config=config,cases=out,plan_SHA256=p.sha(pf),seconds=time.time()-start,
 scope='Deterministic grid exploration of actual multi-table coupled quantity; not outward integral/Cn certificate')
 p.dump(target,result);print(json.dumps(result,indent=2))
if __name__=='__main__':main()
