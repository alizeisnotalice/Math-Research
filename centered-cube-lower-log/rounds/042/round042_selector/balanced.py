"""Exactly balanced cyclic type assignments, frozen extra exploration after initial pilot."""
from pathlib import Path
import json,sys,time,math
import numpy as np
import probe as p
HERE=Path(__file__).resolve().parent

def orbits(types):
 index={a:i for i,a in enumerate(types)};seen=set();out=[]
 for a in types:
  if a in seen:continue
  orbit=[a[r:]+a[:r] for r in range(4)]
  assert len(set(orbit))==4 # odd n rules out periods 1 and 2
  seen.update(orbit);out.append([index[t]for t in orbit])
 return out

def fit(feat,pr,types,c):
 orb=orbits(types);rng=np.random.default_rng(c['train_seed']+1);records=[];best=None
 for lam in c['thresholds']:
  choices=rng.integers(0,4,size=len(orb));s=np.empty(len(types),dtype=int)
  for ids,b in zip(orb,choices):s[ids]=(np.arange(4)+b)%4
  v=p.selector_value(feat,s);changes=0
  for tau in c['temperatures']:
   for sweep in range(c['optimizer_passes']):
    changed=0
    for oi in rng.permutation(len(orb)):
     ids=orb[oi];old=int(choices[oi]);base=v.copy()
     for idx in ids:base-=feat[s[idx],:,:,idx]
     additions=[];scores=[]
     for b in range(4):
      add=sum(feat[(r+b)%4,:,:,idx]for r,idx in enumerate(ids))
      additions.append(add);scores.append(p.utility(base+add,lam,tau))
     b=int(np.argmax(scores))
     if scores[b]>scores[old]+1e-14:
      choices[oi]=b;s[ids]=(np.arange(4)+b)%4;v=base+additions[b];changed+=1
    changes+=changed
    if not changed:break
  pi=np.bincount(s,weights=p.weights(types,c),minlength=4)
  assert np.array_equal(pi,np.full(4,.25)),pi
  score=float(lam*np.mean(v.max(axis=1)>lam));pv=pr.mean(axis=0).max(axis=1)
  records.append(dict(lambda_=lam,training_weak=score,training_prior_weak=float(lam*np.mean(pv>lam)),prior_mass=pi.tolist(),changes=changes))
  if best is None or score>best[0]:best=(score,s.copy(),lam,pi)
 return best,records

def main():
 globalstart=time.time();root=HERE/'balanced';root.mkdir(exist_ok=True)
 for n in [7,15]:
  folder=root/f'n{n}';folder.mkdir(exist_ok=True)
  assert not (folder/'RESULTS.json').exists()
  c=p.CONFIG.copy();c.update(n=n,train_seed=202610290430+n,test_seed=202610290440+n)
  plan=dict(config=c,driver_SHA256=p.sha(__file__),core_SHA256=p.sha(p.__file__),
      alpha_each=.05/24,scope='two preregistered odd dimensions, six models each, both CI directions; balanced type orbits; extension after n8 pilot')
  pf=folder/'PLAN.json';assert not pf.exists();pf.write_text(json.dumps(plan,indent=2)+'\n')
  started=time.time();rng=np.random.default_rng(c['train_seed']);feat,pr,types=p.features(rng.random((c['train_N'],n)),c)
  best,training=fit(feat,pr,types,c);_,s,lam,pi=best
  sf=folder/'FROZEN_SELECTOR.json';sf.write_text(json.dumps(dict(selector=s.tolist(),types=types,lambda_=lam,prior_mass=pi.tolist(),training=training,plan_SHA256=p.sha(pf)),indent=2)+'\n')
  del feat,pr;rng=np.random.default_rng(c['test_seed']);counts=np.zeros(6,dtype=int)
  for start in range(0,c['test_N'],c['test_batch']):
   N=min(c['test_batch'],c['test_N']-start);feat,pr,_=p.features(rng.random((N,n)),c)
   mx=np.concatenate([p.selector_value(feat,s).max(axis=1)[None,:],pr.mean(axis=0).max(axis=1)[None,:],pr.max(axis=2)],axis=0)
   counts+=(mx>lam).sum(axis=1)
  rows=[]
  for name,k in zip(['gated','prior','fixed0','fixed1','fixed2','fixed3'],counts):
   rows.append(dict(model=name,successes=int(k),N=c['test_N'],lambda_=lam,weak_estimate=float(lam*k/c['test_N']),weak_CI=[lam*t for t in p.ci(int(k),c['test_N'],plan['alpha_each'])]))
  out=dict(n=n,rows=rows,prior_mass=pi.tolist(),orbits=len(orbits(types)),types=len(types),
    plan_SHA256=p.sha(pf),selector_SHA256=p.sha(sf),seconds=time.time()-started,
    gated_minus_prior=float(lam*(counts[0]-counts[1])/c['test_N']),
    simultaneous_gap=rows[0]['weak_CI'][0]-rows[1]['weak_CI'][1],
    qualification='Float finite shared radii; balanced extension; no full FP certificate or C_n order claim',python=sys.version,numpy=np.__version__)
  (folder/'RESULTS.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
 print(json.dumps(dict(total_seconds=time.time()-globalstart)))
if __name__=='__main__':main()
