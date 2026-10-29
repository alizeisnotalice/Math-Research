"""Frozen independent test replay: retain every near case for every preset threshold."""
import probe as p
import balanced_probe as b
import numpy as np
import json,time
H=p.HERE
records=[];started=time.time()
for balanced in [False,True]:
 c=b.C if balanced else p.CFG
 for n in c['dimensions']:
  f=json.loads((H/f'{"BALANCED_"if balanced else ""}FROZEN_n{n}.json').read_text());rng=np.random.default_rng(c['test_seed']+n);lam=np.array(c['thresholds']);models=4 if balanced else 6
  near=np.zeros((models,len(lam)),int);raw=np.zeros_like(near);guard=np.zeros_like(near);mind=np.full(near.shape,np.inf);examples=[]
  for at in range(0,c['test_N'],c['batch']):
   x=rng.random((min(c['batch'],c['test_N']-at),n));feat,fixed,_=p.features(x,c)
   if balanced:
    v=b.coeff(feat,np.array(f['paired_choice']),n);prior=fixed.mean(0)
   else:
    v=p.value(feat,np.array(f['selector']));prior=np.einsum('b,bnj->nj',f['prior_mass'],fixed)
   vals=np.concatenate([v.max(-1)[None,:],prior.max(-1)[None,:],fixed.max(-1)])
   d=vals[:,:,None]-lam[None,None,:];near+=(np.abs(d)<=c['guard']).sum(1);raw+=(d>0).sum(1);guard+=(d>c['guard']).sum(1);mind=np.minimum(mind,np.abs(d).min(1))
   for mod,q,j in np.argwhere(np.abs(d)<=c['guard']):examples.append(dict(model=int(mod),query_index=int(at+q),threshold=float(lam[j]),value=float(vals[mod,q]),x=x[q].tolist()))
  records.append(dict(balanced=balanced,n=n,thresholds=lam.tolist(),near_counts=near.tolist(),raw_counts=raw.tolist(),guard_counts=guard.tolist(),minimum_absolute_distance=mind.tolist(),all_near_examples=examples))
p.dump(H/'GUARD_REPLAY.json',dict(records=records,guard=1e-9,seconds=time.time()-started,scope='Same frozen seeds and selector. All preset thresholds; floats only, not outward intervals.'))
print(json.dumps(dict(seconds=time.time()-started,total_near=sum(sum(map(sum,r['near_counts']))for r in records),minimum_margin=min(min(map(min,r['minimum_absolute_distance']))for r in records))))
