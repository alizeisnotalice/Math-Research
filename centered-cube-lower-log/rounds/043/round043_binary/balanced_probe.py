"""Second, preregistered-after-primary pilot. Count-reflection balanced R2 soft selector."""
import probe as p
from scipy.stats import beta
import numpy as np
import time,json,math
HERE=p.HERE
C=p.CFG.copy();C.update(R=2,phases=p.CFG['phases'][:2],train_seed=202610294320,test_seed=202610294321)
def coeff(feat,s,n):
 v=.5*(feat[0,:,:,n//2]+feat[1,:,:,n//2])
 for k,b in enumerate(s):v+=feat[b,:,:,k]+feat[1-b,:,:,n-k]
 return v

def main():
 assert not (HERE/'BALANCED_SUMMARY.json').exists()
 p.dump(HERE/'BALANCED_PLAN.json',dict(config=C,script_SHA256=p.sha(__file__),parent_probe_SHA256=p.sha(p.__file__),constraint='For k<n/2 choose b(k), at n-k choose 1-b(k); at midpoint soft weights=(1/2,1/2). Binomial symmetry gives prior=(1/2,1/2) exactly.',scope='Added after primary collapse observed. New independent train/test seeds; no test-guided retuning. No interval/cell certificate.'))
 summary=[];started=time.time()
 for n in C['dimensions']:
  st=time.time();rng=np.random.default_rng(C['train_seed']+n);feat,fixed,err=p.features(rng.random((C['train_N'],n)),C);best=None;records=[]
  for lam in C['thresholds']:
   s=rng.integers(2,size=n//2);v=coeff(feat,s,n)
   for tau in C['temperatures']:
    for sweep in range(C['sweeps']):
     changed=0
     for k in rng.permutation(n//2):
      b=s[k];other=v-feat[b,:,:,k]-feat[1-b,:,:,n-k]+feat[1-b,:,:,k]+feat[b,:,:,n-k]
      util=lambda vv:np.mean(1/(1+np.exp(-np.clip((vv.max(-1)-lam)/(tau*lam),-40,40))))
      if util(other)>util(v)+1e-13:s[k]=1-b;v=other;changed+=1
     if not changed:break
   score=lam*np.mean(v.max(-1)>lam+C['guard']);records.append(dict(lambda_=lam,training_weak=float(score)))
   if best is None or score>best[0]:best=(score,s.copy(),lam)
  _,s,lam=best;weights=np.zeros((2,n+1));weights[:,n//2]=.5
  for k,b in enumerate(s):weights[b,k]=1;weights[1-b,n-k]=1
  w=np.array([math.comb(n,k)*2.**(-n)for k in range(n+1)]);assert np.max(np.abs(weights@w-.5))<2e-15
  p.dump(HERE/f'BALANCED_FROZEN_n{n}.json',dict(n=n,lambda_=lam,selector_weights=weights.tolist(),paired_choice=s.tolist(),prior_exact=['1/2','1/2'],training=records,plan_SHA256=p.sha(HERE/'BALANCED_PLAN.json')))
  del feat,fixed
  counts=np.zeros(4,dtype=int);raw=counts.copy();near=counts.copy();top=np.zeros(4);diagnostic=np.zeros((4,len(C['thresholds'])),dtype=int);near_examples=[];rng=np.random.default_rng(C['test_seed']+n)
  for at in range(0,C['test_N'],C['batch']):
   x=rng.random((min(C['batch'],C['test_N']-at),n));feat,fixed,e=p.features(x,C);err=max(err,e)
   vals=np.concatenate([coeff(feat,s,n).max(-1)[None,:],fixed.mean(0).max(-1)[None,:],fixed.max(-1)])
   counts+=(vals>lam+C['guard']).sum(-1);raw+=(vals>lam).sum(-1);near+=(np.abs(vals-lam)<=C['guard']).sum(-1);top=np.maximum(top,vals.max(-1));diagnostic+=(vals[:,:,None]>np.array(C['thresholds'])[None,None,:]+C['guard']).sum(1)
   for mod,q in np.argwhere(np.abs(vals-lam)<=C['guard']):near_examples.append(dict(model=int(mod),sample_index=int(at+q),value=float(vals[mod,q]),x=x[q].tolist()))
  rows=[];alpha_each=.05/(2*(6+4)*3);N=C['test_N']
  for i,name in enumerate(['balanced_gate','prior_half','fixed0','fixed1']):
   k=int(counts[i]);lo=0 if k==0 else beta.ppf(alpha_each,k,N-k+1);hi=1 if k==N else beta.ppf(1-alpha_each,k+1,N-k)
   rows.append(dict(model=name,k_guard=k,k_raw=int(raw[i]),near_guard=int(near[i]),weak=float(lam*k/N),ideal_iid_CP_weak_CI=[float(lam*lo),float(lam*hi)],largest_observed=float(top[i])))
  rec=dict(n=n,lambda_=lam,rows=rows,prior_exact=['1/2','1/2'],gate_minus_prior=rows[0]['weak']-rows[1]['weak'],gate_minus_bestfixed=rows[0]['weak']-max(r['weak']for r in rows[2:]),near_examples=near_examples,DP_partition_error=err,all_thresholds_diagnostic=dict(thresholds=C['thresholds'],counts=diagnostic.tolist()),seconds=time.time()-st)
  p.dump(HERE/f'BALANCED_RESULT_n{n}.json',rec);summary.append(rec);print(json.dumps(dict(n=n,lambda_=lam,seconds=rec['seconds'],weak=[r['weak']for r in rows],near=near.tolist())),flush=True)
 p.dump(HERE/'BALANCED_SUMMARY.json',dict(config=C,results=summary,seconds=time.time()-started,full_interval_certified=False,Cn_gain_proved=False,CP_scope='Ideal iid finite floating-classifier Bernoulli probability; alpha=.05/60 directions for 3 n times 10 models including primary, but primary published bounds use their original .05/36. No continuum claim.'))
if __name__=='__main__':main()
