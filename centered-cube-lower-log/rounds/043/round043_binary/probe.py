"""Real periodic K2 source, count selector, positive polynomial DP. Float pilot only."""
from pathlib import Path
import json, math, time, hashlib, argparse
import numpy as np
from scipy.stats import beta
HERE=Path(__file__).resolve().parent
CFG=dict(K=2,R=4,dimensions=[32,64,128],epsilon=.08,delta=.04,
 phases=[[.15,.35],[.35,.85],[.65,.15],[.85,.65]],
 lengths=[.25,.4,.6,.8,1.,1.04,1.08,1.15,1.25,1.4,1.6,1.8,2.,2.1,2.3,2.5,2.7,2.9],
 thresholds=[1.05,1.25,1.5,2.,3.,5.,8.],train_N=512,test_N=8192,batch=256,
 train_seed=202610294310,test_seed=202610294311,temperatures=[.2,.07],sweeps=3,
 guard=1.e-9,alpha_total=.05)
def dump(path,obj):path.write_text(json.dumps(obj,indent=2)+'\n')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def overlap(x,L,left,right):
 out=np.zeros_like(x)
 for shift in range(-3,4):out+=np.maximum(0.,np.minimum(x+L/2,right+shift)-np.maximum(x-L/2,left+shift))
 return out

def features(x,c):
 N,n=x.shape;J=len(c['lengths']);R=c['R'];feat=np.empty((R,N,J,n+1));fixed=np.empty((R,N,J));dp_error=0.
 for b in range(R):
  z=np.empty((n,N,J,2))
  for j in range(2):
   center=(j+c['phases'][b][j])/2
   assert j/2<center-c['delta']/2<center+c['delta']/2<(j+1)/2
   for li,L in enumerate(c['lengths']):
    if float(L).is_integer():z[:,:,li,j]=.5
    else:z[:,:,li,j]=(c['epsilon']*overlap(x,L,j/2,(j+1)/2)+(1-c['epsilon'])/(2*c['delta'])*overlap(x,L,center-c['delta']/2,center+c['delta']/2)).T/L
  a=z.sum(-1);fixed[b]=a.prod(0)
  # Direct nonnegative coefficient recurrence: coefficient k is real window
  # integral restricted to count(label 1)=k. This is PB DP after factoring a.
  cur=np.ones((N,J,1))
  for i in range(n):
   nxt=np.zeros((N,J,i+2));nxt[:,:,:-1]+=cur*z[i,:,:,0,None];nxt[:,:,1:]+=cur*z[i,:,:,1,None];cur=nxt
  feat[b]=cur
  dp_error=max(dp_error,float(np.max(np.abs(cur.sum(-1)-fixed[b])/(1+fixed[b]))))
 return feat,fixed,dp_error

def value(feat,s):
 out=np.zeros(feat.shape[1:3])
 for b in range(feat.shape[0]):out+=feat[b][...,s==b].sum(-1)
 return out

def fit(feat,fixed,c,n):
 rng=np.random.default_rng(c['train_seed']+n);w=np.array([math.comb(n,k)*2.**(-n)for k in range(n+1)]);best=None;records=[]
 for lam in c['thresholds']:
  s=rng.integers(c['R'],size=n+1);v=value(feat,s)
  for tau in c['temperatures']:
   for sweep in range(c['sweeps']):
    changed=0
    for k in rng.permutation(n+1):
     old=s[k];base=v-feat[old,:,:,k]
     candidates=base[None,:,:]+feat[:,:,:,k]
     scores=(1/(1+np.exp(-np.clip((candidates.max(-1)-lam)/(tau*lam),-40,40)))).mean(-1)
     b=int(np.argmax(scores))
     if scores[b]>scores[old]+1e-13:s[k]=b;v=base+feat[b,:,:,k];changed+=1
    if not changed:break
  pi=np.bincount(s,weights=w,minlength=c['R']);vv=value(feat,s).max(-1)
  score=lam*np.mean(vv>lam+c['guard'])
  row=dict(lambda_=lam,weak=float(score),prior=float(lam*np.mean(np.einsum('b,bnj->nj',pi,fixed).max(-1)>lam+c['guard'])),best_fixed=float(max(lam*np.mean(q.max(-1)>lam+c['guard'])for q in fixed)))
  records.append(row)
  if best is None or score>best[0]:best=(score,s.copy(),lam,pi)
 return best,records

def run():
 ap=argparse.ArgumentParser();ap.add_argument('--cost',action='store_true');args=ap.parse_args();c=CFG.copy()
 if args.cost:
  tt=time.time();x=np.random.default_rng(9900).random((64,32));_,_,err=features(x,c);out=dict(n=32,N=64,seconds=time.time()-tt,DP_relative_partition_error=err);dump(HERE/'COST.json',out);print(out);return
 assert not (HERE/'SUMMARY.json').exists()
 dump(HERE/'PLAN.json',dict(config=c,script_SHA256=sha(__file__),scope='Uniform pilot; all lambda > 1. CP intervals only ideal iid point-query Bernoulli; no interval FP or cell certificate. No Cn growth claim.',names=['gated','actual_prior','fixed0','fixed1','fixed2','fixed3'],selection='largest guarded training weak among frozen threshold list; deterministic smooth coordinate search'))
 allrows=[];start=time.time()
 for n in c['dimensions']:
  nt=time.time();trainx=np.random.default_rng(c['train_seed']+n).random((c['train_N'],n));feat,fixed,err=features(trainx,c);best,records=fit(feat,fixed,c,n);_,s,lam,pi=best
  # Include exact binomial rational prior numerators; denominator=2^n.
  numerator=[sum(math.comb(n,k)for k in range(n+1)if s[k]==b)for b in range(c['R'])]
  frozen=dict(n=n,lambda_=lam,selector=s.tolist(),prior_mass=pi.tolist(),prior_exact_numerators=numerator,prior_exact_denominator=2**n,training=records,config_SHA256=sha(HERE/'PLAN.json'))
  frozenpath=HERE/f'FROZEN_n{n}.json';dump(frozenpath,frozen)
  del trainx,feat,fixed
  rng=np.random.default_rng(c['test_seed']+n);counts=np.zeros(6,dtype=int);raw=np.zeros(6,dtype=int);near=np.zeros(6,dtype=int);all_lam=np.zeros((6,len(c['thresholds'])),dtype=int);top=np.zeros(6);sumv=np.zeros(6);near_examples=[]
  for st in range(0,c['test_N'],c['batch']):
   N=min(c['batch'],c['test_N']-st);x=rng.random((N,n));feat,fixed,e=features(x,c);err=max(err,e)
   maxima=np.concatenate([value(feat,s).max(-1)[None,:],np.einsum('b,bnj->nj',pi,fixed).max(-1)[None,:],fixed.max(-1)])
   counts+=(maxima>lam+c['guard']).sum(-1);raw+=(maxima>lam).sum(-1);near+=(np.abs(maxima-lam)<=c['guard']).sum(-1);top=np.maximum(top,maxima.max(-1));sumv+=maxima.sum(-1)
   all_lam+=(maxima[:,:,None]>np.array(c['thresholds'])[None,None,:]+c['guard']).sum(1)
   for mod,q in np.argwhere(np.abs(maxima-lam)<=c['guard']):near_examples.append(dict(model=int(mod),sample_index=int(st+q),value=float(maxima[mod,q]),x=x[q].tolist()))
  rows=[];alpha_each=c['alpha_total']/(2*6*len(c['dimensions']));N=c['test_N']
  for i,name in enumerate(['gated','actual_prior']+[f'fixed{b}'for b in range(4)]):
   k=int(counts[i]);lo=0. if k==0 else float(beta.ppf(alpha_each,k,N-k+1));hi=1. if k==N else float(beta.ppf(1-alpha_each,k+1,N-k))
   rows.append(dict(model=name,k_guard=k,k_raw=int(raw[i]),near_guard=int(near[i]),weak=float(lam*k/N),ideal_iid_CP_weak_CI=[lam*lo,lam*hi],largest_observed=float(top[i]),mean_max=float(sumv[i]/N)))
  row=dict(n=n,lambda_=lam,rows=rows,guard=c['guard'],exact_binomial_prior_numerators=numerator,prior_mass=pi.tolist(),best_fixed=max(r['weak']for r in rows[2:]),gate_minus_bestfixed=rows[0]['weak']-max(r['weak']for r in rows[2:]),gate_minus_prior=rows[0]['weak']-rows[1]['weak'],all_thresholds_diagnostic=dict(thresholds=c['thresholds'],successes=all_lam.tolist(),scope='postfreeze diagnostic; no selection claim or simultaneous CI'),DP_partition_error=err,near_examples=near_examples,seconds=time.time()-nt,frozen_SHA256=sha(frozenpath))
  allrows.append(row);dump(HERE/f'RESULT_n{n}.json',row);print(json.dumps(dict(n=n,lambda_=lam,seconds=row['seconds'],weak=[r['weak']for r in rows],near=near.tolist())),flush=True)
 dump(HERE/'SUMMARY.json',dict(config=c,results=allrows,seconds=time.time()-start,full_interval_certified=False,Cn_gain_proved=False,source_mass=1.,plan_SHA256=sha(HERE/'PLAN.json'),meaning='Exploration only. CP covers ideal iid guarded float-classifier probability; FP integration and continuum query-cell uncertainty unverified.'))
if __name__=='__main__':run()
