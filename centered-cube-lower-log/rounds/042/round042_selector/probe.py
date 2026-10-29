"""Finite-radius exact histogram polynomial exploration, genuine periodic box sources.
Float arithmetic, independent heldout queries. No asymptotic C_n claim.
"""
from pathlib import Path
import argparse,hashlib,json,math,time,sys
import numpy as np
HERE=Path(__file__).resolve().parent
CONFIG=dict(n=8,K=4,R=4,epsilon=.08,delta=.05,
 lengths=[.3,.45,.6,.8,1.2,1.5,1.8],train_N=2048,test_N=16384,
 train_seed=202610290420,test_seed=202610290421,
 thresholds=[.7,1.,1.5,2.,3.,5.,8.],optimizer_passes=6,
 temperatures=[.25,.12,.06],test_batch=256)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def states(n,K):
 if K==1:return [(n,)]
 return [(a,)+s for a in range(n+1) for s in states(n-a,K-1)]
def transitions(n,K):
 layers=[states(i,K) for i in range(n+1)];out=[]
 for k in range(1,n+1):
  prev={s:i for i,s in enumerate(layers[k-1])};tab=[]
  for j in range(K):
   dst=[];src=[]
   for a,s in enumerate(layers[k]):
    if s[j]:
     t=list(s);t[j]-=1;dst.append(a);src.append(prev[tuple(t)])
   tab.append((np.array(dst),np.array(src)))
  out.append(tab)
 return layers,out
def overlap(x,L,left,right):
 # Our L<2: all contributing period shifts lie in [-2,2].
 result=np.zeros_like(x)
 for shift in range(-2,3):
  result+=np.maximum(0.,np.minimum(x+L/2,right+shift)-np.maximum(x-L/2,left+shift))
 return result

def kernels(x,c):
 N,n=x.shape;K=c['K'];R=c['R'];ls=np.array(c['lengths'])
 out=np.empty((R,n,N,len(ls),K))
 for b in range(R):
  for j in range(K):
   # Equal cell mass; alternate phase order to avoid mere global translates.
   phase=((b+.5)/R if j%2==0 else (R-b-.5)/R)
   center=(j+phase)/K
   assert j/K<=center-c['delta']/2 and center+c['delta']/2<=(j+1)/K
   for li,L in enumerate(ls):
    out[b,:,:,li,j]=(c['epsilon']*overlap(x,L,j/K,(j+1)/K)
      +(1-c['epsilon'])/(K*c['delta'])*overlap(x,L,center-c['delta']/2,center+c['delta']/2)).T/L
 return out

def features(x,c):
 z=kernels(x,c);layers,tr=transitions(c['n'],c['K']);N=len(x);J=len(c['lengths'])
 out=np.empty((c['R'],N,J,len(layers[-1])))
 for b in range(c['R']):
  cur=np.ones((N,J,1))
  for i,tab in enumerate(tr):
   nxt=np.zeros((N,J,len(layers[i+1])))
   for j,(dst,src) in enumerate(tab):nxt[:,:,dst]+=cur[:,:,src]*z[b,i,:,:,j,None]
   cur=nxt
  out[b]=cur
 a=z.sum(axis=-1).prod(axis=1)
 assert np.allclose(out.sum(axis=-1),a,rtol=2e-13,atol=2e-13)
 return out,a,layers[-1]
def weights(types,c):
 n=c['n'];K=c['K'];return np.array([math.factorial(n)/math.prod(math.factorial(k)for k in a)/K**n for a in types])
def selector_value(feat,s):
 # Avoid NumPy mixed-advanced indexing reordering.
 result=np.zeros(feat.shape[1:3])
 for b in range(feat.shape[0]):result+=feat[b][...,s==b].sum(axis=-1)
 return result

def utility(v,lam,tau):
 t=(v.max(axis=1)-lam)/(tau*lam)
 return np.mean(1/(1+np.exp(-np.clip(t,-40,40))))
def fit(feat,prior,types,c):
 rng=np.random.default_rng(c['train_seed']+1);records=[];best=None
 for lam in c['thresholds']:
  # Balanced initial type assignment; no target/query-specific selector thereafter.
  s=rng.integers(0,c['R'],size=len(types));v=selector_value(feat,s)
  changes=0
  for tau in c['temperatures']:
   for sweep in range(c['optimizer_passes']):
    changed=0
    for t in rng.permutation(len(s)):
     old=int(s[t]);base=v-feat[old,:,:,t];scores=[]
     for b in range(c['R']):scores.append(utility(base+feat[b,:,:,t],lam,tau))
     b=int(np.argmax(scores))
     if scores[b]>scores[old]+1e-14:s[t]=b;v=base+feat[b,:,:,t];changed+=1
    changes+=changed
    if not changed:break
  score=float(lam*np.mean(v.max(axis=1)>lam));pi=np.bincount(s,weights=weights(types,c),minlength=c['R'])
  pv=np.einsum('b,bnj->nj',pi,prior).max(axis=1)
  row=dict(lambda_=lam,training_weak=score,training_prior_weak=float(lam*np.mean(pv>lam)),
     type_counts=np.bincount(s,minlength=c['R']).tolist(),prior_mass=pi.tolist(),changes=changes)
  records.append(row)
  if best is None or score>best[0]:best=(score,s.copy(),lam,pi)
 return best,records

def bernoulli_kl(p,q):
 if p==0:return -math.log1p(-q) if q<1 else math.inf
 if p==1:return -math.log(q) if q>0 else math.inf
 if q<=0 or q>=1:return math.inf
 return p*math.log(p/q)+(1-p)*math.log((1-p)/(1-q))
def ci(k,N,alpha):
 p=k/N;r=math.log(1/alpha)/N
 if k==0:lo=0.
 else:
  a,b=0.,p
  for _ in range(90):
   m=(a+b)/2
   if bernoulli_kl(p,m)>r:a=m
   else:b=m
  lo=b
 if k==N:hi=1.
 else:
  a,b=p,1.
  for _ in range(90):
   m=(a+b)/2
   if bernoulli_kl(p,m)>r:b=m
   else:a=m
  hi=a
 return [lo,hi]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',default=str(HERE/'RESULTS.json'));args=ap.parse_args()
 target=Path(args.output);assert not target.exists();c=CONFIG.copy();started=time.time()
 plan=HERE/'PLAN.json'
 payload=dict(config=c,script_SHA256=sha(__file__),scope='Finite n histogram selector; freeze before train/test. 2 test models+4 fixed tables at selected training lambda, 12 directions, alpha=.05/12. No rates.')
 if plan.exists():assert json.loads(plan.read_text())==payload
 else:plan.write_text(json.dumps(payload,indent=2)+'\n')
 rng=np.random.default_rng(c['train_seed']);x=rng.random((c['train_N'],c['n']))
 feat,pr,types=features(x,c);best,train=fit(feat,pr,types,c)
 _,s,lam,pi=best
 frozen=dict(lambda_=lam,selector=s.tolist(),types=[list(a)for a in types],prior_mass=pi.tolist(),training=train,
    plan_SHA256=sha(plan),selection='Maximum training weak value among seven fixed thresholds after deterministic coordinate search')
 sel=HERE/'FROZEN_SELECTOR.json';assert not sel.exists();sel.write_text(json.dumps(frozen,indent=2)+'\n')
 del feat,pr,x
 rng=np.random.default_rng(c['test_seed']);counts=np.zeros(2+c['R'],dtype=np.int64);sumv=np.zeros_like(counts,dtype=float)
 for start in range(0,c['test_N'],c['test_batch']):
  N=min(c['test_batch'],c['test_N']-start);x=rng.random((N,c['n']));feat,pr,_=features(x,c)
  maxima=np.concatenate([selector_value(feat,s).max(axis=1)[None,:],np.einsum('b,bnj->nj',pi,pr).max(axis=1)[None,:],pr.max(axis=2)],axis=0)
  counts+=(maxima>lam).sum(axis=1);sumv+=maxima.sum(axis=1)
 rows=[]
 for name,k,sv in zip(['gated','prior']+[f'fixed{b}'for b in range(c['R'])],counts,sumv):
  bounds=ci(int(k),c['test_N'],.05/(2*len(counts)))
  rows.append(dict(model=name,successes=int(k),N=c['test_N'],lambda_=lam,weak_estimate=lam*k/c['test_N'],weak_CI=[lam*t for t in bounds],mean_finite_max=float(sv/c['test_N'])))
 out=dict(config=c,plan_SHA256=sha(plan),selector_SHA256=sha(sel),rows=rows,
   gated_minus_prior=float(lam*(counts[0]-counts[1])/c['test_N']),
   simultaneous_gated_lower_minus_prior_upper=rows[0]['weak_CI'][0]-rows[1]['weak_CI'][1],
   confidence_scope='12 Chernoff-KL directions conditional on independent frozen training, finite shared radii, floating kernels not interval certified',
   periodic_mass=1.,L1_cutoff_transfer='See independent proof, large integer boxes, finite L boundary discarded',
   Cn_log_gain_proved=False,full_FP_certificate=False,python=sys.version,numpy=np.__version__,seconds=time.time()-started)
 target.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
