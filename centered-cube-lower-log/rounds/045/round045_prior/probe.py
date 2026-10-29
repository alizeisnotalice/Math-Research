"""Whole-vector prior-mixture importance pilot. Genuine periodic positive densities."""
from pathlib import Path
import json,time,math,hashlib,argparse
import numpy as np
from scipy.special import logsumexp
H=Path(__file__).resolve().parent
C=dict(dimensions=[32,128,512],epsilon=.02,uniform_proposal=.2,train_N=2048,test_N=32768,batch=256,train_seed=202610294510,test_seed=202610294511,
 log_thresholds=[0.,.5,1.,2.,4.,6.,8.,12.,16.,24.,32.],
 bank=[dict(m=1,phase=.07,width_factor=1.),dict(m=1,phase=.37,width_factor=.25),dict(m=2,phase=.11,width_factor=1.),dict(m=2,phase=.29,width_factor=.25)],
 mixture_candidates=[[.25,.25,.25,.25],[.5,.5,0.,0.],[0.,0.,.5,.5],[.5,0.,.5,0.],[0.,.5,0.,.5]],guard_log=1e-8,alpha=.05)
def dump(name,obj):(H/name).write_text(json.dumps(obj,indent=2)+'\n')
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def lengths(n):
 # Include small-grid anchors and radii approaching one from both directions.
 d=1/math.sqrt(n)
 return sorted(set([.4,.5,.6,.8,1-d,1.,1+d,1+2*d,1.4,1.8,2.2,2.6]))
def primitive(t,left,width):
 z=t-left;q=np.floor(z);return q*width+np.minimum(np.maximum(z-q,0.),width)
def kernel(x,L,b,n):
 delta=b['width_factor']/n;m=b['m'];out=np.zeros_like(x)
 assert delta<1/m
 for j in range(m):
  left=(b['phase']+j/m)%1
  out+=primitive(x+L/2,left,delta)-primitive(x-L/2,left,delta)
 return C['epsilon']+(1-C['epsilon'])*out/(m*delta*L)
def log_scores(x,n,ls):
 out=np.empty((len(x),len(C['bank']),len(ls)))
 for b,bank in enumerate(C['bank']):
  for j,L in enumerate(ls):
   a=kernel(x,L,bank,n);assert np.min(a)>0
   if float(L).is_integer():a=np.ones_like(a)
   out[:,b,j]=np.log(a).sum(1)
 return out

def sample(rng,N,n,ls,source_n=None):
 source_n=n if source_n is None else source_n
 R=len(C['bank']);J=len(ls);gamma=C['uniform_proposal'];comp=rng.choice(1+R*J,size=N,p=[gamma]+[(1-gamma)/(R*J)]*(R*J));x=rng.random((N,n))
 for c in np.unique(comp):
  if c==0:continue
  ids=np.flatnonzero(comp==c);b,j=divmod(int(c)-1,J);bank=C['bank'][b];m=bank['m'];delta=bank['width_factor']/source_n;shape=(len(ids),n)
  bg=rng.random(shape)<C['epsilon'];peak=(bank['phase']+rng.integers(m,size=shape)/m+delta*rng.random(shape))%1
  y=np.where(bg,rng.random(shape),peak);x[ids]=(y+ls[j]*(rng.random(shape)-.5))%1
 return x,comp

def log_q(scores):
 gamma=C['uniform_proposal'];R,J=scores.shape[1:]
 return np.logaddexp(math.log(gamma),math.log((1-gamma)/(R*J))+logsumexp(scores,axis=(1,2)))
def log_model(scores,pi):
 with np.errstate(divide='ignore'):lp=np.log(np.array(pi))
 return logsumexp(scores+lp[None,:,None],axis=1).max(1)
def bounded_ci(mean,variance,N,B,alpha):
 # Two-sided empirical Bernstein (Maurer-Pontil), union charged per direction conservatively.
 rad=math.sqrt(2*variance*math.log(2/alpha)/N)+7*B*math.log(2/alpha)/(3*(N-1))
 return [max(0.,mean-rad),min(B,mean+rad)],rad

def run():
 ap=argparse.ArgumentParser();ap.add_argument('--cost',action='store_true');args=ap.parse_args()
 if args.cost:
  start=time.time();n=128;ls=lengths(n);x,comp=sample(np.random.default_rng(4500),256,n,ls);s=log_scores(x,n,ls);lq=log_q(s);out=dict(n=n,N=256,R=4,J=len(ls),seconds=time.time()-start,logq_range=[float(lq.min()),float(lq.max())],RN_weight_max=float(np.exp(-lq).max()),scope='Cost only, no accuracy claim');dump('COST.json',out);print(out);return
 assert not(H/'SUMMARY.json').exists()
 dump('PLAN.json',dict(config=C,lengths_by_dimension={str(n):lengths(n)for n in C['dimensions']},script_SHA256=sha(__file__),cost_SHA256=sha(H/'COST.json'),qualification='Exploratory banks; not accepted sqrt(n) constructions. Full-vector mixture q covers every evaluated bank/scale, all source densities integrate 1.',selection='Train finite 5 mixtures and 11 log thresholds; freeze pair maximizing training IS weak. Independent test chosen mixture plus 4 fixed tables, all thresholds with charged family.'))
 results=[];started=time.time();dirs=2*len(C['dimensions'])*(5*len(C['log_thresholds'])+1);alpha=C['alpha']/dirs
 for n in C['dimensions']:
  st=time.time();ls=lengths(n);R=4;J=len(ls);hs=np.array(C['log_thresholds']);lam=np.exp(hs);B=R*J/(1-C['uniform_proposal']);rng=np.random.default_rng(C['train_seed']+n);means=np.zeros((5,len(hs)));near=np.zeros_like(means,dtype=int)
  for at in range(0,C['train_N'],C['batch']):
   x,_=sample(rng,min(C['batch'],C['train_N']-at),n,ls);s=log_scores(x,n,ls);lq=log_q(s)
   for m,pi in enumerate(C['mixture_candidates']):
    v=log_model(s,pi);d=v[:,None]-hs[None,:];y=np.exp(hs[None,:]-lq[:,None])*(d>C['guard_log']);assert y.max()<=B*(1+1e-10);means[m]+=y.sum(0);near[m]+=(np.abs(d)<=C['guard_log']).sum(0)
  means/=C['train_N'];mi,hi=np.unravel_index(np.argmax(means),means.shape);pi=C['mixture_candidates'][mi];selected_h=float(hs[hi]);dump(f'FROZEN_n{n}.json',dict(n=n,pi=pi,log_threshold=selected_h,training_IS_grid=means.tolist(),training_near_counts=near.tolist(),train_selected=float(means[mi,hi]),plan_SHA256=sha(H/'PLAN.json')))
  models=[pi]+np.eye(R).tolist();rng=np.random.default_rng(C['test_seed']+n);sums=np.zeros((5,len(hs)));sums2=sums.copy();hits=np.zeros_like(sums,dtype=int);rawhits=hits.copy();nearc=hits.copy();mind=np.full_like(sums,np.inf);maxy=np.zeros_like(sums);sumrn=0.;sumrn2=0.;maxrn=0.;examples={};hist=np.zeros(1+R*J,dtype=int)
  for at in range(0,C['test_N'],C['batch']):
   x,comp=sample(rng,min(C['batch'],C['test_N']-at),n,ls);hist+=np.bincount(comp,minlength=len(hist));s=log_scores(x,n,ls);lq=log_q(s);rn=np.exp(-lq);sumrn+=rn.sum();sumrn2+=(rn*rn).sum();maxrn=max(maxrn,float(rn.max()))
   for m,pp in enumerate(models):
    v=log_model(s,pp);d=v[:,None]-hs[None,:];z=d>C['guard_log'];yy=np.exp(hs[None,:]-lq[:,None])*z;assert yy.max()<=B*(1+1e-9);sums[m]+=yy.sum(0);sums2[m]+=(yy*yy).sum(0);hits[m]+=z.sum(0);rawhits[m]+=(d>0).sum(0);nearc[m]+=(np.abs(d)<=C['guard_log']).sum(0);mind[m]=np.minimum(mind[m],np.abs(d).min(0));maxy[m]=np.maximum(maxy[m],yy.max(0))
    for q,j in np.argwhere(np.abs(d)<=C['guard_log']):examples.setdefault(f'{m}_{j}',[]).append(int(at+q))
  N=C['test_N'];mean=sums/N;variance=np.maximum(0.,(sums2-sums*sums/N)/(N-1));rows=[]
  for m in range(5):
   bb=B*max(models[m]);ci=[bounded_ci(float(mean[m,j]),float(variance[m,j]),N,bb,alpha)[0]for j in range(len(hs))]
   rows.append(dict(model='mixture'if m==0 else f'fixed{m-1}',pi=models[m],IS_weak=mean[m].tolist(),sample_variance=variance[m].tolist(),EB_CI=ci,integrand_cap=bb,hits=hits[m].tolist(),raw_hits=rawhits[m].tolist(),near=nearc[m].tolist(),min_log_margin=mind[m].tolist(),max_integrand=maxy[m].tolist()))
  rnmean=sumrn/N;rnvar=(sumrn2-sumrn*sumrn/N)/(N-1);rnci=bounded_ci(rnmean,rnvar,N,1/C['uniform_proposal'],alpha)[0]
  row=dict(n=n,pi=pi,selected_log_threshold=selected_h,train_weak=float(means[mi,hi]),rows=rows,chosen_threshold_mixture_IS=float(mean[0,hi]),chosen_threshold_bestfixed_IS=float(mean[1:,hi].max()),RN_uniform_mass_mean=rnmean,RN_uniform_mass_EB_CI=rnci,RN_max=maxrn,proposal_component_hist=hist.tolist(),all_near_examples=examples,seconds=time.time()-st,scope='Floating whole-vector importance exploration, no interval/continuous-cell certificate; iid PRNG assumption.')
  dump(f'NEAR_INDICES_n{n}.json',dict(indices=examples,scope='All near query indices, regenerated from frozen seed; no query/features cache.'));row['all_near_examples']=dict(file=f'NEAR_INDICES_n{n}.json',count=sum(map(len,examples.values())));dump(f'RESULT_n{n}.json',row);results.append(row);print(json.dumps(dict(n=n,pi=pi,h=selected_h,mix=float(mean[0,hi]),fixed=float(mean[1:,hi].max()),RNmass=rnmean,seconds=row['seconds'])),flush=True)
 dump('SUMMARY.json',dict(config=C,results=results,seconds=time.time()-started,alpha_per_direction=alpha,directions=dirs,Cn_new_value_certified=False,accepted_sqrt_source=False,full_FP_certificate=False))
if __name__=='__main__':run()
