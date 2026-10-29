"""Table-specific mirror labels; true periodic boxes, mass-normalized count weights."""
from pathlib import Path
import argparse,importlib.util,json,math,hashlib,time,sys
import numpy as np
HERE=Path(__file__).resolve().parent
BASE=HERE.parent/'round043_binary/probe.py'
if not BASE.exists():BASE=Path(__file__).resolve().parents[2]/'043/round043_binary/probe.py'
spec=importlib.util.spec_from_file_location('generic_fit43',BASE);old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
C=dict(R=4,dimensions=[32,64,128],epsilon=.08,delta=.05,
 theta=[0.,.125,.25,.375],offsets=[.125,.1875,.25,.3125],
 lengths=[.25,.4,.6,.8,1.,1.04,1.08,1.15,1.25,1.4,1.5,1.6,1.8,2.,2.1,2.3,2.5,2.7,2.9],
 thresholds=[1.25,1.5,2.,3.,5.,8.],train_N=512,test_N=8192,batch=256,
 train_seed=202610294410,test_seed=202610294411,temperatures=[.2,.07],sweeps=3,guard=1e-9)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(path,value):path.write_text(json.dumps(value,indent=2)+'\n')
def components(b,c=C,broken=False):
 theta=c['theta'][b];d=c['offsets'][b];eps=c['epsilon'];delta=c['delta'];out=[]
 for l,r,j in [(0.,theta,0),(theta,theta+.5,1),(theta+.5,1.,0)]:
  if r>l:out.append((l,r,eps,j))
 for sign,j in [(1,1),(-1,0)]:
  center=(theta+sign*d)%1;w=delta*(2 if broken and sign<0 else 1);height=(1-eps)/(2*w)
  l=center-w/2;r=center+w/2
  segments=[(l,r)] if l>=0 and r<=1 else ([(0.,r),(1+l,1.)] if l<0 else [(l,1.),(0.,r-1)])
  out.extend((ll,rr,height,j)for ll,rr in segments)
 return out
def kernels(x,lengths,c=C,broken=False):
 N,n=x.shape;out=np.zeros((c['R'],n,N,len(lengths),2))
 for b in range(c['R']):
  for li,L in enumerate(lengths):
   if float(L).is_integer():out[b,:,:,li,:]=.5;continue
   for left,right,height,j in components(b,c,broken):
    v=np.zeros_like(x)
    for shift in range(-3,4):v+=np.maximum(0.,np.minimum(x+L/2,right+shift)-np.maximum(x-L/2,left+shift))
    out[b,:,:,li,j]+=height*v.T/L
 return out
def features(x,c=C):
 N,n=x.shape;J=len(c['lengths']);z=kernels(x,c['lengths'],c);feat=np.empty((c['R'],N,J,n+1));fixed=z.sum(-1).prod(1);err=0.
 for b in range(c['R']):
  cur=np.ones((N,J,1))
  for i in range(n):
   nxt=np.zeros((N,J,i+2));nxt[:,:,:-1]+=cur*z[b,i,:,:,0,None];nxt[:,:,1:]+=cur*z[b,i,:,:,1,None];cur=nxt
  feat[b]=cur;err=max(err,float(np.max(abs(cur.sum(-1)-fixed[b])/(1+fixed[b]))))
 return feat,fixed,err
def scan():
 target=HERE/'SCAN.json';assert not target.exists();start=time.time()
 cfg=dict(N=32768,lengths=[.3,.6,.9,1.2,1.5,1.8,2.4],tilts=[-32,-8,-1,0,1,8,32],broken_controls=True)
 pf=HERE/'SCAN_PLAN.json';assert not pf.exists();dump(pf,dict(config=cfg,source=C,script_SHA256=sha(__file__),fit_helper_SHA256=sha(BASE)))
 x=((np.arange(cfg['N'])+.5)/cfg['N'])[:,None];records=[]
 for broken in [False,True]:
  z=kernels(x,cfg['lengths'],broken=broken)[:,0];a=z.sum(-1);p=z/a[...,None]
  for b in range(4):
   theta=C['theta'][b];xr=((2*theta-x)%1);zr=kernels(xr,cfg['lengths'],broken=broken)[b,0]
   ar=zr.sum(-1);pr=zr/ar[...,None]
   for li,L in enumerate(cfg['lengths']):
    for t in cfg['tilts']:
     logs=(t+1)*np.log(a[b,:,li]);w=np.exp(logs-logs.max());w/=w.sum();pp=p[b,:,li]
     mu=w@pp;H=(pp.T*w)@pp*2;ev=np.linalg.eigvalsh(H)
     records.append(dict(b=b,broken=broken,L=L,t=t,mu_plus=float(mu[1]),gamma=float(ev[-1]),
       old_trace=float(np.trace(H)),second_eigen=float(ev[0]),
       max_score_reflection_error=float(np.max(abs(a[b,:,li]-ar[:,li]))),
       max_posterior_swap_error=float(np.max(abs(pp-pr[:,li,::-1])))))
 good=[r for r in records if not r['broken']];bad=[r for r in records if r['broken']]
 summary=dict(plan_SHA256=sha(pf),cases=len(records),mirror_max_mu_error=max(abs(r['mu_plus']-.5)for r in good),
  mirror_max_gamma_error=max(abs(r['gamma']-1)for r in good),
  mirror_max_score_reflection_error=max(r['max_score_reflection_error']for r in good),
  mirror_max_posterior_swap_error=max(r['max_posterior_swap_error']for r in good),
  broken_max_mu_error=max(abs(r['mu_plus']-.5)for r in bad),broken_max_gamma=max(r['gamma']for r in bad),
  representatives=sorted(good,key=lambda r:r['second_eigen'],reverse=True)[:8]+sorted(bad,key=lambda r:abs(r['mu_plus']-.5),reverse=True)[:4],
  seconds=time.time()-start,python=sys.version,numpy=np.__version__,scope='Deterministic physical-source exploration; no C_n estimate or outward quadrature certificate')
 dump(target,summary);print(json.dumps(summary,indent=2))
def search():
 target=HERE/'SUMMARY.json';assert not target.exists();start=time.time();pf=HERE/'PLAN.json'
 assert not pf.exists();dump(pf,dict(config=C,script_SHA256=sha(__file__),fit_helper_SHA256=sha(BASE),scope='Table-specific fixed mirror class labels, same count weights; true source. Exploratory floats only.'))
 rows=[]
 for n in C['dimensions']:
  nt=time.time();x=np.random.default_rng(C['train_seed']+n).random((C['train_N'],n));feat,fixed,err=features(x);best,training=old.fit(feat,fixed,C,n);_,s,lam,pi=best
  nums=[sum(math.comb(n,k)for k in range(n+1)if s[k]==b)for b in range(4)]
  frozen=dict(n=n,lambda_=lam,selector=s.tolist(),prior_mass=pi.tolist(),prior_exact_numerators=nums,prior_exact_denominator=2**n,training=training,plan_SHA256=sha(pf))
  fp=HERE/f'FROZEN_n{n}.json';assert not fp.exists();dump(fp,frozen);del x,feat,fixed
  rng=np.random.default_rng(C['test_seed']+n);counts=np.zeros(6,dtype=int);raw=counts.copy();near=counts.copy();minimum=np.full(6,np.inf)
  for at in range(0,C['test_N'],C['batch']):
   x=rng.random((min(C['batch'],C['test_N']-at),n));feat,fixed,e=features(x);err=max(err,e)
   mx=np.concatenate([old.value(feat,s).max(-1)[None,:],np.einsum('b,bnj->nj',pi,fixed).max(-1)[None,:],fixed.max(-1)])
   counts+=(mx>lam+C['guard']).sum(-1);raw+=(mx>lam).sum(-1);near+=(abs(mx-lam)<=C['guard']).sum(-1);minimum=np.minimum(minimum,abs(mx-lam).min(-1))
  values=(lam*counts/C['test_N']).tolist()
  row=dict(n=n,lambda_=lam,N=C['test_N'],model_order=['gated','actual_prior','fixed0','fixed1','fixed2','fixed3'],counts=counts.tolist(),raw_counts=raw.tolist(),near_counts=near.tolist(),minimum_margin=minimum.tolist(),weak_float_estimates=values,
    gate_minus_bestfixed=values[0]-max(values[2:]),gate_minus_prior=values[0]-values[1],DP_error=err,prior_mass=pi.tolist(),frozen_SHA256=sha(fp),seconds=time.time()-nt)
  rows.append(row);print(json.dumps(row),flush=True)
 dump(target,dict(cases=rows,plan_SHA256=sha(pf),seconds=time.time()-start,full_FP_or_continuum_certificate=False,Cn_log_gain_proved=False))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['scan','search']);a=ap.parse_args();scan()if a.mode=='scan'else search()
