"""Query reweighting of true periodic-box class posteriors. Auxiliary, not C_n."""
from pathlib import Path
import importlib.util,json,hashlib,time,sys
import numpy as np
HERE=Path(__file__).resolve().parent
BASE=HERE.parent/'round042_channel/probe.py'
if not BASE.exists():BASE=Path(__file__).resolve().parents[2]/'042/round042_channel/probe.py'
spec=importlib.util.spec_from_file_location('channel42',BASE);base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
CONFIG=dict(K=[4,8,16,32],q=[1,2,3],ell=[.1,1/3,.5,2/3,.9],kinds=['uniform','positive_box'],shifts=[0,.125,.5],
 score_tilts=[-32,-8,-2,-1,0,1,2,8,32],random_measures_per_channel=8,point_pairs_per_channel=8,
 seed=202610290430,quadrature_order=16)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def posterior(K,q,ell,kind,shift,order=16):
 c=base.source(K,kind,shift);x,w=base.nodes(c,ell,order)
 r=base.integral(c,K,x+ell/2)-base.integral(c,K,x-ell/2)
 a=(q+r.sum(axis=1))/(q+ell);p=(q/K+r)/((q+ell)*a[:,None]);return x,w,a,p

def trace(p,w):
 w=w/w.sum();eta=w@p;v=w@(p*p)
 assert eta.min()>0 and abs(eta.sum()-1)<3e-13
 return float(np.sum(v/eta)-1),eta

def main():
 target=HERE/'SUMMARY.json';assert not target.exists();start=time.time()
 plan=dict(config=CONFIG,script_SHA256=sha(__file__),channel42_SHA256=sha(BASE),scope='Pre-frozen numerical exploration of arbitrary query laws and score tilts. Not an asymptotic or C_n estimate.')
 pf=HERE/'PLAN.json';assert not pf.exists();pf.write_text(json.dumps(plan,indent=2)+'\n')
 rng=np.random.default_rng(CONFIG['seed']);groups={};worst=[];refine=[];count=0;max_violation=-np.inf
 for K in CONFIG['K']:
  for q in CONFIG['q']:
   for ell in CONFIG['ell']:
    for kind in CONFIG['kinds']:
     for shift in CONFIG['shifts']:
      x,w,a,p=posterior(K,q,ell,kind,shift);spec=dict(K=K,q=q,ell=ell,kind=kind,shift=shift)
      cases=[]
      for t in CONFIG['score_tilts']:
       logs=np.log(w)+(1+t)*np.log(a);ww=np.exp(logs-logs.max());cases.append(('score_tilt',t,ww))
      for j in range(CONFIG['random_measures_per_channel']):
       logs=8*rng.standard_normal(len(w));ww=np.exp(logs-logs.max());cases.append(('random_query',j,ww))
      for j in range(CONFIG['point_pairs_per_channel']):
       idx=rng.choice(len(w),2,replace=False);ww=np.zeros(len(w));z=float(rng.uniform(.01,.99));ww[idx]=[z,1-z];cases.append(('two_points',j,ww))
      local=[]
      for name,key,ww in cases:
       value,eta=trace(p,ww);count+=1;group=groups.setdefault(name,dict(count=0,max_trace=-1.,min_trace=1.,max_case=None))
       group['count']+=1;group['min_trace']=min(group['min_trace'],value)
       if value>group['max_trace']:group.update(max_trace=value,max_case=dict(**spec,parameter=key,eta=eta.tolist()))
       local.append((value,name,key));max_violation=max(max_violation,value-1/(q*(q+1)))
      value,name,key=max(local);worst.append(dict(**spec,trace=value,measure=name,parameter=key))
      if q==1 and K in [4,32]and kind=='positive_box'and shift==0 and ell in [1/3,.5]:
       records=[]
       for order in [8,16,32,64]:
        xx,ww,aa,pp=posterior(K,q,ell,kind,shift,order)
        for t in [-32,0,32]:
         logs=np.log(ww)+(1+t)*np.log(aa);weight=np.exp(logs-logs.max());val,eta=trace(pp,weight)
         records.append(dict(order=order,t=t,trace=val))
       refine.append(dict(**spec,records=records))
 out=dict(config=CONFIG,plan_SHA256=sha(pf),measures_evaluated=count,groups=groups,
    worst_cases=sorted(worst,key=lambda r:r['trace'],reverse=True)[:8],
    maximum_excess_over_q_Kantorovich_candidate=float(max_violation),refinement=refine,
    python=sys.version,numpy=np.__version__,seconds=time.time()-start,
    proof_status='Numerical evidence precedes formal proof work; no C_n log gain or continuous-parameter numeric certificate')
 target.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
