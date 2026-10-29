"""MILP pilot on frozen round043 source; true nonnegative count DP, no Gaussian proxy."""
from pathlib import Path
import importlib.util,json,time,math,hashlib,argparse
import sys
sys.dont_write_bytecode=True
import numpy as np
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix
from scipy.stats import beta
H=Path(__file__).resolve().parent
SOURCE=H.parent/'round043_binary'/'probe.py'
spec=importlib.util.spec_from_file_location('round043_source',SOURCE);source=importlib.util.module_from_spec(spec);spec.loader.exec_module(source)
C=source.CFG.copy();C.update(n=32,train_N=128,test_N=8192,train_seed=202610294410,test_seed=202610294411,thresholds=[1.5,2.,3.],milp_margin=1e-7,guard=1e-9,time_limit=30.)
def dump(name,obj):(H/name).write_text(json.dumps(obj,indent=2)+'\n')
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def solve(feat,lam,limit):
 R,N,J,T=feat.shape;W=R*T;Y=N*J;target=lam+C['milp_margin'];coeff=np.minimum(1.,feat/target).transpose(1,2,0,3).reshape(Y,W)
 # Exact conceptual row: sum capped feature*w - y>=0. All coefficients nonnegative.
 r=np.concatenate([np.repeat(np.arange(Y),W),np.arange(Y),Y+np.repeat(np.arange(T),R),Y+T+np.repeat(np.arange(N),J)])
 col=np.concatenate([np.tile(np.arange(W),Y),W+np.arange(Y),np.tile(np.arange(R),T)*T+np.repeat(np.arange(T),R),W+np.arange(Y)])
 data=np.concatenate([coeff.ravel(),-np.ones(Y),np.ones(W),np.ones(Y)])
 A=coo_matrix((data,(r,col)),shape=(Y+T+N,W+Y)).tocsc();lo=np.concatenate([np.zeros(Y),np.ones(T),np.full(N,-np.inf)]);hi=np.concatenate([np.full(Y,np.inf),np.ones(T),np.ones(N)])
 objective=np.concatenate([np.zeros(W),-np.ones(Y)]);start=time.time()
 res=milp(objective,integrality=np.ones(W+Y),bounds=Bounds(np.zeros(W+Y),np.ones(W+Y)),constraints=LinearConstraint(A,lo,hi),options=dict(time_limit=limit,mip_rel_gap=0.,presolve=True))
 rec=dict(lambda_=lam,seconds=time.time()-start,status=int(res.status),message=str(res.message),success=bool(res.success),time_limit=limit,variables=W+Y,constraints=Y+T+N,nonzeros=int(A.nnz),incumbent_available=res.x is not None,
   mip_gap=None if getattr(res,'mip_gap',None) is None else float(res.mip_gap),mip_node_count=None if getattr(res,'mip_node_count',None) is None else int(res.mip_node_count),mip_dual_bound=None if getattr(res,'mip_dual_bound',None) is None else float(res.mip_dual_bound),objective=None if res.fun is None else float(res.fun))
 if res.x is None:return rec,None
 w=res.x[:W].reshape(R,T);s=w.argmax(0);rounded=np.eye(R)[s].T;v=source.selector_value(feat,s) if hasattr(source,'selector_value') else source.value(feat,s);mx=v.max(-1)
 y=res.x[W:].reshape(N,J);claimed=y.max(-1)>.5
 # Recompute coefficient rows and original, uncapped real features using rounded map.
 gate_rows=np.einsum('inrt,rt->in',coeff.reshape(N,J,R,T),rounded)
 claimed_rows=y>.5
 rec.update(selector=s.tolist(),integrality_max_error=float(np.max(np.abs(res.x-np.rint(res.x)))),rounded_onehot_max_error=float(np.max(np.abs(w-rounded))),claimed_queries=int(claimed.sum()),
  recomputed_train_successes=int((mx>lam+C['guard']).sum()),recomputed_margin_successes=int((mx>=target).sum()),recomputed_train_weak=float(lam*np.mean(mx>lam+C['guard'])),
  claimed_query_fail_at_lambda_guard=int(np.sum(claimed & (mx<=lam+C['guard']))),claimed_rows_below_target=int(np.sum(claimed_rows & (v<target))),minimum_claimed_original_slack=None if not claimed_rows.any()else float(np.min(v[claimed_rows]-target)),minimum_claimed_cap_slack=None if not claimed_rows.any()else float(np.min(gate_rows[claimed_rows]-1)),train_near_count=int((np.abs(mx-lam)<=C['guard']).sum()),
  minimum_train_distance=float(np.min(np.abs(mx-lam))))
 return rec,s

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--cost',action='store_true');args=ap.parse_args()
 if args.cost:
  x=np.random.default_rng(4400).random((64,32));ft,_,err=source.features(x,C);r,s=solve(ft,2.,3.);r['DP_partition_error']=err;dump('COST.json',r);print(json.dumps(r));return
 assert not(H/'SUMMARY.json').exists()
 dump('PLAN.json',dict(config=C,script_SHA256=sha(__file__),source_path=str(SOURCE),source_SHA256=sha(SOURCE),cost_SHA256=sha(H/'COST.json'),objective='Maximize train query witnesses over shared onehot count selector. Select frozen lambda and rounded selector by recomputed guarded train weak, not solver claimed objective.',test_scope='Ideal iid uniform floating point exploration; no outward intervals or query-cell certificate.',formal_train_after_cost=True))
 start=time.time();rng=np.random.default_rng(C['train_seed']);x=rng.random((C['train_N'],32));feat,fixed,err=source.features(x,C);records=[];best=None
 for lam in C['thresholds']:
  r,s=solve(feat,lam,C['time_limit']);records.append(r);dump('TRAINING.json',records);print(json.dumps(r),flush=True)
  if s is not None and (best is None or r['recomputed_train_weak']>best[0]):best=(r['recomputed_train_weak'],s.copy(),lam)
 assert best is not None,'No incumbent; training status retained, need explicit fallback'
 _,s,lam=best;prior_num=[sum(math.comb(32,k)for k in range(33)if s[k]==b)for b in range(4)];pi=np.array(prior_num)/2**32
 frozen=dict(n=32,lambda_=lam,selector=s.tolist(),prior_mass=pi.tolist(),prior_exact_numerators=prior_num,prior_exact_denominator=2**32,training_SHA256=sha(H/'TRAINING.json'),plan_SHA256=sha(H/'PLAN.json'),selection='Largest recomputed train weak at lambda+guard, ties first preset lambda')
 dump('FROZEN.json',frozen);del feat,fixed,x
 rng=np.random.default_rng(C['test_seed']);counts=np.zeros(6,dtype=int);raw=counts.copy();near=counts.copy();top=np.zeros(6);mind=np.full((6,3),np.inf);allnear=np.zeros((6,3),int);allraw=np.zeros((6,3),int);allguard=np.zeros((6,3),int);examples=[]
 for at in range(0,C['test_N'],C['batch']):
  x=rng.random((min(C['batch'],C['test_N']-at),32));feat,fixed,e=source.features(x,C);err=max(err,e)
  maxima=np.concatenate([source.value(feat,s).max(-1)[None,:],np.einsum('b,bnj->nj',pi,fixed).max(-1)[None,:],fixed.max(-1)])
  counts+=(maxima>lam+C['guard']).sum(-1);raw+=(maxima>lam).sum(-1);near+=(np.abs(maxima-lam)<=C['guard']).sum(-1);top=np.maximum(top,maxima.max(-1))
  dif=maxima[:,:,None]-np.array(C['thresholds'])[None,None,:];mind=np.minimum(mind,np.abs(dif).min(1));allnear+=(np.abs(dif)<=C['guard']).sum(1);allraw+=(dif>0).sum(1);allguard+=(dif>C['guard']).sum(1)
  for mod,q,j in np.argwhere(np.abs(dif)<=C['guard']):examples.append(dict(model=int(mod),sample_index=int(at+q),threshold=C['thresholds'][j],value=float(maxima[mod,q]),x=x[q].tolist()))
 N=C['test_N'];rows=[];a=.05/12
 for i,name in enumerate(['gate','actual_prior']+[f'fixed{b}'for b in range(4)]):
  k=int(counts[i]);lo=0. if k==0 else float(beta.ppf(a,k,N-k+1));hi=1. if k==N else float(beta.ppf(1-a,k+1,N-k));rows.append(dict(model=name,k=k,N=N,weak=float(lam*k/N),CP_weak=[lam*lo,lam*hi],raw_count=int(raw[i]),near_guard=int(near[i]),max_observed=float(top[i])))
 out=dict(config=C,lambda_=lam,rows=rows,gate_minus_bestfixed=rows[0]['weak']-max(r['weak']for r in rows[2:]),gate_minus_prior=rows[0]['weak']-rows[1]['weak'],training_best=best[0],gate_train_minus_test=best[0]-rows[0]['weak'],DP_partition_error=err,guard_diagnostics=dict(thresholds=C['thresholds'],raw_counts=allraw.tolist(),guard_counts=allguard.tolist(),near_counts=allnear.tolist(),minimum_distance=mind.tolist(),all_near_examples=examples),seconds=time.time()-start,source_mass=1.,full_interval_certificate=False,Cn_growth_proved=False,CI_scope='12 CP directions; ideal iid point-query frozen guarded FP classifier only; no continuum cell or outward integral certification',frozen_SHA256=sha(H/'FROZEN.json'))
 dump('SUMMARY.json',out);print(json.dumps(dict(lambda_=lam,weak=[r['weak']for r in rows],seconds=out['seconds'])),flush=True)
if __name__=='__main__':main()
