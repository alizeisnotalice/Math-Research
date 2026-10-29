"""Same frozen samples, analytic half-period identities; no retraining or reselection."""
import sys
sys.dont_write_bytecode=True
import probe as p
import numpy as np,json,math
H=p.H;C=p.C;results=[]
for n in C['dimensions']:
 f=json.loads((H/f'FROZEN_n{n}.json').read_text());old=json.loads((H/f'RESULT_n{n}.json').read_text());ls=p.lengths(n);models=[f['pi']]+np.eye(4).tolist();hs=np.array(C['log_thresholds']);rng=np.random.default_rng(C['test_seed']+n);ss=np.zeros((5,len(hs)));counts=np.zeros_like(ss,dtype=int);raw=counts.copy();maxdiff=0.;qdiff=0.
 for at in range(0,C['test_N'],C['batch']):
  x,_=p.sample(rng,min(C['batch'],C['test_N']-at),n,ls);s=p.log_scores(x,n,ls);s2=s.copy()
  for b,bank in enumerate(C['bank']):
   for j,L in enumerate(ls):
    if float(L*bank['m']).is_integer():s2[:,b,j]=0.
  maxdiff=max(maxdiff,float(np.max(np.abs(s2-s))));lq=p.log_q(s2);qdiff=max(qdiff,float(np.max(np.abs(lq-p.log_q(s)))))
  for m,pi in enumerate(models):
   v=p.log_model(s2,pi);d=v[:,None]-hs[None,:];yes=d>C['guard_log'];ss[m]+=(np.exp(hs[None,:]-lq[:,None])*yes).sum(0);counts[m]+=yes.sum(0);raw[m]+=(d>0).sum(0)
 assert counts.tolist()==[r['hits']for r in old['rows']]
 results.append(dict(n=n,guarded_classifications_unchanged=True,corrected_raw_counts=raw.tolist(),corrected_IS_weak=(ss/C['test_N']).tolist(),max_logscore_identity_correction=maxdiff,max_logq_change=qdiff,max_weak_change=float(np.max(np.abs(ss/C['test_N']-np.array([r['IS_weak']for r in old['rows']])))),scope='Frozen-only same-seed identity audit. No interval or query-cell certificate.'))
p.dump('PLATEAU_REPLAY.json',dict(results=results,identity='rho_b has period 1/m; if m*L is integer then a_b,L=1 exactly.'))
print(json.dumps([dict(n=r['n'],same=r['guarded_classifications_unchanged'],weak_change=r['max_weak_change'])for r in results]))
