import sys
sys.dont_write_bytecode=True
import probe as p
import is_engine as eng
from fractions import Fraction as F
import numpy as np,json,math
H=p.H;C=p.C
assert p.sha(H/'probe.py')==json.loads((H/'PLAN.json').read_text())['script_SHA256']
# Exact source masses and rational 1d interval integrals, including long windows.
err=0.;checks=0
for n in C['dimensions']:
 for b in C['bank']:
  m=b['m'];delta=F(str(b['width_factor']))/n;eps=F(1,50)
  assert eps+(1-eps)*m*delta/(m*delta)==1
  for x in [F(1,16),F(7,16),F(15,16)]:
   for LL in [F(2,5),F(1),F(13,5),F(217,10)]:
    integral=eps*LL
    for j in range(m):
     left=(F(str(b['phase']))+F(j,m))%1
     integral+=(1-eps)/(m*delta)*sum(max(F(0),min(x+LL/2,left+delta+k)-max(x-LL/2,left+k))for k in range(-20,21))
    exact=integral/LL;fl=float(p.kernel(np.array([float(x)]),float(LL),b,n)[0]);err=max(err,abs(float(exact)-fl));checks+=1
assert err<1e-11
class Adapter:
 def sample_1d(self,rng,shape,bank_id,n):
  bank=C['bank'][bank_id];m=bank['m'];delta=bank['width_factor']/n;bg=rng.random(shape)<C['epsilon'];peak=(bank['phase']+rng.integers(m,size=shape)/m+delta*rng.random(shape))%1;return np.where(bg,rng.random(shape),peak)
 def average_1d(self,x,L,bank_id,n):
  return np.ones_like(x)if float(L).is_integer()else p.kernel(x,L,C['bank'][bank_id],n)
a=Adapter();ee=eng.WholeVectorIS(a,range(4),32,p.lengths(32));x,_=ee.draw(np.random.default_rng(450150),128);scores=ee.scores(x);scoreerror=float(np.max(np.abs(scores-p.log_scores(x,32,p.lengths(32)))));qerror=float(np.max(np.abs(ee.log_density(scores)-p.log_q(scores))));assert scoreerror<1e-11 and qerror<1e-11
# Independently evaluate complete mixture density in longdouble on a random n2 point set.
xx=np.random.default_rng(450151).random((100,2));s=p.log_scores(xx,32,p.lengths(32));q_direct=np.longdouble(C['uniform_proposal'])
for b,bank in enumerate(C['bank']):
 for L in p.lengths(32):q_direct=q_direct+np.longdouble((1-C['uniform_proposal'])/(4*len(p.lengths(32))))*np.prod(p.kernel(xx,L,bank,32).astype(np.longdouble),axis=1)
qerr=float(np.max(np.abs(np.log(q_direct)-p.log_q(s))));assert qerr<1e-12
plateau=json.loads((H/'PLATEAU_REPLAY.json').read_text());assert all(r['guarded_classifications_unchanged']and r['max_weak_change']<1e-12 for r in plateau['results'])
summary=json.loads((H/'SUMMARY.json').read_text());allcaps=True
for row in summary['results']:
 for rr in row['rows']:
  assert max(rr['max_integrand'])<=rr['integrand_cap']*(1+1e-9)
 assert row['RN_max']<=1/C['uniform_proposal']*(1+1e-10)
 assert sum(row['proposal_component_hist'])==C['test_N']
p.dump('VERIFICATION.json',dict(passed=True,rational_kernel_cases=checks,max_kernel_float_error=err,generic_engine_score_maxerror=scoreerror,generic_engine_logq_maxerror=qerror,independent_longdouble_complete_q_maxerror=qerr,all_observed_integrand_caps_passed=True,analytic_plateau_replay_guarded_counts_unchanged=True,RN_cap_passed=True,scope='Actual algebra and numerical crosschecks; no full interval/continuous query certificate. Highdim holdout aggregates not independently rerun.'))
print(json.dumps(dict(passed=True,rational_cases=checks,kernel_error=err,engine_qerror=qerror,direct_qerror=qerr)))
