"""Independent NumPy-only frozen-test audit. Writes only REPLAY_SUMMARY.json.
Reconstructs source convolution by full-period plus residual-arc overlap,
not the production primitive-difference implementation; no training rerun.
"""
from pathlib import Path
import json,time,math,hashlib
import numpy as np
H=Path(__file__).resolve().parent

def read(name):return json.loads((H/name).read_text())
def digest(name):return hashlib.sha256((H/name).read_bytes()).hexdigest()
def lse(a,axis):
 m=np.max(a,axis=axis,keepdims=True)
 return np.squeeze(m+np.log(np.exp(a-m).sum(axis=axis,keepdims=True)),axis=axis)
def draw(rng,N,n,ls,C):
 R=len(C['bank']);J=len(ls);gamma=C['uniform_proposal']
 component=rng.choice(1+R*J,N,p=[gamma]+[(1-gamma)/(R*J)]*(R*J));xx=rng.random((N,n))
 for code in np.unique(component):
  if code==0:continue
  ids=np.flatnonzero(component==code);b,j=divmod(int(code)-1,J);bank=C['bank'][b];shape=(len(ids),n)
  background=rng.random(shape)<C['epsilon']
  peak=(bank['phase']+rng.integers(bank['m'],size=shape)/bank['m']+(bank['width_factor']/n)*rng.random(shape))%1
  yy=np.where(background,rng.random(shape),peak)
  xx[ids]=(yy+ls[j]*(rng.random(shape)-.5))%1
 return xx,component

def score(xx,n,ls,C):
 out=np.empty((len(xx),len(C['bank']),len(ls)))
 for b,bank in enumerate(C['bank']):
  m=bank['m'];width=bank['width_factor']/n
  for j,L in enumerate(ls):
   whole=math.floor(L);tail=L-whole
   if tail==0:out[:,b,j]=0.;continue
   covered=np.zeros_like(xx)
   for peak in range(m):
    left=(bank['phase']+peak/m)%1
    start=np.remainder(xx-L/2-left,1.)
    # Residual arc [start,start+tail], with at most one wrap.
    covered+=whole*width+np.minimum(tail,np.maximum(width-start,0.))+np.minimum(width,np.maximum(start+tail-1.,0.))
   aa=C['epsilon']+(1-C['epsilon'])*covered/(m*width*L)
   assert np.min(aa)>0
   out[:,b,j]=np.log(aa).sum(1)
 return out

def moments(N,ss,ss2,caps,alpha):
 mean=ss/N;var=np.maximum(0.,(ss2-ss*ss/N)/(N-1));logterm=math.log(2/alpha)
 radius=np.sqrt(2*var*logterm/N)+7*caps[:,None]*logterm/(3*(N-1))
 ci=np.stack([np.maximum(0.,mean-radius),np.minimum(caps[:,None],mean+radius)],axis=-1)
 return mean,var,radius,ci

def run():
 start=time.time();plan=read('PLAN.json');original=read('SUMMARY.json');C=plan['config'];assert C==original['config']
 alpha=original['alpha_per_direction'];directions=2*len(C['dimensions'])*(5*len(C['log_thresholds'])+1);assert directions==original['directions']==336;assert alpha==C['alpha']/directions
 if (H/'probe.py').exists():assert digest('probe.py')==plan['script_SHA256']
 all_results=[];hs=np.array(C['log_thresholds']);N=C['test_N'];R=len(C['bank'])
 for n in C['dimensions']:
  st=time.time();freeze=read(f'FROZEN_n{n}.json');assert freeze['plan_SHA256']==digest('PLAN.json')
  train=np.array(freeze['training_IS_grid']);mi,hi=np.unravel_index(np.argmax(train),train.shape)
  assert freeze['pi']==C['mixture_candidates'][mi]and freeze['log_threshold']==C['log_thresholds'][hi]
  old=next(r for r in original['results']if r['n']==n);ls=plan['lengths_by_dimension'][str(n)];J=len(ls);gamma=C['uniform_proposal'];v=(1-gamma)/(R*J)
  models=np.vstack([freeze['pi'],np.eye(R)]);caps=models.max(1)/v
  with np.errstate(divide='ignore'):logpi=np.log(models)
  sums=np.zeros((5,len(hs)));sums2=sums.copy();correct_sums=sums.copy();correct_sums2=sums.copy();hits=np.zeros_like(sums,int);raw=hits.copy();near=hits.copy();correct_hits=hits.copy();correct_raw=hits.copy();mind=np.full_like(sums,np.inf);maxz=np.zeros_like(sums)
  nearids={};hist=np.zeros(1+R*J,int);rng=np.random.default_rng(C['test_seed']+n);sumrn=sumrn2=maxrn=0.;max_identity_score=max_identity_logq=0.
  for at in range(0,N,C['batch']):
   xx,component=draw(rng,min(C['batch'],N-at),n,ls,C);hist+=np.bincount(component,minlength=len(hist));scores=score(xx,n,ls,C)
   # Complete whole-vector mixture; no coordinate-wise mixing or selected-component denominator.
   logq=np.logaddexp(math.log(gamma),math.log(v)+lse(scores,axis=(1,2)))
   corrected=scores.copy()
   for b,bank in enumerate(C['bank']):
    for j,L in enumerate(ls):
     if float(L*bank['m']).is_integer():corrected[:,b,j]=0.
   corrected_logq=np.logaddexp(math.log(gamma),math.log(v)+lse(corrected,axis=(1,2)))
   max_identity_score=max(max_identity_score,float(np.max(np.abs(corrected-scores))));max_identity_logq=max(max_identity_logq,float(np.max(np.abs(corrected_logq-logq))))
   rn=np.exp(-logq);sumrn+=rn.sum();sumrn2+=(rn*rn).sum();maxrn=max(maxrn,float(rn.max()))
   base=np.exp(hs[None,:]-logq[:,None]);correct_base=np.exp(hs[None,:]-corrected_logq[:,None])
   for m in range(5):
    value=lse(scores+logpi[m][None,:,None],axis=1).max(1);difference=value[:,None]-hs;guarded=difference>C['guard_log'];z=base*guarded
    sums[m]+=z.sum(0);sums2[m]+=(z*z).sum(0);hits[m]+=guarded.sum(0);raw[m]+=(difference>0).sum(0);near[m]+=(np.abs(difference)<=C['guard_log']).sum(0);mind[m]=np.minimum(mind[m],np.abs(difference).min(0));maxz[m]=np.maximum(maxz[m],z.max(0))
    for qi,hj in np.argwhere(np.abs(difference)<=C['guard_log']):nearids.setdefault(f'{m}_{hj}',[]).append(int(at+qi))
    cv=lse(corrected+logpi[m][None,:,None],axis=1).max(1);cd=cv[:,None]-hs;cy=cd>C['guard_log'];cz=correct_base*cy
    correct_sums[m]+=cz.sum(0);correct_sums2[m]+=(cz*cz).sum(0);correct_hits[m]+=cy.sum(0);correct_raw[m]+=(cd>0).sum(0)
  mean,var,radius,ci=moments(N,sums,sums2,caps,alpha);cmean,cvar,cradius,cci=moments(N,correct_sums,correct_sums2,caps,alpha)
  expected=np.array([r['IS_weak']for r in old['rows']]);expected_var=np.array([r['sample_variance']for r in old['rows']]);expected_ci=np.array([r['EB_CI']for r in old['rows']]);expected_sums=N*expected;expected_sums2=(N-1)*expected_var+N*expected*expected
  expected_radius=np.sqrt(2*expected_var*math.log(2/alpha)/N)+7*caps[:,None]*math.log(2/alpha)/(3*(N-1))
  nearname=f'NEAR_INDICES_n{n}.json';near_frozen=read('NEAR_SHA256.json')['files'][nearname]
  canonical=(json.dumps(dict(indices=nearids,scope='All near query indices, regenerated from frozen seed; no query/features cache.'),indent=2)+'\n').encode('utf-8')
  rebuilt_near_sha=hashlib.sha256(canonical).hexdigest();near_sha_matches=(rebuilt_near_sha==near_frozen['sha256']and len(canonical)==near_frozen['canonical_bytes'])
  assert near_sha_matches,'Regenerated near canonical hash differs from frozen original'
  if (H/nearname).exists():
   assert digest(nearname)==near_frozen['sha256']
   near_exact=(nearids==read(nearname)['indices']);near_mode='array_and_canonical_sha'
  else:near_exact=near_sha_matches;near_mode='canonical_sha_only'
  counts_exact=(hits.tolist()==[r['hits']for r in old['rows']]);near_exact=near_exact and near.tolist()==[r['near']for r in old['rows']]
  deltas=dict(mean_max_absolute=float(np.max(np.abs(mean-expected))),sums_max_absolute=float(np.max(np.abs(sums-expected_sums))),sums2_max_absolute=float(np.max(np.abs(sums2-expected_sums2))),variance_max_absolute=float(np.max(np.abs(var-expected_var))),EB_radius_max_absolute=float(np.max(np.abs(radius-expected_radius))),EB_endpoint_max_absolute=float(np.max(np.abs(ci-expected_ci))),raw_count_max_difference=int(np.max(np.abs(raw-np.array([r['raw_hits']for r in old['rows']])))))
  assert counts_exact and near_exact and deltas['mean_max_absolute']<1e-9 and deltas['EB_endpoint_max_absolute']<1e-9 and deltas['sums2_max_absolute']<1e-4
  assert np.array_equal(correct_hits,hits)and np.max(maxz-caps[:,None])<=1e-8 and maxrn<=1/gamma*(1+1e-12)
  assert hist.tolist()==old['proposal_component_hist']
  profile=[]
  for m,rr in enumerate(old['rows']):
   j=int(mean[m].argmax());jj=int(expected[m].argmax());profile.append(dict(model=rr['model'],max_IS_weak=float(mean[m,j]),max_log_threshold=float(hs[j]),saved_argmax_matches=j==jj))
  fixed_flat=int(mean[1:].argmax());fb,fh=np.unravel_index(fixed_flat,mean[1:].shape);mj=int(mean[0].argmax())
  rnmean=sumrn/N;rnvar=(sumrn2-sumrn*sumrn/N)/(N-1);rnrad=math.sqrt(2*rnvar*math.log(2/alpha)/N)+7*(1/gamma)*math.log(2/alpha)/(3*(N-1))
  rnci=[max(0.,rnmean-rnrad),min(1/gamma,rnmean+rnrad)]
  all_results.append(dict(n=n,seconds=time.time()-st,selected_pi=freeze['pi'],selected_h=freeze['log_threshold'],training_choice_matches_frozen=True,guarded_counts_exact=counts_exact,near_ids_exact=near_exact,near_verification_mode=near_mode,near_canonical_bytes=len(canonical),near_rebuilt_SHA256=rebuilt_near_sha,near_total=sum(map(len,nearids.values())),comparison_errors=deltas,
   sums=sums.tolist(),sums2=sums2.tolist(),mean=mean.tolist(),sample_variance=var.tolist(),EB_radius=radius.tolist(),EB_CI=ci.tolist(),per_model_full_profile_max=profile,
   whole_profile_mixture_max=float(mean[0,mj]),whole_profile_mixture_h=float(hs[mj]),whole_profile_bestfixed_max=float(mean[1+fb,fh]),whole_profile_bestfixed_model=f'fixed{fb}',whole_profile_bestfixed_h=float(hs[fh]),mixture_fullmax_minus_bestfixed_fullmax=float(mean[0,mj]-mean[1+fb,fh]),
   plateau_audit=dict(corrected_guarded_counts_unchanged=True,max_score_identity_change=max_identity_score,max_logq_identity_change=max_identity_logq,corrected_mean=cmean.tolist(),corrected_sums=correct_sums.tolist(),corrected_sums2=correct_sums2.tolist(),corrected_EB_radius=cradius.tolist(),corrected_raw_counts=correct_raw.tolist(),max_mean_change=float(np.max(np.abs(cmean-mean)))),
   RN_audit=dict(mean=rnmean,sums=sumrn,sums2=sumrn2,max=maxrn,EB_radius=rnrad,EB_CI=rnci,mean_vs_saved_absolute=abs(rnmean-old['RN_uniform_mass_mean']),CI_vs_saved_max_absolute=float(np.max(np.abs(np.array(rnci)-old['RN_uniform_mass_EB_CI'])))),
   raw_count_scope='Raw h0 may differ across independent algebra due exact plateaus; guarded counts and all near IDs must match.'))
 out=dict(passed=True,seconds=time.time()-start,numpy=np.__version__,dimensions=C['dimensions'],all_thresholds=C['log_thresholds'],results=all_results,alpha_per_direction=alpha,directions=directions,implementation='Independent NumPy-only RNG replay, full-period plus residual overlap kernel, local stable logsumexp, complete whole-vector q.',scope='Floating reproducibility audit, not a strict probability/interval/continuous-cell certificate.',portable_minimum_files=['replay_summary.py','PLAN.json','SUMMARY.json','NEAR_SHA256.json']+[f'FROZEN_n{n}.json'for n in C['dimensions']],input_SHA256={name:digest(name)for name in ['PLAN.json','SUMMARY.json','NEAR_SHA256.json']+[f'FROZEN_n{n}.json'for n in C['dimensions']]})
 (H/'REPLAY_SUMMARY.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(passed=True,seconds=round(out['seconds'],3),dimensions=C['dimensions'],near_modes=[r['near_verification_mode']for r in all_results],canonical_hashes=len(all_results),canonical_bytes=sum(r['near_canonical_bytes']for r in all_results),max_mean_error=max(r['comparison_errors']['mean_max_absolute']for r in all_results))))
if __name__=='__main__':run()
