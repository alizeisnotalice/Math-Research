"""Frozen count-space quantile-source tail pilot; zero n-length query vectors.

PB/Cantelli envelopes pay analytic approximation errors. Ordinary FP and its
empirical guard are not an interval certificate. No asymptotic C_n claim.
"""
import os
for k in ['VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS']:os.environ.setdefault(k,'4')
from pathlib import Path
import argparse,hashlib,json,math,time,sys,resource
from fractions import Fraction as F
import numpy as np
from build_cache import build_cache

HERE=Path(__file__).resolve().parent
SEED=202610290381
HGRID=np.arange(0,10.00001,.5)
MODELS=['fixed0','fixed1','fixed2','fixed3','joint2','joint4','prior2','prior4']

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def lse(x,axis=-1):
 m=np.max(x,axis=axis,keepdims=True);safe=np.where(np.isfinite(m),m,0.)
 with np.errstate(divide='ignore',invalid='ignore',under='ignore'):out=safe+np.log(np.exp(x-safe).sum(axis=axis,keepdims=True))
 return np.squeeze(out,axis)
def phi(z):
 z=np.asarray(z);a=np.fromiter((.5*math.erfc(-float(v)/math.sqrt(2))for v in z.ravel()),dtype=float,count=z.size)
 return a.reshape(z.shape)
def cdf_enclosure(cut,mean,var,error,eps):
 vm=var-error;vp=var+error
 assert np.min(vm)>0
 sm=np.sqrt(vm);sp=np.sqrt(vp)
 nm=cut+.5-(mean+error);np_=cut+.5-(mean-error)
 zl=np.minimum(nm/sm,nm/sp);zu=np.maximum(np_/sm,np_/sp)
 return np.maximum(0,phi(zl)-eps),np.minimum(1,phi(zu)+eps)
def log_prob(p):
 with np.errstate(divide='ignore'):return np.log(p)
def log_mixture(logl,logu,logr):
 lv=2*lse(logu+.5*logr,axis=1)
 return lse(logl,axis=1),lse(logu,axis=1),lv

def intervals(z,pbmean,pbvar,logr,errors,cuts,n,pi):
 # Inputs all (batch,4,J). Errors are (4,J) parameter enclosures.
 inf=errors['loginflate'][None,:,:];pu=errors['p_sum_error'][None,:,:]
 logal=z-inf;logau=z+inf
 cdfs=[cdf_enclosure(c,pbmean,pbvar,pu,1/math.sqrt(n))for c in cuts]
 prob4=[]
 for b in range(4):
  lo_l=0. if b==0 else cdfs[b-1][0][:,b,:]
  lo_u=0. if b==0 else cdfs[b-1][1][:,b,:]
  hi_l=1. if b==3 else cdfs[b][0][:,b,:]
  hi_u=1. if b==3 else cdfs[b][1][:,b,:]
  prob4.append((np.maximum(0,hi_l-lo_u),np.minimum(1,hi_u-lo_l)))
 p2=[(cdfs[1][0][:,0,:],cdfs[1][1][:,0,:]),
     (1-cdfs[1][1][:,2,:],1-cdfs[1][0][:,2,:])]
 out={}
 for b in range(4):out['fixed'+str(b)]=(logal[:,b,:],logau[:,b,:],2*logau[:,b,:]+logr[:,b,:])
 l=np.stack([logal[:,b,:]+log_prob(prob4[b][0])for b in range(4)],axis=1)
 u=np.stack([logau[:,b,:]+log_prob(prob4[b][1])for b in range(4)],axis=1)
 out['joint4']=log_mixture(l,u,logr)
 ix=[0,2]
 l=np.stack([logal[:,b,:]+log_prob(p2[t][0])for t,b in enumerate(ix)],axis=1)
 u=np.stack([logau[:,b,:]+log_prob(p2[t][1])for t,b in enumerate(ix)],axis=1)
 out['joint2']=log_mixture(l,u,logr[:,ix,:])
 for r,indices in [(2,[0,2]),(4,[0,1,2,3])]:
  l=logal[:,indices,:]+np.array([log_prob(a)for a,b in pi[r]])[None,:,None]
  u=logau[:,indices,:]+np.array([log_prob(b)for a,b in pi[r]])[None,:,None]
  out['prior'+str(r)]=log_mixture(l,u,logr[:,indices,:])
 return out,{'minimum_class_lower_probability':min(float(a.min())for a,b in prob4),
             'mean_class_enclosure_width':float(np.mean([np.mean(b-a)for a,b in prob4]))}

def tail_enclosures(logl,logu,logv,h):
 with np.errstate(over='ignore',under='ignore',divide='ignore',invalid='ignore'):
  d=np.maximum(0,-np.expm1(h-logl));r=np.exp(logv-2*logl)
  low=np.where(logl>h,d*d/(r+d*d),0.)
  t=np.maximum(0,-np.expm1(logu-h));rv=np.exp(logv-2*h)
  up=np.where(logu<h,rv/(rv+t*t),1.)
  up=np.minimum(up,np.minimum(1,np.exp(logu-h)))
 low=np.nan_to_num(low,nan=0.);up=np.nan_to_num(up,nan=1.)
 return low.max(axis=1),np.minimum(1,up.sum(axis=1))

def bounded_ci(mean,cap,n,alpha,upper):
 x=max(0,min(1,mean/cap));target=math.log(1/alpha)/n
 if (upper and x==1)or(not upper and x==0):return cap if upper else 0.
 a,b=(x,1.)if upper else(0.,x)
 for _ in range(80):
  p=(a+b)/2
  div=(x*math.log(x/p)if x else 0)+((1-x)*math.log((1-x)/(1-p))if x<1 else 0)
  if upper:
   if div>target:b=p
   else:a=p
  else:
   if div>target:a=p
   else:b=p
 return cap*((a+b)/2)
def summarize(values,cap,alpha,upper):
 n=len(values);s=float(values.sum());ss=float(np.dot(values,values))
 if s==0:ess=0.
 else:
  positive=values[values>0];ll=np.log(positive);ess=math.exp(2*float(lse(ll))-float(lse(2*ll)))
 return {'mean':float(values.mean()),'SE_descriptive':float(values.std(ddof=1)/math.sqrt(n)),
         'ESS':ess,'CI':bounded_ci(float(values.mean()),cap,n,alpha,upper),
         'max_contribution':float(values.max())}

def draw_counts(src,features,norms,J,rng,N):
 counts=np.empty((N,src.Q),dtype=np.int64);tags=[]
 for i in range(N):
  if rng.random()<.5:
   p=np.full(src.Q,1/src.Q);tags.append([-1,-1])
  else:
   b=int(rng.integers(4));j=int(rng.integers(J));p=np.exp(features[0,b*J+j]);p/=p.sum();tags.append([b,j])
  counts[i]=rng.multinomial(src.n,p)
 return counts,tags

def micro_checks():
 # Strict threshold handling and parameter uncertainty corner directions.
 l=np.array([[math.log(2.),math.log(.5)],[0.,0.]])
 low,up=tail_enclosures(l,l,np.full_like(l,-math.inf),0.)
 assert np.array_equal(low,[1.,0.])and np.array_equal(up,[1.,1.])
 for cut in [1,4,9]:
  a,b=cdf_enclosure(cut,np.array([5.]),np.array([2.]),np.array([.01]),.01)
  assert 0<=a[0]<=b[0]<=1
 # C=.5606 and p in [.3,.7] imply CDF error <=1/sqrt(n).
 assert F(2803,5000)**2*F(29,50)**2/F(21,100)<1
 return {'Cantelli_zero_variance_strict_threshold_checks':True,'BE_unit_over_sqrt_n_rational_check':True}

def main():
 target=HERE/'HIGH_DIM_RESULTS.json'
 if target.exists()or(HERE/'HIGH_DIM_REGISTRATION.json').exists():raise RuntimeError('Refusing repeated pilot')
 micro=micro_checks();started=time.time()
 plan={'seed':SEED,'candidate_N':[128,64],'batch':16,'N_choice':'pure matrix cost:128 if predicted<=120sec else64 if<=180sec else abort',
       'peak_RSS_limit_bytes':2_000_000_000,'log_thresholds':HGRID.tolist(),'models':MODELS,
       'mean_method':'true positive-polynomial conditional mean, PB BE enclosures; parameter error corners',
       'floating_qualification':'analytic midpoint error plus empirical FP guard only; no interval numerical certificate',
       'script_SHA256':sha(Path(__file__)),'source_SHA256':sha(HERE/'source_tables.py'),'builder_SHA256':sha(HERE/'build_cache.py')}
 (HERE/'HIGH_DIM_PLAN.json').write_text(json.dumps(plan,indent=2)+'\n')
 src,grid,features,norms,errors=build_cache();J=len(grid)
 assert .49<src.pchi<.51
 errors['meta'].update({'alphabet_B':4,'B_ge_v':False,'v_ge_16logn':bool(src.v>=16*math.log(src.n)),'original_matching_asymptotic_source_claimed':False})
 # Deterministic cost timing, no RNG/observed-tail choice.
 counts=np.full((8,src.Q),src.n//src.Q,dtype=float)
 counts[:,:src.n%src.Q]+=1
 t=time.time();_cost=counts@features.reshape(4*4*J,src.Q).T;bench=time.time()-t;del counts,_cost
 estimate=bench*16*1.5
 N=128 if estimate<=120 else 64 if estimate/2<=180 else 0
 if N==0:raise RuntimeError('Matrix cost gate failed before random sampling')
 if errors['meta']['ru_maxrss_after_build']['platform']=='darwin'and errors['meta']['ru_maxrss_after_build']['raw']>2_000_000_000:
  raise RuntimeError('RSS resource gate failed before sampling')
 mu0=src.n*src.pchi;sd0=math.sqrt(src.n*src.pchi*(1-src.pchi));cuts=[int(round(mu0+a*sd0))for a in [-1,0,1]]
 c0=[cdf_enclosure(c,np.array([mu0]),np.array([sd0*sd0]),np.array([src.n*1e-13]),1/math.sqrt(src.n))for c in cuts]
 pi4=[]
 for b in range(4):
  ll=0 if b==0 else c0[b-1][0][0];lu=0 if b==0 else c0[b-1][1][0]
  hl=1 if b==3 else c0[b][0][0];hu=1 if b==3 else c0[b][1][0]
  pi4.append((max(0,hl-lu),min(1,hu-ll)))
 pi={4:pi4,2:[(c0[1][0][0],c0[1][1][0]),(1-c0[1][1][0],1-c0[1][0][0])]}
 family=2*len(MODELS)*len(HGRID);alpha=.05/family
 cap_all=2*4*J*errors['C_Z']*errors['C_E']
 reg={'plan_SHA256':sha(HERE/'HIGH_DIM_PLAN.json'),'source':errors['meta'],'N':N,'batch':16,'J':J,
      'cost_benchmark_8_rows_seconds':bench,'predicted_128_rows_matrix_seconds':estimate,
      'cuts':cuts,'cuts_definition':'explicit fixed integer cuts chosen at nominal mu+[-1,0,1]sigma; not exact quantiles',
      'prior_probability_enclosures':pi,'CDF_BE_bound':1/math.sqrt(src.n),
      'normalizer_C_Z':errors['C_Z'],'kernel_mean_C_E':errors['C_E'],'all_threshold_contribution_cap':cap_all,
      'max_logmean_inflation':float(errors['loginflate'].max()),'max_PB_parameter_sum_error':float(errors['p_sum_error'].max()),
      'max_analytic_mean_error_abs':float(errors['analytic_mean_error_abs'].max()),'empirical_FP_guard_abs':float(errors['floating_point_guard_abs'].max()),
      'alpha_total':.05,'alpha_each':alpha,'CI_family_count':family,'micro_checks':micro,
      'python':sys.version,'numpy':np.__version__,'status':'Frozen cache/cuts/N/thresholds before any random count sampling'}
 (HERE/'HIGH_DIM_REGISTRATION.json').write_text(json.dumps(reg,indent=2)+'\n')
 print('REGISTRATION',N,'cache_seconds',errors['meta']['build_runtime_seconds'],'bench8',bench,flush=True)
 rng=np.random.default_rng(SEED);lower={m:[]for m in MODELS};upper={m:[]for m in MODELS};weights=[];sample_tags=[];datahash=hashlib.sha256();rb=[];extra=[]
 sample_start=time.time()
 for first in range(0,N,16):
  count,tags=draw_counts(src,features,norms,J,rng,min(16,N-first));sample_tags.extend(tags)
  datahash.update(count.tobytes());a=count@features.reshape(4*4*J,src.Q).T
  a=a.reshape(len(count),4,4,J);z=a[:,0];pbmu=a[:,1];pbvar=a[:,2];logr=log_prob(np.expm1(a[:,3]));del count,a
  logD=np.logaddexp(math.log(.5),math.log(.5)-math.log(4*J)+lse((z-src.n*np.log(norms)[None,:,:]).reshape(len(z),-1),axis=1))
  weight=np.exp(-logD);assert weight.max()<=2+1e-12;weights.extend(weight.tolist())
  bounds,diag=intervals(z,pbmu,pbvar,logr,errors,cuts,src.n,pi);extra.append(diag)
  rb.extend(np.exp(logr).ravel().tolist())
  for model in MODELS:
   lo=[];up=[]
   for h in HGRID:
    B,U=tail_enclosures(*bounds[model],float(h));x=np.exp(h-logD)*B;y=np.exp(h-logD)*U
    cap=min(2*math.exp(h),cap_all);assert max(x.max(),y.max())<=cap*(1+1e-8)
    lo.append(x);up.append(y)
   lower[model].append(np.array(lo).T);upper[model].append(np.array(up).T)
  print('sampled',first+len(z),flush=True)
 models={}
 for model in MODELS:
  lo=np.concatenate(lower[model]);up=np.concatenate(upper[model]);rows=[]
  for j,h in enumerate(HGRID):
   cap=min(2*math.exp(h),cap_all)
   rows.append({'h':float(h),'lambda':math.exp(h),'cap':cap,
                'lower':summarize(lo[:,j],cap,alpha,False),'upper':summarize(up[:,j],cap,alpha,True)})
  best=max(rows,key=lambda x:x['lower']['mean']);best_ci=max(rows,key=lambda x:x['lower']['CI'])
  models[model]={'profile':rows,'best_lower_MC':best,'best_lower_CI':max(1.,best_ci['lower']['CI']),
                 'lower_CI_anchor_includes_L1':True}
 result={'registration_SHA256':sha(HERE/'HIGH_DIM_REGISTRATION.json'),'models':models,'N':N,'seed':SEED,
 'mass_weight_mean':float(np.mean(weights)),'mass_weight_SE_descriptive':float(np.std(weights,ddof=1)/math.sqrt(N)),
 'weight_max':max(weights),'uniform_rows':sum(b<0 for b,j in sample_tags),
 'tilted_b_counts':np.bincount([b for b,j in sample_tags if b>=0],minlength=4).tolist(),
 'count_sample_stream_SHA256':datahash.hexdigest(),'relative_variance_cap_quantiles':np.quantile(rb,[0,.5,.9,.99,1]).tolist(),
 'class_probability_diagnostics':extra,'sampling_and_bounds_seconds':time.time()-sample_start,
 'total_seconds':time.time()-started,'ru_maxrss_raw':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
 'ru_maxrss_unit':'bytes'if sys.platform=='darwin'else'KiB',
 'scope':'finite natural-grid true-source conditional tail envelopes as numerical evidence; ordinary FP not interval certified; not original matching asymptotic source and no C_n/log gain claim'}
 target.write_text(json.dumps(result,indent=2)+'\n')
 for m in MODELS:
  b=models[m]['best_lower_MC'];print(m,'best_lower',b['lower']['mean'],'h',b['h'],'ESS',b['lower']['ESS'],'maxLowerCI',models[m]['best_lower_CI'],flush=True)
 print('DONE mass',result['mass_weight_mean'],'total',result['total_seconds'],'RSS',result['ru_maxrss_raw'],flush=True)

if __name__=='__main__':main()
