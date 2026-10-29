"""True continuous-torus count-tail gates and matched-amplitude mixtures."""
from pathlib import Path
import hashlib,json,time,math
import numpy as np
from scipy.special import gammaln,logsumexp
HERE=Path(__file__).resolve().parent

def pb(p):
 N,n=p.shape;dp=np.zeros((N,n+1));dp[:,0]=1
 for j in range(n):
  old=dp[:,:j+1].copy();dp[:,1:j+2]=dp[:,1:j+2]*(1-p[:,j,None])+old*p[:,j,None];dp[:,0]*=1-p[:,j]
 return dp

def profile(logz,hs):
 return np.array([np.exp(h)*np.mean(logz>h) for h in hs])

def selected(logz,hs,Ntrain):
 w=profile(logz[:Ntrain],hs);j=int(w.argmax());h=hs[j];sample=np.exp(h)*(logz[Ntrain:]>h)
 return dict(h=float(h),train_W=float(w[j]),test_W=float(sample.mean()),test_SE=float(sample.std(ddof=1)/np.sqrt(len(sample))),test_hits=int((logz[Ntrain:]>h).sum()))

def main():
 plan=json.loads((HERE/'FROZEN_PLAN.json').read_text());start=time.perf_counter();rows=[];checks=[]
 Ntrain=plan['N_train'];N=Ntrain+plan['N_test'];hs=np.arange(-4,5.001,.125);pi=.5;c=.05;tau=.2
 for n in plan['n']:
  rng=np.random.default_rng(plan['base_seed']+n);x=rng.random((N,n));k=np.arange(n+1)
  logprior=gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)-n*np.log(2);prior=np.exp(logprior)
  specs=[('tail_zero',k==0),('lower_eighth',k<=n//8),('lower_quarter',k<=n//4),('lower_half',k<=n//2),
         ('upper_eighth',k>=n-n//8),('two_tails',(k<=n//8)|(k>=n-n//8)),('endpoint_pair',(k==0)|(k==n))]
  gates=[]
  for name,mask in specs:
   lw=np.where(mask,0.,-np.inf);lm=logsumexp(logprior+lw);nu=np.exp(logprior+lw-lm)
   gates.append((name,lw,lm,nu))
  lw=-logprior-np.log(n+1);gates.append(('uniform_count_prior',lw,logsumexp(logprior+lw),np.ones(n+1)/(n+1)))
  alpha=1+tau*(k/n-pi);mean_alpha=np.array([nu@alpha for name,lw,lm,nu in gates])
  bank=np.empty((n+1,N,3));gatez=np.empty((len(gates),N,3));mixz=np.empty_like(gatez);meanz=np.empty_like(gatez)
  for li,Ltext in enumerate(plan['scales']):
   num,den=map(int,Ltext.split('/'));L=num/den;d=.9*np.sinc(L)*np.cos(2*np.pi*x);Y=1+d;p=pi+c*d/Y
   dp=pb(p);base=np.log1p(d).sum(1)
   for ki,a in enumerate(alpha):bank[ki,:,li]=np.log1p(a*d).sum(1)
   for gi,(name,lw,lm,nu) in enumerate(gates):
    gatez[gi,:,li]=base+logsumexp(np.log(dp)+lw,axis=1)-lm
    mixz[gi,:,li]=logsumexp(bank[:,:,li]+np.log(nu,where=nu>0,out=np.full_like(nu,-np.inf))[:,None],axis=0)
    meanz[gi,:,li]=np.log1p(mean_alpha[gi]*d).sum(1)
   checks.append(dict(n=n,L=Ltext,max_PB_mass_error=float(abs(dp.sum(1)-1).max()),p_min=float(p.min()),p_max=float(p.max()),
                      extreme_posterior_min=float(dp[:,0].min())))
  maxbank=bank.max(axis=2);train_single=np.array([profile(z[:Ntrain],hs) for z in maxbank]);bi,bh=np.unravel_index(train_single.argmax(),train_single.shape)
  bestsingle=selected(maxbank[bi],hs,Ntrain);bestsingle.update(k=int(bi),alpha=float(alpha[bi]))
  for gi,(name,lw,lm,nu) in enumerate(gates):
   mg=gatez[gi].max(axis=1);mm=mixz[gi].max(axis=1);ma=meanz[gi].max(axis=1)
   selg=selected(mg,hs,Ntrain);selm=selected(mm,hs,Ntrain);sela=selected(ma,hs,Ntrain)
   mu=float(nu@k);var=float(nu@((k-mu)**2))
   rows.append(dict(n=n,gate=name,source_log_weight_normalizer=float(lm),prior_mean_count=mu,prior_count_variance=var,
      mean_amplitude=float(mean_alpha[gi]),amplitude_std=float(np.sqrt(nu@((alpha-mean_alpha[gi])**2))),
      gate_selected=selg,matched_mixture_selected=selm,mean_product_selected=sela,best_single_grid=bestsingle,
      gate_over_best_single_test=float(selg['test_W']/bestsingle['test_W']),
      gate_over_matched_mixture_test=float(selg['test_W']/selm['test_W']),
      max_query_log_gate_over_matched_mixture=float((gatez[gi]-mixz[gi]).max()),
      full_profile=dict(h=hs.tolist(),gate_train=profile(mg[:Ntrain],hs).tolist(),gate_test=profile(mg[Ntrain:],hs).tolist(),
                        mixture_test=profile(mm[Ntrain:],hs).tolist(),mean_product_test=profile(ma[Ntrain:],hs).tolist()),
      MC_mass_gate_per_scale=np.exp(gatez[gi]).mean(axis=0).tolist(),MC_mass_matched_mixture_per_scale=np.exp(mixz[gi]).mean(axis=0).tolist()))
   rows[-1]['MC_L1_mixture_minus_mean_product_per_scale']=np.abs(np.exp(mixz[gi])-np.exp(meanz[gi])).mean(axis=0).tolist()
  # k0 andkn are true endpoint products, with no approximation.
  assert abs(gatez[0]-bank[0]).max()<1e-11
  print(json.dumps({'completed_n':n,'seconds':time.perf_counter()-start}),flush=True)
 out=dict(scope=plan['scope'],plan_sha256=hashlib.sha256((HERE/'FROZEN_PLAN.json').read_bytes()).hexdigest(),
          script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),rows=rows,checks=checks,
          source_legality=dict(rho_min=.1,rho1_min=.005,rho0_min=.095,pi=pi,c=c,tau=tau),elapsed_seconds=time.perf_counter()-start)
 (HERE/'RESULTS.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'rows':len(rows),'seconds':out['elapsed_seconds'],'max_test_gate_over_best_single':max(r['gate_over_best_single_test'] for r in rows)}))
if __name__=='__main__':main()
