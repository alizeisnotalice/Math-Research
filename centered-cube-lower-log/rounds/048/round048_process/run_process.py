"""Stationary OU Gaussian proxy, exact whole-path anchor mixture RN."""
from pathlib import Path
import hashlib,json,time,sys
import numpy as np
from scipy.special import log_ndtr
from scipy.optimize import brentq
HERE=Path(__file__).resolve().parent

def fixed_best(v):
 s=np.sqrt(v);t=brentq(lambda x:np.exp(-x*x/2-.5*np.log(2*np.pi)-log_ndtr(-x))-s,0,s+1)
 h=s*t-v/2;return float(h),float(np.exp(h+log_ndtr(-(h+v/2)/s)))

def profile(logx,lw,hs,part):
 N=len(logx);order=np.argsort(logx);x=logx[order];w=np.exp(lw[order]);p=part[order]
 tail=np.r_[np.cumsum(w[::-1])[::-1],0.];sq=np.r_[np.cumsum((w*w)[::-1])[::-1],0.]
 pt=np.r_[np.cumsum((w*p)[::-1])[::-1],0.];k=np.searchsorted(x,hs,side='right')
 m=np.exp(hs)*tail[k]/N;var=np.exp(2*hs)*np.maximum(0,sq[k]-tail[k]**2/N)/(N-1)
 return dict(mean=m.tolist(),variance=var.tolist(),hits=(N-k).tolist(),
             ess=np.divide(tail[k]**2,sq[k],out=np.zeros_like(m),where=sq[k]>0).tolist(),
             participation=np.divide(pt[k],tail[k],out=np.zeros_like(m),where=tail[k]>0).tolist())

def sample(v,c,R,T,dt,N,seed,hs):
 rng=np.random.default_rng(seed);J=round(T/dt)+1;rho=np.exp(-dt/v);innov=np.sqrt(-np.expm1(-2*dt/v))
 mode=(rng.random(N)>=.3);bank=rng.integers(R,size=N);anchor=rng.integers(J,size=N);idx=np.arange(N)
 z0=rng.standard_normal(N)+mode*np.sqrt(c*v);u=rng.standard_normal((N,R))
 common=-v/2+np.sqrt(c*v)*z0;s=np.sqrt((1-c)*v)
 lm=np.full(N,-np.inf);lsum=np.full(N,-np.inf);partmax=np.ones(N)
 for j in range(J):
  if j:u=rho*u+innov*rng.standard_normal((N,R))
  logits=s*u;logits[idx,bank]+=mode*((1-c)*v)*np.exp(-abs(j-anchor)*dt/v)
  maximum=logits.max(axis=1);ew=np.exp(logits-maximum[:,None]);sw=ew.sum(axis=1)
  ls=common+maximum+np.log(sw)-np.log(R)
  part=(ew*ew).sum(axis=1)/(sw*sw)
  new=ls>lm;partmax[new]=part[new];lm=np.maximum(lm,ls);lsum=np.logaddexp(lsum,ls)
  if j==0:l0=ls.copy();part0=part.copy()
 lB=lsum-np.log(J);lq=np.logaddexp(np.log(.3),np.log(.7)+lB);lw=-lq
 assert np.exp(lw+lm).max()<=J/.7*(1+1e-11)
 def ctl(a):return dict(mean=float(a.mean()),SE=float(a.std(ddof=1)/np.sqrt(N)))
 return dict(local=profile(lm,lw,hs,partmax),node=profile(l0,lw,hs,part0),
             RN_mass=ctl(np.exp(lw)),mean_B=ctl(np.exp(lw+lB)),sampled_max_wM=float(np.exp(lw+lm).max()))

def main():
 plan=json.loads((HERE/'FROZEN_PLAN.json').read_text());start=time.perf_counter();rows=[];N=plan['samples_per_seed'];family=96;L=np.log(4*family/.01)
 for vi,v in enumerate(plan['v']):
  for ci,c in enumerate(plan['c']):
   for ri,R in enumerate(plan['R']):
    for ti,T in enumerate(plan['T']):
     J=round(T/plan['dt'])+1;B=J/.7;h0=v/2-2
     hs=np.unique(np.r_[np.arange(-4,v/2+6.001,.5),[-1,-.5,-.1,0,.1],h0-np.log(R)+np.arange(-3,4.001,.5)])
     seeds=[s+10000*vi+1000*ci+100*ri+ti for s in plan['base_seeds']]
     runs=[sample(v,c,R,T,plan['dt'],N,s,hs) for s in seeds]
     selected={}
     for target in ['local','node']:
      a,b=runs[0][target],runs[1][target];k=int(np.argmax(a['mean']));m=b['mean'][k];var=b['variance'][k]
      radius=float(np.sqrt(2*var*L/N)+7*B*L/(3*(N-1)))
      selected[target]=dict(h=float(hs[k]),mean=m,variance=var,empirical_SE=float(np.sqrt(var/N)),
               simultaneous99_radius=radius,lower=max(0,m-radius),upper=m+radius,
               participation=b['participation'][k],weighted_ESS=b['ess'][k])
     bh,bw=fixed_best(v)
     rows.append(dict(v=v,c=c,R=R,T=T,J=J,seeds=seeds,h=hs.tolist(),runs=runs,selected=selected,
               fixed_single_node_best_h=bh,fixed_single_node_best_W=bw,
               local_gain_over_own_node=selected['local']['mean']/selected['node']['mean'],
               local_gain_over_fixed_single=selected['local']['mean']/bw))
   print(json.dumps({'configs':len(rows),'elapsed_seconds':time.perf_counter()-start}),flush=True)
 result=dict(scope='Gaussian OU proxy only, not actual HL; finite-grid extrema only.',rows=rows,
             elapsed_seconds=time.perf_counter()-start,plan_sha256=hashlib.sha256((HERE/'FROZEN_PLAN.json').read_bytes()).hexdigest(),
             script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             confidence='96 seedA-selected targets tested independently by seedB, two-sided empirical Bernstein union bound alpha=.01; capJ/.7; logfactor log(4*96/.01)',
             runtime=dict(python=sys.version,numpy=np.__version__))
 (HERE/'RESULTS.json').write_text(json.dumps(result,separators=(',',':'))+'\n')
 print(json.dumps({'configs':len(rows),'elapsed_seconds':result['elapsed_seconds']}))
if __name__=='__main__':main()
