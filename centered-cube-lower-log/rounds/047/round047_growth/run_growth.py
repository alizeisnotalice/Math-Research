"""Gaussian proxy only: defensive whole-vector anchor-mixture importance sampling."""
from pathlib import Path
import json,time,hashlib,sys
import numpy as np
from scipy.special import logsumexp,log_ndtr
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent

def fixed(h,v):
    return np.exp(np.asarray(h)+log_ndtr(-(np.asarray(h)+v/2)/np.sqrt(v)))

def fixed_best(v):
    sv=np.sqrt(v)
    t=brentq(lambda x:np.exp(-x*x/2-.5*np.log(2*np.pi)-log_ndtr(-x))-sv,0,sv+1)
    h=sv*t-v/2
    return float(h),float(fixed(h,v))

def simulate(v,R,c,N,seed,hs):
    rng=np.random.default_rng(seed);sv=np.sqrt(v)
    ts=np.array([0.,sv]);pis=np.array([.3,.7])
    mode=rng.choice(2,size=N,p=pis);anchor=rng.integers(R,size=N);shift=ts[mode]
    z0=rng.standard_normal(N)+shift*np.sqrt(c)
    zb=rng.standard_normal((N,R));zb[np.arange(N),anchor]+=shift*np.sqrt(1-c)
    y=np.sqrt(c)*z0[:,None]+np.sqrt(1-c)*zb
    logS=logsumexp(-v/2+sv*y,axis=1)-np.log(R)
    logq=np.logaddexp(np.log(.3),np.log(.7)+logS)
    lw=-logq;w=np.exp(lw)
    assert w.max()<=1/.3*(1+1e-12)
    # Exactly: q/P=.3+.7*S; the shifted part sums every anchor, not the sampled anchor alone.
    order=np.argsort(logS);sx=logS[order];sw=w[order]
    tw=np.r_[np.cumsum(sw[::-1])[::-1],0.];tw2=np.r_[np.cumsum((sw*sw)[::-1])[::-1],0.]
    k=np.searchsorted(sx,hs,side='right');tail=tw[k];tail2=tw2[k]
    estimate=np.exp(hs)*tail/N
    se=np.exp(hs)*np.sqrt(np.maximum(0,tail2/N-(tail/N)**2)/(N-1))
    ess=np.divide(tail*tail,tail2,out=np.zeros_like(tail),where=tail2>0)
    ws=np.exp(lw+logS)
    assert ws.max()<=1/.7*(1+1e-12)
    # e^h*w*1{S>e^h} <= w*S <= 10/7, for every threshold, without sampled-max assumptions.
    def meanse(a):return {'mean':float(a.mean()),'se':float(a.std(ddof=1)/np.sqrt(N))}
    return dict(W=estimate.tolist(),SE=se.tolist(),tail_ess=ess.tolist(),
                proposal_tail_hits=(N-k).tolist(),RN_mass=meanse(w),mean_one_control=meanse(ws))

def main():
    plan=json.loads((HERE/'FROZEN_PLAN.json').read_text());start=time.perf_counter();rows=[]
    # Full fixed set selected after the saved cost pilot. No configuration is selected by outcome.
    for vi,v in enumerate(plan['v']):
        bh,bw=fixed_best(v)
        for ci,c in enumerate(plan['c']):
            for ri,R in enumerate(plan['R']):
                h0=v/2-2
                hs=np.unique(np.r_[np.arange(-6,h0+6.001,.5),[-1,-.5,-.25,-.1,0,.1,.25],
                                          h0-np.log(R)+np.arange(-4,4.001,.5),h0])
                seeds=[s+1000*vi+100*ci+ri for s in plan['base_seeds']]
                tests=[simulate(v,R,c,plan['samples_per_seed'],s,hs) for s in seeds]
                a,b=tests;aa=np.array(a['W']);bb=np.array(b['W']);sa=np.array(a['SE']);sb=np.array(b['SE'])
                combined=(aa+bb)/2;se=np.sqrt(sa*sa+sb*sb)/2
                discovery=int(aa.argmax());combinedbest=int(combined.argmax());hi=int(np.argmin(abs(hs-h0)))
                exactfixed=fixed(hs,v)
                sj=np.exp(hs+np.log(R)+log_ndtr(-(hs+v/2+np.log(R))/np.sqrt(v)))
                approx=np.exp(-np.log(R)**2/(2*v));ratio=float(combined[hi]/exactfixed[hi])
                rows.append(dict(v=v,R=R,c=c,seeds=seeds,h=hs.tolist(),fixed_W=exactfixed.tolist(),
                    single_big_jump_W=sj.tolist(),seed_runs=tests,W=combined.tolist(),SE=se.tolist(),
                    ratio_to_fixed_same_h=(combined/exactfixed).tolist(),fixed_global_best_h=bh,fixed_global_best_W=bw,
                    discovery_h=float(hs[discovery]),independent_W_at_discovery_h=float(bb[discovery]),
                    independent_SE_at_discovery_h=float(sb[discovery]),independent_gain_over_fixed_global=float(bb[discovery]/bw),
                    combined_grid_best_h=float(hs[combinedbest]),combined_grid_best_W=float(combined[combinedbest]),
                    combined_grid_gain_over_fixed_global=float(combined[combinedbest]/bw),
                    high_baseline=dict(h=h0,W=float(combined[hi]),SE=float(se[hi]),
                        ratio_to_fixed=ratio,ratio_SE=float(se[hi]/exactfixed[hi]),
                        single_big_jump_exact_ratio=float(sj[hi]/exactfixed[hi]),
                        exp_logR_squared_approx=float(approx),observed_to_approx=float(ratio/approx),
                        observed_to_single_big_jump_exact=float(combined[hi]/sj[hi]))))
        print(json.dumps({'completed_configs':len(rows),'elapsed_seconds':time.perf_counter()-start}),flush=True)
    out=dict(scope='Gaussian equicorrelated mean-one lognormal proxy; no actual HL coordinates or C_n claim.',
             proposal='q=.3 p+.7 R^-1 sum_anchor N(mu_anchor(sqrt(v)),I); exact full q/P=.3+.7*S.',
             deterministic_target_cap=1/.7,
             event='S=R^-1 sum_b exp(-v/2+sqrt(v)*(sqrt(c)Z0+sqrt(1-c)Zb)); strict S>exp(h)',
             plan_sha256=hashlib.sha256((HERE/'FROZEN_PLAN.json').read_bytes()).hexdigest(),
             script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             elapsed_seconds=time.perf_counter()-start,runtime=dict(python=sys.version,numpy=np.__version__),rows=rows)
    (HERE/'RESULTS.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
    summary=[{k:r[k] for k in ['v','R','c','combined_grid_best_h','combined_grid_best_W','combined_grid_gain_over_fixed_global',
               'discovery_h','independent_gain_over_fixed_global','high_baseline']} for r in rows]
    (HERE/'SUMMARY.json').write_text(json.dumps({'elapsed_seconds':out['elapsed_seconds'],'rows':summary},indent=2)+'\n')
    print(json.dumps({'completed_configs':len(rows),'elapsed_seconds':out['elapsed_seconds'],'output':str(HERE/'RESULTS.json')}))

if __name__=='__main__':main()
