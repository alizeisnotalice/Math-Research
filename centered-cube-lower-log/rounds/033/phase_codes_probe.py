"""Frozen local-vs-global phase pilot; nine lengths, not asymptotic evidence."""
from pathlib import Path
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','4')
os.environ.setdefault('OMP_NUM_THREADS','4')
import sys,importlib.util,json,math,time
import numpy as np
OUT=Path(__file__).resolve().parent
SOURCE=Path(sys.argv[1]) if len(sys.argv)>1 else OUT.parent/'032/multinomial_probe.py'
spec=importlib.util.spec_from_file_location('source32',SOURCE)
s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
N=512;SEED=202610290331;R=8
THRESHOLDS=np.exp(np.array([0,.5,1,1.5,2,3,4]))

def main():
    start=time.time();rng=np.random.default_rng(SEED);src=s.Source(8,0x11d,4,R)
    J=len(s.LENGTHS);Q=src.q
    codes={'global':np.arange(R)*(Q//R),
           'local_quarter':np.arange(R)*6,
           'local_one':np.arange(R)*24,
           'local_four':np.arange(R)*96}
    pool=np.concatenate(list(codes.values()));P=len(pool)
    base=[src.base_query(j) for j in range(J)]
    # All 4 codebooks are included, including their duplicated zero phases.
    counts=np.empty((N,Q));labels=[]
    for row in range(N):
        if rng.random()<.5:
            probs=np.full(Q,1/Q);labels.append([-1,-1])
        else:
            j=int(rng.integers(J));r=int(rng.integers(P))
            m=np.roll(base[j][0],int(pool[r]));probs=m/m.sum();labels.append([j,r])
        counts[row]=rng.multinomial(src.n,probs)
    logs=np.empty((N,J,P));rel=np.empty_like(logs);norms=np.empty((J,P))
    for j,(m,tr) in enumerate(base):
        means=np.column_stack([np.roll(m,int(t)) for t in pool])
        trans=np.column_stack([np.roll(tr,int(t)) for t in pool])
        norms[j]=np.array([np.roll(m,int(t)).mean() for t in pool])
        logs[:,j]=counts@np.log(means)
        number=counts@trans.astype(float)
        rel[:,j]=np.expm1(src.n*np.log1p(src.t0)+number*(np.log1p(src.t1)-np.log1p(src.t0)))
    logq=np.logaddexp(math.log(.5),math.log(.5)+s.lse(
                   (logs-src.n*np.log(norms)[None,:,:]).reshape(N,-1))-math.log(J*P))
    weights=np.exp(-logq);assert weights.max()<=2*(1+1e-12)
    inflation=math.exp(max(0.,float(src.n*np.log(norms).max())))
    # Product uses the zero phase. Union of conditional Cantelli upper bounds.
    lp=logs[:,:,0];rp=rel[:,:,0]
    product_upper={}
    for t in THRESHOLDS:
        a=np.clip(math.log(t)-lp,0,350);gap=np.expm1(a)
        upper=np.where(lp>=math.log(t),1.,np.divide(rp,rp+gap*gap,
                       out=np.ones_like(rp),where=(rp+gap*gap)>0))
        product_upper[float(t)]=np.minimum(1.,upper.sum(axis=1))
    result={}
    for ci,(name,phases) in enumerate(codes.items()):
        z=logs[:,:,ci*R:(ci+1)*R];rv=rel[:,:,ci*R:(ci+1)*R]
        lm=s.lse(z,axis=2)-math.log(R)
        phase_weights=np.exp(z-s.lse(z,axis=2)[:,:,None])
        cv=(phase_weights*np.sqrt(rv)).sum(axis=2)**2
        best=np.argmax(z,axis=1)
        selected=np.take_along_axis(z,best[:,None,:],axis=1)[:,0,:]
        logu=s.lse(selected)-math.log(R);logv=lm.max(axis=1)
        assert np.min(logu-logv)>-1e-10
        table=[]
        for t in THRESHOLDS:
            gap=np.maximum(0.,-np.expm1(np.minimum(700.,math.log(t)-lm)))
            lower=np.divide(gap*gap,cv+gap*gap,out=np.zeros_like(gap),where=(gap>0)).max(axis=1)
            score=t*weights*lower
            cap=min(2*t,2*J*P/R*inflation)
            assert score.max()<=cap*(1+1e-10)
            pu=t*weights*product_upper[float(t)];diff=score-pu
            table.append(dict(threshold=float(t),shared_lower_MC=float(score.mean()),
                shared_lower_SE=float(score.std(ddof=1)/np.sqrt(N)),
                shared_contribution_ESS=float(score.sum()**2/np.dot(score,score)) if np.any(score) else 0.,
                product_upper_MC=float(pu.mean()),
                paired_lower_minus_product_upper_MC=float(diff.mean()),
                paired_difference_SE=float(diff.std(ddof=1)/np.sqrt(N)),
                maximum_contribution=float(score.max()),contribution_cap=float(cap)))
        # Post-pilot structural diagnostics on the SAME frozen sample.
        # Rank distributions are descriptive self-normalized weighted empirical laws.
        rank=np.empty_like(selected)
        for r in range(R):
            order=np.argsort(selected[:,r])
            rank[order,r]=(np.cumsum(weights[order])-.5*weights[order])/weights.sum()
        rank_mean=(weights[:,None]*rank).sum(axis=0)/weights.sum()
        centered=rank-rank_mean
        cov=(centered*weights[:,None]).T@centered/weights.sum()
        corr=cov/np.sqrt(np.diag(cov)[:,None]*np.diag(cov)[None,:])
        common=z[np.arange(N),np.argmax(lm,axis=1)]
        component_weights=np.exp(common-s.lse(common)[:,None])
        effective=1/np.sum(component_weights*component_weights,axis=1)
        result[name]=dict(phase_bin_shifts=phases.tolist(),
            phase_displacements=(phases/Q).tolist(),
            proposal_log_Ubar_over_Vbar_quantiles=np.quantile(logu-logv,[0,.5,.9,1]).tolist(),
            target_weighted_mean_conditional_log_gap=float(np.mean(weights*(logu-logv))),
            descriptive_weighted_phase_max_rank_correlation=float(corr[~np.eye(R,dtype=bool)].mean()),
            descriptive_weighted_effective_phases_at_mean_winner=float(np.dot(weights,effective)/weights.sum()),
            fixed_thresholds=table)
    out=dict(seed=SEED,N=N,numpy_version=np.__version__,source=src.metadata(),
       lengths=s.LENGTHS.tolist(),proposal=dict(component_pool=P,target_R=R,
       uniform_rows=sum(j<0 for j,r in labels),weight_mean=float(weights.mean()),
       weight_mean_SE=float(weights.std(ddof=1)/np.sqrt(N)),weight_max=float(weights.max()),
       known_normalizer_inflation=inflation),results=result,
       pointwise_95pct_signed_Hoeffding_coefficient_factor=float(4*np.sqrt(np.log(20)/ (2*N))),
       runtime_seconds=time.time()-start,
       scope='Frozen nine-length code pilot. MC bounds lack interval certification; no log lower bound.')
    (OUT/'phase_codes_results.json').write_text(json.dumps(out,indent=2))
    print(json.dumps({k:dict(log_gap=d['target_weighted_mean_conditional_log_gap'],
                   best_fixed_lower=max(t['shared_lower_MC'] for t in d['fixed_thresholds']),
                   largest_paired_difference=max(t['paired_lower_minus_product_upper_MC'] for t in d['fixed_thresholds']))
                   for k,d in result.items()},indent=2))

if __name__=='__main__':main()
