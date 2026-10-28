"""Bin-conditional moment interface: analytic piece integrals, exact count sampling.
The floating moment calculations are diagnostics, not interval-certified inputs.
"""
from pathlib import Path
import importlib.util
import numpy as np
import json, math, sys

OUT=Path(__file__).resolve().parent
SOURCE=(Path(sys.argv[1]) if len(sys.argv)>1 else
        OUT.parent/'031/probe.py')
spec=importlib.util.spec_from_file_location('cube_source31',SOURCE)
s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
R=len(s.PHASES);omega=2*np.pi

def atom_density(t):
    u=np.mod(t,1)
    return sum((np.abs(u-a)<s.DELTA/2).astype(float)
               for a in s.ATOMS)/(s.K*s.DELTA)

def moments(Q,L):
    knots=[np.arange(Q+1)/Q]
    for phase in s.PHASES:
        for sign in [-1,1]:
            for edge in [-1,1]:
                knots.append(np.mod(s.ATOMS+phase+sign*L/2
                                    +edge*s.DELTA/2,1))
    knots=np.unique(np.concatenate(knots));means=np.zeros((Q,R))
    second=np.zeros((Q,R,R))
    d=s.P/s.K*s.C1*np.sin(np.pi*L)/L
    for left,right in zip(knots[:-1],knots[1:]):
        h=right-left
        if h<1e-15:continue
        mid=(left+right)/2;bin_id=min(int(mid*Q),Q-1)
        x=mid-s.PHASES;psi=omega*x;z=omega*h/2
        a=s.average_from_shifted(x,L)-d*np.cos(psi)
        b=s.P/L*(atom_density(x+L/2)-atom_density(x-L/2))
        sinc=np.sinc(z/np.pi)
        usin=(sinc-np.cos(z))/omega
        cmean=np.cos(psi)*sinc
        ucos=-np.sin(psi)*usin
        pm=a[:,None]*a[None,:]+b[:,None]*b[None,:]*h*h/12
        pc=d*(a[:,None]*cmean[None,:]+b[:,None]*ucos[None,:])
        cc=.5*d*d*(np.cos(psi[:,None]-psi[None,:])
                         +np.cos(psi[:,None]+psi[None,:])*np.sinc(2*z/np.pi))
        means[bin_id]+=(a+d*cmean)*h*Q
        second[bin_id]+=(pm+pc+pc.T+cc)*h*Q
    assert np.all(means>0) and np.all(second>0)
    return means,second

def conditional(counts,means,second):
    lm=counts@np.log(means)
    lm2=counts@np.log(second.reshape(len(means),-1))
    m=np.exp(lm).mean(axis=1)
    e2=np.exp(lm2).mean(axis=1)
    vv=e2-m*m
    assert np.min(vv/(m*m))>-1e-10
    return m,np.maximum(vv,0)

def actual(x,lengths):
    shifted=x[:,:,None]-s.PHASES
    vals=[]
    for L in lengths:
        vals.append(np.exp(np.log(s.average_from_shifted(shifted,L)).sum(axis=1)).mean(axis=1))
    return np.column_stack(vals)

def main():
    rng=np.random.default_rng(20261029032)
    lengths=np.linspace(10/9,11/9,9);N=4096;rows=[]
    for n in [8,32,64]:
        x=rng.random((N,n));g=actual(x,lengths);mx=g.max(axis=1)
        for Q in [16,64,256]:
            ids=np.minimum((x*Q).astype(int),Q-1)
            counts=np.zeros((N,Q),dtype=np.int64)
            np.add.at(counts,(np.repeat(np.arange(N),n),ids.ravel()),1)
            mm=[];vv=[]
            for L in lengths:
                a,b=moments(Q,L);m,v=conditional(counts,a,b)
                mm.append(m);vv.append(v)
            mm=np.column_stack(mm);vv=np.column_stack(vv)
            # Fixed thresholds. Select L only from counts, via its conditional bound.
            for t in [0.8,1.,1.2,1.5,2.]:
                gap=np.maximum(mm-t,0)
                lower=np.divide(gap*gap,vv+gap*gap,out=np.zeros_like(gap),
                                where=(gap>0)).max(axis=1)
                indicator=(mx>t).astype(float);diff=indicator-lower
                rows.append(dict(n=n,Q=Q,threshold=t,N=N,
                    conditional_probability_lower_estimate=float(lower.mean()),
                    actual_grid_probability_estimate=float(indicator.mean()),
                    paired_difference=float(diff.mean()),
                    paired_difference_se=float(diff.std(ddof=1)/np.sqrt(N)),
                    median_conditional_relative_variance=float(np.median(vv/(mm*mm)))))
    # Independently validate one exact bin pattern via fresh within-bin samples.
    check=[];n=8;Q=16;ids=np.array([0,1,3,5,7,10,12,15]);counts=np.bincount(ids,minlength=Q)[None,:]
    M=65536;x=(ids[None,:]+rng.random((M,n)))/Q
    for L in [10/9,7/6,11/9]:
        a,b=moments(Q,L);m,v=conditional(counts,a,b)
        g=actual(x,[L])[:,0]
        check.append(dict(L=float(L),analytic_conditional_mean=float(m[0]),
                    analytic_conditional_variance=float(v[0]),MC_mean=float(g.mean()),
                    MC_variance=float(g.var()),MC_mean_se=float(g.std(ddof=1)/np.sqrt(M))))
    # Sample bin counts directly from a mixture of all (L,phase) product tilts.
    n=64;Q=64;Nis=8192;J=len(lengths)
    tabs=[moments(Q,L) for L in lengths]
    probs=np.stack([a/Q for a,b in tabs])
    normalization_error=float(np.max(abs(probs.sum(axis=1)-1)))
    probs/=probs.sum(axis=1,keepdims=True)
    labels=rng.integers(J*R,size=Nis);counts=np.zeros((Nis,Q),dtype=np.int64)
    for label in range(J*R):
        take=np.flatnonzero(labels==label)
        counts[take]=rng.multinomial(n,probs[label//R,:,label%R],size=len(take))
    mm=[];vv=[];log_parts=[]
    for j,(a,b) in enumerate(tabs):
        m,v=conditional(counts,a,b);mm.append(m);vv.append(v)
        log_parts.append(counts@np.log(probs[j]*Q))
    mm=np.column_stack(mm);vv=np.column_stack(vv)
    lp=np.concatenate(log_parts,axis=1)
    log_q=s.lse(lp)-np.log(J*R)
    importance=[]
    for t in [0.8,1.,1.2,1.5,2.]:
        gap=np.maximum(mm-t,0)
        lower=np.divide(gap*gap,vv+gap*gap,out=np.zeros_like(gap),where=(gap>0)).max(axis=1)
        contribution=t*lower*np.exp(-log_q)
        assert contribution.max()<=J*(1+1e-10)
        importance.append(dict(threshold=t,N=Nis,weak_coefficient_lower_estimate=float(contribution.mean()),
                         coefficient_se=float(contribution.std(ddof=1)/np.sqrt(Nis)),
                         maximum_contribution=float(contribution.max()),bound=J,
                         pointwise_95pct_Hoeffding_radius=J*np.sqrt(np.log(40)/(2*Nis))))
    out=dict(seed=20261029032,numpy_version=np.__version__,lengths=lengths.tolist(),
             rows=rows,within_bin_checks=check,
             count_importance=importance,proposal_normalization_float_error=normalization_error,
             status='Numerical conditional-moment diagnostics; no interval or asymptotic certificate.')
    (OUT/'conditional_checks.json').write_text(json.dumps(out,indent=2))
    print(json.dumps({'within_bin_checks':check,'comparison_rows':len(rows),
         'minimum_paired_difference':min(r['paired_difference'] for r in rows)},indent=2))

if __name__=='__main__':main()
