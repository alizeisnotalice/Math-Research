#!/usr/bin/env python3
"""Exact finite-law checks and read-only numerical replay of round043_tilt."""
from fractions import Fraction as F
from itertools import product
from math import factorial,prod
from pathlib import Path
import hashlib,importlib.util,json,time
import numpy as np
HERE=Path(__file__).resolve().parent
SRC=HERE.parent/'round043_tilt'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def comps(n,K):
    if K==1:yield(n,)
    else:
        for a in range(n+1):
            for z in comps(n-a,K-1):yield(a,)+z
def dist(nums):
    s=sum(nums);return[F(v,s) for v in nums]
def posterior(eta,u,q):
    return [[eta[c]*(q+row[c])/(q+sum(e*v for e,v in zip(eta,row)))
             for c in range(len(eta))]for row in u]
def moments(p,nu):
    K=len(p[0]);mu=[sum(w*row[c]for w,row in zip(nu,p))for c in range(K)]
    M=[[sum(w*row[c]*row[d]for w,row in zip(nu,p))for d in range(K)]for c in range(K)]
    tau=sum(M[c][c]/mu[c]-mu[c]for c in range(K))
    return mu,M,tau
def pi(alpha,mu):return F(factorial(sum(alpha)),prod(factorial(a)for a in alpha))*prod(v**a for v,a in zip(mu,alpha))
def poly(rows):
    K=len(rows[0]);d={(0,)*K:F(1)}
    for row in rows:
        nd={}
        for a,v in d.items():
            for c,w in enumerate(row):
                key=a[:c]+(a[c]+1,)+a[c+1:];nd[key]=nd.get(key,0)+v*w
        d=nd
    return d
def matmul(a,b):return[[sum(x*y for x,y in zip(row,col))for col in zip(*b)]for row in a]
def complete_trace(B,n):
    K=len(B);A=[[F(int(i==j))for j in range(K)]for i in range(K)];tr=[]
    for j in range(n):A=matmul(A,B);tr.append(sum(A[i][i]for i in range(K)))
    h=[F(1)]
    for m in range(1,n+1):h.append(sum(tr[j-1]*h[m-j]for j in range(1,m+1))/m)
    return h[n]
def hist_trace(p,nu,mu,n,extra=None,oldmu=None):
    total=F(0)
    for inds in product(range(len(p)),repeat=n):
        w=prod(nu[i]for i in inds)
        if extra is not None:w*=prod(extra[i]for i in inds)
        d=poly([p[i]for i in inds])
        total+=w*sum(v*v/pi(a,mu if oldmu is None else oldmu)for a,v in d.items())
    return total

def physical_full_window_example():
    # Direct length5/4 centered integrals, WITHOUT q+ell kernel reduction.
    classes=[[(F(0),F(1,2),F(1,10)),(F(9,40),F(11,40),F(9))],
             [(F(1,2),F(1),F(1))]]
    L=F(5,4)
    def z(segments,x):
        total=F(0)
        for a,b,h in segments:
            for shift in range(-2,3):
                total+=h*max(F(0),min(x+L/2,b+shift)-max(x-L/2,a+shift))
        return total/L
    knots={F(0),F(1)}
    for segments in classes:
        for a,b,h in segments:
            for e in (a,b):
                for sign in (-1,1):knots.add((e+sign*L/2)%1)
    knots=sorted(knots);mass=F(0);m=F(0);nums=[F(0),F(0)];zz=[[F(0)]*2 for _ in range(2)]
    for left,right in zip(knots,knots[1:]):
        zl=[z(s,left)for s in classes];zr=[z(s,right)for s in classes]
        al,ar=sum(zl),sum(zr);h=right-left
        mass+=h*(al+ar)/2
        m+=h*(al*al+al*ar+ar*ar)/3
        for c in range(2):nums[c]+=h*(2*al*zl[c]+al*zr[c]+ar*zl[c]+2*ar*zr[c])/6
        for c in range(2):
            for d in range(2):zz[c][d]+=h*(2*zl[c]*zl[d]+zl[c]*zr[d]+zr[c]*zl[d]+2*zr[c]*zr[d])/6
    assert mass==1 and m==F(505697,500000)
    mu=[num/m for num in nums]
    assert nums==[F(128591,250000),F(49703,100000)]
    assert mu==[F(257182,505697),F(248515,505697)]
    factor=sum(v*v/F(1,2)for v in mu)
    assert factor==F(255804572698,255729455809)>1
    source=json.loads((SRC/'EXACT_TILT_EXAMPLE.json').read_text())
    assert source['script_SHA256']==sha(SRC/'exact_tilt_example.py')
    assert source['m_t_exact']==str(m) and source['tilted_numerator_exact']==list(map(str,nums))
    gram=[[v/F(1,2)/m for v in row]for row in zz]
    price=json.loads((SRC/'PRICE.json').read_text())
    assert price['script_SHA256']==sha(SRC/'price_probe.py')
    assert price['exact_source_sha']==sha(SRC/'EXACT_TILT_EXAMPLE.json')
    assert gram==[[F(v)for v in row]for row in price['original_prior_Gram_exact']]
    chi=factor-1;trace=gram[0][0]+gram[1][1];det=gram[0][0]*gram[1][1]-gram[0][1]**2
    lower=1+chi;upper=1+2*chi
    assert lower**2-trace*lower+det<=0
    assert upper>trace/2 and upper**2-trace*upper+det>=0
    return dict(method='Direct physical centered length5/4 window, own rational overlaps and source-edge breakpoints',
        source_mass=str(mass),m1=str(m),numerators=list(map(str,nums)),mu=list(map(str,mu)),
        exact_Jensen_factor=str(factor),intervals=len(knots)-1,
        original_prior_Gram_exact=[[str(v)for v in row]for row in gram],
        gamma_lower=str(lower),gamma_upper=str(upper),exact_characteristic_polynomial_bounds_pass=True)

def main():
    started=time.time();rng=np.random.default_rng(430043)
    scalar=[];finite=0
    for q in (1,2,3,7):
        a,b=F(q,q+1),F(q+1,q);c=F(1,q*(q+1))
        whigh=(1-a)/(b-a);vmean=(1-whigh)*a+whigh*b
        variance=(1-whigh)*(a-vmean)**2+whigh*(b-vmean)**2
        assert vmean==1 and variance/vmean==c
        scalar.append(dict(q=q,a=str(a),b=str(b),extremal_high_probability=str(whigh),sharp_scalar_ratio=str(c)))
    saved=[];spectral_checks=0
    for K in (2,3,4,5):
        for q in (1,2,3):
            for trial in range(12):
                eta=dist([int(x)for x in rng.integers(1,10,size=K)])
                nu=dist([int(x)for x in rng.integers(1,10,size=3)])
                u=[[F(int(x),7)for x in rng.integers(0,8,size=K)]for j in range(3)]
                p=posterior(eta,u,q);mu,M,tau=moments(p,nu)
                assert all(sum(row)==1 for row in p)and sum(mu)==1
                assert 0<=tau<=F(1,q*(q+1));finite+=1
                kappa=F(1,q*(q+1));chi=sum(v*v/e for v,e in zip(mu,eta))-1
                oldtrace=sum(M[c][c]/eta[c]for c in range(K))
                assert oldtrace<=1+kappa and chi>=0
                h=np.array([[float(M[c][d]/(eta[c]*eta[d])) for d in range(K)]for c in range(K)])
                h*=np.sqrt(np.asarray(eta,dtype=float))[:,None]*np.sqrt(np.asarray(eta,dtype=float))[None,:]
                ev=np.linalg.eigvalsh(h);gamma=ev[-1]
                assert 1+float(chi)-1e-12<=gamma<=1+float(chi/(1-kappa))+1e-12
                assert ev[:-1].sum()<=float(kappa)+1e-12
                spectral_checks+=1
                if trial==0 and K<=3:saved.append((eta,nu,p,q))
    histchecks=0;mgfchecks=0
    for eta,nu,p,q in saved:
        mu,M,tau=moments(p,nu);K=len(eta)
        B=[[M[c][d]/mu[d]for d in range(K)]for c in range(K)]
        for n in (1,2,3,4):
            observed=hist_trace(p,nu,mu,n)
            assert observed==complete_trace(B,n)
            assert observed<=1/(1-F(1,q*(q+1)))
            histchecks+=1
        score=[F(2,3),F(1),F(5,4)]
        for t in (-3,-1,0,1,3):
            weights=[v*a**t for v,a in zip(nu,score)];m=sum(weights)
            tilted=[w/m for w in weights]
            mut,Mt,taut=moments(p,tilted)
            for n in (1,2,3):
                lhs=hist_trace(p,nu,mut,n,extra=[a**t for a in score])
                rhs=m**n*hist_trace(p,tilted,mut,n)
                assert lhs==rhs
                old=hist_trace(p,tilted,mut,n,oldmu=eta)
                ratio=max(pi(a,mut)/pi(a,eta)for a in comps(n,K))
                assert old<=ratio*hist_trace(p,tilted,mut,n)
                mismatch=sum(v*v/e for v,e in zip(mut,eta))
                assert old>=mismatch**n
                assert sum(pi(a,mut)**2/pi(a,eta)for a in comps(n,K))==mismatch**n
                mgfchecks+=1
    # Explicit demonstration that the new histogram prior differs from eta.
    eta=[F(1,3),F(2,3)];nu=[F(1,4),F(3,4)]
    p=posterior(eta,[[F(0),F(1)],[F(1),F(0)]],1)
    mu,M,tau=moments(p,nu)
    assert mu==[F(17,40),F(23,40)]
    toy=[]
    for n in (1,2,4,8,16,32,64):
        old=sum(pi(a,[F(2,3),F(1,3)])**2/pi(a,[F(1,2)]*2)for a in comps(n,2))
        assert old==F(10,9)**n
        toy.append(dict(n=n,wrong_trace=str(old),correct_trace='1'))

    # Recreate all source/query-law inputs read only, use centered variances.
    spec=importlib.util.spec_from_file_location('tilt43',SRC/'probe.py')
    mainprobe=importlib.util.module_from_spec(spec);spec.loader.exec_module(mainprobe)
    plan=json.loads((SRC/'PLAN.json').read_text());summary=json.loads((SRC/'SUMMARY.json').read_text())
    assert plan['script_SHA256']==sha(SRC/'probe.py')
    assert plan['channel42_SHA256']==sha(mainprobe.BASE)
    assert summary['plan_SHA256']==sha(SRC/'PLAN.json')
    c=mainprobe.CONFIG;rng=np.random.default_rng(c['seed']);groups={};count=0;maxerr=0.
    for K,q,ell,kind,shift in product(c['K'],c['q'],c['ell'],c['kinds'],c['shifts']):
        x,w,a,pp=mainprobe.posterior(K,q,ell,kind,shift)
        cases=[]
        for t in c['score_tilts']:
            logs=np.log(w)+(1+t)*np.log(a);cases.append(('score_tilt',np.exp(logs-logs.max())))
        for j in range(c['random_measures_per_channel']):
            logs=8*rng.standard_normal(len(w));cases.append(('random_query',np.exp(logs-logs.max())))
        for j in range(c['point_pairs_per_channel']):
            idx=rng.choice(len(w),2,replace=False);ww=np.zeros(len(w));z=float(rng.uniform(.01,.99));ww[idx]=[z,1-z];cases.append(('two_points',ww))
        for name,ww in cases:
            ww=ww/ww.sum();newprior=ww@pp
            value=float(np.sum((ww[:,None]*(pp-newprior)**2).sum(axis=0)/newprior))
            original,_=mainprobe.trace(pp,ww);maxerr=max(maxerr,abs(value-original))
            assert value<=1/(q*(q+1))+1e-12
            group=groups.setdefault(name,dict(count=0,max_trace=0.))
            group['count']+=1;group['max_trace']=max(group['max_trace'],value);count+=1
    assert count==summary['measures_evaluated']==9000
    for name,row in groups.items():
        assert row['count']==summary['groups'][name]['count']
        assert abs(row['max_trace']-summary['groups'][name]['max_trace'])<2e-12
    norm=json.loads((SRC/'NORMALIZATION.json').read_text())
    assert norm['script_SHA256']==sha(SRC/'normalization_probe.py')
    for row in norm['actual_score_tilt_examples']:
        spec={key:row[key]for key in ('K','q','ell','kind','shift')}
        x,w,a,pp=mainprobe.posterior(**spec)
        logs=np.log(w)+(1+row['t'])*np.log(a);ww=np.exp(logs-logs.max());ww/=ww.sum()
        mm=ww@pp
        assert max(abs(mm-np.array(row['eta_t'])))<1e-12
        assert abs(float(np.sum(mm*mm)*row['K'])-row['old_prior_Jensen_factor'])<1e-12
    out=dict(status='PASS',scalar_extremal_checks=scalar,finite_exact_posterior_law_checks=finite,
        exact_histogram_symmetric_trace_checks=histchecks,exact_score_tilt_MGF_checks=mgfchecks,
        exact_old_prior_Jensen_checks=mgfchecks,exact_exponential_toy=toy,
        actual_normalization_examples_replayed=len(norm['actual_score_tilt_examples']),
        old_prior_spectral_price_checks=spectral_checks,
        exact_physical_score_tilt=physical_full_window_example(),
        changed_prior_example=dict(eta=[str(x)for x in eta],mu=[str(x)for x in mu],tau=str(tau)),
        numerical_replay_count=count,groups=groups,max_centered_trace_difference=maxerr,
        source_SHA256=sha(SRC/'probe.py'),helper42_SHA256=sha(mainprobe.BASE),
        scope='Self-contained bounds and exact finite-law identities; source quadrature replay remains floating exploration, no C_n inference',
        seconds=time.time()-started,numpy_version=np.__version__)
    (HERE/'CHECKS.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
