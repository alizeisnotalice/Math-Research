#!/usr/bin/env python3
"""Exact signed overlap and exchange-price guard, not actual FIRST."""
from fractions import Fraction as F
from pathlib import Path
from math import prod
import hashlib
import json
HERE=Path(__file__).resolve().parent
PREFIX='late_multiface_source_budget'
checks=[]
rounds=[]

def rec(label, ok):
    checks.append({'label':label,'pass':bool(ok)})

def Q(u,d):
    return [sum((2*u[x]+u[x^(1<<i)])/3 for i in range(d)) for x in range(len(u))]

def transform(u):
    N=len(u)
    return [sum(a*(-1 if (x&m).bit_count()%2 else 1) for x,a in enumerate(u)) for m in range(N)]

def inverse(v):
    N=len(v)
    return [sum(a*(-1 if (x&m).bit_count()%2 else 1) for m,a in enumerate(v))/N for x in range(N)]

def resolvent(u,d,t):
    v=transform(u)
    return inverse([a/(1+t*F(2*m.bit_count(),3)) for m,a in enumerate(v)])

def ihash(value):
    return hashlib.sha256(value.to_bytes((value.bit_length()+7)//8,'big')).hexdigest()

for d,n,j in [(3,512,64),(4,32768,1024),(5,2097152,16384)]:
    start=len(checks)
    u=[F(2+x%7,1+x%3) if x%4==0 or x==1 else F(0) for x in range(1<<d)]
    qu=Q(u,d)
    sigma=[d*a-b for a,b in zip(u,qu)]
    kappa=max(-a for a in sigma)
    nu=[kappa+sigma[x] if u[x]>0 else F(0) for x in range(len(u))]
    mu=[nu[x]-sigma[x] for x in range(len(u))]
    W=sum(nu)
    tag=f'n{n}/finite_d{d}'
    rec(tag+'/original-cap', all(F(0)<=a<=kappa for a in mu))
    rec(tag+'/original-saturation', all(mu[x]==kappa for x,a in enumerate(u) if a>0))
    rec(tag+'/source-mass', sum(mu)==W and all(a>=0 for a in nu))
    L=(n+1).bit_length()
    receipts=[]
    for t in [F(1),F(L)]:
        rv=resolvent(nu,d,t)
        rm=resolvent(mu,d,t)
        h=[a-b for a,b in zip(rv,rm)]
        hp=[max(a,F(0)) for a in h]
        hm=[max(-a,F(0)) for a in h]
        b=[min(a,c) for a,c in zip(rv,rm)]
        tt=tag+f'/t{t.numerator}'
        rec(tt+'/overlap-positive-cap', all(F(0)<=a<=kappa for a in b))
        rec(tt+'/overlap-source-mass', sum(b)<=W)
        rec(tt+'/original-exterior-holding', all(hp[x]==0 for x,a in enumerate(u) if a==0))
        rec(tt+'/overlap-identity', all(b[x]==rm[x]-hm[x] and hp[x]==rv[x]-b[x] for x in range(len(u))))
        rec(tt+'/positive-source-mass', sum(hp)<=W)
        hph,nh,bh=transform(hp),transform(nu),transform(b)
        for r in [F(0),F(1,3),F(2,3),F(1)]:
            for m in range(1<<d):
                s=F(2*m.bit_count(),3)
                km=(1-F(2,3)*r)**m.bit_count()
                # Exact coefficient pair (constant, common exp(-r*s)).
                lhs=(km*(1+t*s)*hph[m],-(1+t*s)*hph[m])
                rhs=(km*nh[m]-km*(1+t*s)*bh[m],-nh[m]+(1+t*s)*bh[m])
                rec(tt+f'/r{r}/mode{m}', lhs==rhs)
        C=23*t+1+(2*t+F(1,n))/(n-1)
        rec(tt+'/signed-one-face-total', C<=26*t+1)
        cutoff=12
        partial=n*sum((2*t+F(1,n))/F(n**(m-1)) for m in range(3,cutoff+1))
        tail=n*(2*t+F(1,n))/F(n**cutoff)/(1-F(1,n))
        rec(tt+'/signed-one-face-infinite-tail', partial+tail==(2*t+F(1,n))/(n-1))
        rec(tt+'/degree2-ledger', 2*t/n+2*(1+t*n)/(n*n)<=5*t/n)
        receipts.append({'t':str(t),'W_b':str(W),'kappa':str(kappa),'overlap_mass':str(sum(b)),'one_face_C_t':str(C)})
    numerator=n**j
    denominator=prod(range(n-j+1,n+1))
    exponent=j*(j-1)//(2*n)
    rec(tag+'/without-replacement-loss', numerator>=denominator*(1<<exponent))
    rec(tag+'/late-degree-scope', j*j>n and j<n)
    rec(tag+'/direct-D-small-r-constant', F(873,64)<14)
    rec(tag+'/direct-D-large-r-constant', n*n<=(1<<(n+1)))
    rec(tag+'/direct-D-source-envelope', n*F(14,n*n)==F(14,n))
    rec(tag+'/direct-D-ledger', 4*F(14,n)==F(56,n))
    rounds.append({'finite_algebra_dimension':d,'compressed_dimension':n,'degree':j,'checks':len(checks)-start,
                   'source_receipts':receipts,'exchange':{'lower_power_of_two':exponent,'numerator_bits':numerator.bit_length(),
                   'denominator_bits':denominator.bit_length(),'numerator_sha256':ihash(numerator),'denominator_sha256':ihash(denominator)},
                   'scope':'Finite algebra / compressed exchange price only; not Rn obstacle samples.'})
registration=HERE/f'{PREFIX}_registration_20261007.json'
result={'status':'PASS' if all(c['pass'] for c in checks) else 'FAIL','seed':None,'total_checks':len(checks),
        'failed_checks':[c for c in checks if not c['pass']],'rounds':rounds,'checks':checks,
        'scope':'New signed overlap-source identities and full-one-face fee; compressed Q^j comparison loss. Not actual FIRST or a weak counterexample.',
        'registration_sha256':hashlib.sha256(registration.read_bytes()).hexdigest(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(HERE/f'{PREFIX}_results_20261007.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['status','total_checks','failed_checks','script_sha256']},ensure_ascii=False))
raise SystemExit(0 if result['status']=='PASS' else 1)
