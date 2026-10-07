#!/usr/bin/env python3
"""New rare internal-clock selector certificate. Not a physical weak endpoint."""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import sys,json,hashlib
sys.set_int_max_str_digits(0)
BASE=Path(__file__).resolve().parent
REG=BASE/'high_support_source_selector_registration_20261007.json'
OUT=BASE/'high_support_source_selector_results_20261007.json'
registration=json.loads(REG.read_text())
assert registration['status']=='registered_before_execution'
checks=0
def check(value):
    global checks
    assert value
    checks+=1
def fq(x):
    return {'numerator':str(x.numerator),'denominator':str(x.denominator)}
def exp_negative(x,m):
    """Positive Taylor lower/upper for exp(x), then exact reciprocal."""
    term=Q(1); total=term
    for j in range(1,m+1):
        term=term*x/j
        total+=term
    ratio=x/Q(m+2)
    check(ratio<1)
    nextterm=term*x/Q(m+1)
    upper=total+nextterm/(1-ratio)
    return 1/upper,1/total
def gamma_tail(k,x,m,lower):
    el,eu=exp_negative(x,m)
    term=Q(1); series=term
    for j in range(1,k):
        term=term*x/j
        series+=term
    if lower:
        lo=max(Q(0),1-eu*series);hi=min(Q(1),1-el*series)
    else:
        lo=max(Q(0),el*series);hi=min(Q(1),eu*series)
    check(0<=lo<=hi<=1)
    return lo,hi
rounds=[]
for spec in registration['rounds']:
    n,K,m=spec['n'],spec['K'],spec['exp_terms']
    before=checks
    check(K**3<=n<(K+1)**3)
    check(K>=2)
    C=3*Q(2,3)**(K+1)+Q(3,2)*Q(1,3)**(K+1)
    check(C<=Q(17,18))
    rows=[];masslo=Q(0);masshi=Q(0);qsumhi=Q(0)
    for k in range(K+1,n+1):
        ll,lu=gamma_tail(k,Q(k,4),m,True)
        ul,uu=gamma_tail(k,Q(4*k),m,False)
        check(lu<=Q(2,3)**k)
        check(uu<=Q(1,3)**k)
        p=Q(k,n)**k*Q(n-k,n)**(n-k)
        binpeak=comb(n,k)*p
        check(0<binpeak<=1)
        check(0<=ll+ul<=lu+uu<=1)
        masslo+=binpeak*(ll+ul)
        masshi+=binpeak*(lu+uu)
        qsumhi+=lu+uu
        rows.append({'k':k,'gamma_lower_tail_interval':[fq(ll),fq(lu)],
          'gamma_upper_tail_interval':[fq(ul),fq(uu)],
          'binomial_peak_mass':fq(binpeak)})
    check(masslo<=masshi<=qsumhi<=C)
    # Exact Chernoff constants: e^(1/4)<=4/3 and e^2>6.
    lo,hi=exp_negative(Q(1,4),m)
    check(1/lo<=Q(4,3))
    lo2,hi2=exp_negative(Q(2),m)
    check(1/hi2>6)
    rounds.append({'n':n,'K':K,'exp_terms':m,'checks':checks-before,
       'C_Gamma':fq(C),'envelope_mass_interval':[fq(masslo),fq(masshi)],
       'sum_gamma_tail_upper':fq(qsumhi),'rows':rows})
    print('round',n,'PASS',checks-before,'checks',flush=True)
result={'status':'PASS_EXACT_RATIONAL_INTERVAL','total_checks':checks,
 'registration_sha256':hashlib.sha256(REG.read_bytes()).hexdigest(),
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'rounds':rounds,
 'scope':'New original Exp/Gamma latent-law rare-sector budget. Uniform c follows analytic Markov pushforward, not a sampled spatial reconstruction. Regular-clock selector remains unpaid.'}
OUT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print('TOTAL',checks,'PASS',flush=True)

