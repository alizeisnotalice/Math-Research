"""Small deterministic/Fraction checks only; no round40 payload transformation or MC."""
from pathlib import Path
from fractions import Fraction as F
import itertools,json,math
HERE=Path(__file__).resolve().parent
def phi(b):return 3*b/(1+2*b)
frac=0
for factor in [F(1,100),F(1,2),F(1),F(10)]:
    for b in [F(0),F(1,100),F(1,4),F(1,2),F(3,4),F(1)]:
        x=factor*b;recovered=x/factor
        assert recovered==b
        transformed=factor*phi(recovered)
        direct=3*x/(1+2*x/factor)
        assert transformed==direct and x<=transformed<=factor
        old_upper=factor*(phi(b)+1)/2
        assert transformed<=old_upper
        frac+=1
for vals in itertools.product([F(0),F(1,10),F(1,2),F(1)],repeat=3):
    assert phi(max(vals))==max(map(phi,vals))
FP=[]
for h in [0.,10.]:
    for logD in [-math.log(2),0.,100.,700.,800.]:
        factor=math.exp(h-logD)
        for b in [0.,1e-200,1e-20,.01,.5,1.]:
            x=factor*b
            if not x:
                output=0.;recovered=0.
            else:
                logb=math.log(x)+logD-h
                recovered=math.exp(min(0.,logb))
                output=3*x/(1+2*recovered)
            FP.append({'h':h,'logD':logD,'B_generated':b,'stored_x':x,
                       'B_recovered_clipped_descriptive':recovered,'new_x_descriptive':output})
def endpoint(x,n,alpha,upper):
    target=math.log(1/alpha)/n
    if(upper and x==1)or(not upper and x==0):return 1. if upper else 0.
    a,b=(x,1.)if upper else(0.,x)
    for _ in range(80):
        p=(a+b)/2
        kl=(x*math.log(x/p)if x else 0)+((1-x)*math.log((1-x)/(1-p))if x<1 else 0)
        if upper:
            if kl>target:b=p
            else:a=p
        else:
            if kl>target:a=p
            else:b=p
    return(a+b)/2
# Illustrate sample-wise domination with deliberately data-dependent valid bounds.
n=8;events=0;true_fail=0;bound_fail=0
for sample in itertools.product([0.,1.],repeat=n):
    m=sum(sample)/n
    # Choose the more favorable of several bounds after seeing all observations.
    lowers=[[a*x for x in sample]for a in [.2,.8,1.]]
    uppers=[[x+a*(1-x)for x in sample]for a in [0.,.2,.8]]
    ml=max(sum(a)/n for a in lowers);mu=min(sum(a)/n for a in uppers)
    assert endpoint(ml,n,.05,False)<=endpoint(m,n,.05,False)+1e-15
    assert endpoint(mu,n,.05,True)>=endpoint(m,n,.05,True)-1e-15
    tf=endpoint(m,n,.05,False)>.5 or endpoint(m,n,.05,True)<.5
    bf=endpoint(ml,n,.05,False)>.5 or endpoint(mu,n,.05,True)<.5
    assert not bf or tf
    true_fail+=tf;bound_fail+=bf;events+=1
out={'Fraction_recovery_and_phi_cases':frac,'max_exchange_cases':64,'FP_stress_cases':FP,
     'FP_qualification':'Clipping/stress observations are not directed-rounding certification; zero cannot recover lost subnormal information.',
     'adaptive_bound_domination_toy':{'n':n,'enumerated_samples':events,
       'two_direction_true_endpoint_failure_probability':true_fail/events,
       'bound_endpoint_failure_probability':bound_fail/events},
     'family':{'targets':4*3*21,'old_one_sided':504,'new_lower':252,'reserved_new_upper':252,
       'reserved_total':1008,'alpha_each':.05/1008,'lower_stage_functions':756,
       'lower_stage_union_budget':756*.05/1008,'unused_reserved_budget':252*.05/1008},
     'scope':'No Gaussian/probability approximation, no old payload transformed, no MC'}
target=HERE/'postprocess_checks.json'
if target.exists():raise RuntimeError('Refusing replacement')
target.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items()if k!='FP_stress_cases'},indent=2))
