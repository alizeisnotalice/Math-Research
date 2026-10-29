"""Deterministic checks of KL monotonicity and method-uniform target coverage."""
from pathlib import Path
import itertools,json,math
HERE=Path(__file__).resolve().parent
def ci(x,n,alpha,upper):
    beta=math.log(1/alpha)/n
    if x==0:return 1-math.exp(-beta)if upper else 0.
    if x==1:return 1. if upper else math.exp(-beta)
    a,b=(x,1.)if upper else(0.,x)
    for _ in range(90):
        p=(a+b)/2
        # Extreme n1 roots can be closer to 0/1 than binary64 resolves.
        # Retain the conservative boundary instead of taking log(0).
        if p<=0:return 0.
        if p>=1:return 1.
        kl=x*math.log(x/p)+(1-x)*math.log((1-x)/(1-p))
        if upper:
            if kl>beta:b=p
            else:a=p
        else:
            if kl>beta:a=p
            else:b=p
    return(a+b)/2
grid_checks=0
for n in [1,8,256]:
    for alpha in [.05,.05/504,.05/1008]:
        for upper in [False,True]:
            vals=[ci(i/64,n,alpha,upper)for i in range(65)]
            assert all(a<=b+1e-15 for a,b in zip(vals,vals[1:]));grid_checks+=64
adaptive=0
for sample in itertools.product([0.,.25,.75,1.],repeat=4):
    m=sum(sample)/4
    # Whole-sample choices; candidate bound observations need not be iid.
    offset=(1-m)/3
    lo=sum(max(0,x-offset)for x in sample)/4
    up=sum(min(1,x+m/3)for x in sample)/4
    for alpha in [.05,.05/1008]:
        assert ci(lo,4,alpha,False)<=ci(m,4,alpha,False)+1e-15
        assert ci(up,4,alpha,True)>=ci(m,4,alpha,True)-1e-15
        adaptive+=1
# Data-tuned cap need not bound the population: n1 Bernoulli(.5), cap=max(Y,1e-6).
tuned_cap_upper_at_zero=1e-6*ci(0.,1,.05,True)
assert tuned_cap_upper_at_zero<.5
# Data-tuned new target: Uniform[100], choose indicator of observed category.
chosen_target_lower=ci(1.,1,.05,False)
assert chosen_target_lower>1/100
out={'monotonic_grid_neighbor_checks':grid_checks,'adaptive_samplewise_sandwich_checks':adaptive,
     'data_tuned_cap_counterexample':{'n':1,'target_mean':.5,'alpha':.05,
        'upper_endpoint_when_Y0':tuned_cap_upper_at_zero,'failure_probability_at_least':.5},
     'data_tuned_target_counterexample':{'categories':100,'n':1,'selected_indicator_mean':.01,
        'alpha':.05,'selected_lower_endpoint':chosen_target_lower,'failure_probability':1.},
     'scope':'Finite exhaustive calculations, no random samples or source replay'}
p=HERE/'method_uniform_ci_checks.json'
if p.exists():raise RuntimeError('Refusing replacement')
p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
