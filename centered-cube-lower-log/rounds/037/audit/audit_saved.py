"""Saved quantile priors/hashes/heldout CI audit; no MC rerun."""
from pathlib import Path
from fractions import Fraction as F
import json,math,hashlib,sys

HERE=Path(__file__).resolve().parent;SRC=HERE.parent/'cell'
input_name=sys.argv[1]if len(sys.argv)>1 else 'RESULTS_stable_ESS.json'
g=json.loads((SRC/'registration.json').read_text());r=json.loads((SRC/input_name).read_text())
hashes={name:hashlib.sha256((SRC/name).read_bytes()).hexdigest()==g[key]
        for name,key in [('probe.py','script_SHA256'),('PLAN.md','plan_SHA256'),('calibration.json','calibration_SHA256')]}
hashes['registration.json']=hashlib.sha256((SRC/'registration.json').read_bytes()).hexdigest()==r['registration_SHA256']
assert all(hashes.values())
def ci(mean,cap,N,alpha,upper):
    x=mean/cap
    if upper and x>=1:return cap
    if not upper and x<=0:return 0.
    lo,hi=(x,1.)if upper else(0.,x);target=math.log(1/alpha)/N
    for _ in range(90):
        z=(lo+hi)/2
        div=(x*math.log(x/z)if x else 0)+((1-x)*math.log((1-x)/(1-z))if x<1 else 0)
        if upper:
            if div>target:hi=z
            else:lo=z
        else:
            if div>target:lo=z
            else:hi=z
    return cap*(hi if upper else lo)
family=0;diff=0.;rows=[]
for case in r['cases']:
    meta=case['metadata'];n,R=meta['n'],meta['R'];den=1<<n
    nums=[math.comb(n,s)for s in range(n+1)]
    cumulative=[];v=0
    for x in nums:v+=x;cumulative.append(v)
    cuts=[-1]+[next(s for s,x in enumerate(cumulative)if x*R>=b*den)for b in range(1,R)]+[n]
    priors=[F(sum(nums[cuts[b]+1:cuts[b+1]+1]),den)for b in range(R)]
    assert sum(priors,F(0))==1
    for b,segment in enumerate(meta['segments']):
        assert segment['lower']==cuts[b]+1 and segment['upper']==cuts[b+1]
        assert F(segment['prior_exact'])==priors[b]
    for name,model in case['models'].items():
        q=model['holdout'];family+=2
        assert q['log_threshold']<model['discovery']['log_threshold']
        for upper,key in [(False,'CI_lower'),(True,'CI_upper')]:
            calc=ci(q['coefficient'],q['contribution_cap'],case['N']//2,r['alpha_each'],upper)
            diff=max(diff,abs(calc-q[key]))
    joint=case['models']['joint']
    rows.append({'n':n,'R':R,'N':case['N'],'J':meta['J'],
                 'prior_exact':list(map(str,priors)),
                 'empty_classes':[b for b,x in enumerate(priors)if x==0],
                 'joint_unrestricted':joint['empirical_unrestricted'],
                 'joint_ESS64':joint['empirical_ESS_ge64'],'joint_holdout':joint['holdout'],
                 'prior_holdout':case['models']['prior']['holdout'],
                 'mass_mean':case['mass_weight_mean'],
                 'PB_mass_residual':case['max_PB_probability_mass_residual']})
assert family==120 and math.isclose(family*r['alpha_each'],.05,rel_tol=1e-14)
out={'scope':'Read saved statistics only; no source/tail MC rerun; no interval FP certification.',
     'input_results':input_name,
     'hash_checks':hashes,'exact_quantile_prior_cases':len(rows),
     'heldout_CI_family_count':family,'alpha_sum':family*r['alpha_each'],
     'max_CI_recalculation_difference':diff,'cases':rows}
(HERE/('audit_stable_results.json'if input_name=='RESULTS_stable_ESS.json'else 'audit_saved_results.json')).write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items()if k!='cases'},indent=2))
