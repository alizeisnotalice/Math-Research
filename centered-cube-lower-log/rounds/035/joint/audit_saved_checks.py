"""Read-only independent recomputation of saved CI summaries; no MC."""
from pathlib import Path
import hashlib,json,math

HERE=Path(__file__).resolve().parent
r=json.loads((HERE/'RESULTS.json').read_text())
reg=json.loads((HERE/'registration.json').read_text())
hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest()==reg[k]
        for n,k in [('probe.py','script_SHA256'),('PLAN.md','plan_SHA256'),('calibration.json','calibration_SHA256')]}
hashes['registration.json']=hashlib.sha256((HERE/'registration.json').read_bytes()).hexdigest()==r['registration_SHA256']
assert all(hashes.values())
def ci(mean,cap,n,alpha,upper):
    x=mean/cap
    if upper and x>=1:return cap
    if not upper and x<=0:return 0.
    a,b=(x,1.)if upper else(0.,x)
    target=math.log(1/alpha)/n
    for _ in range(90):
        z=(a+b)/2
        kl=(x*math.log(x/z)if x else 0)+((1-x)*math.log((1-x)/(1-z))if x<1 else 0)
        if upper:
            if kl>target:b=z
            else:a=z
        else:
            if kl>target:a=z
            else:b=z
    return cap*(b if upper else a)
family=0;diff=0.;table=[]
for case in r['cases']:
    row={'n':case['n'],'K':case['K'],'N':case['N'],'J':case['J'],
         'mass_mean':case['mass_weight_mean'],'models':{}}
    for model,v in case['models'].items():
        curve=v['profile'];family+=2*len(curve)
        for q in curve:
            lo=ci(q['coefficient_MC'],q['contribution_bound'],case['N'],r['alpha_each'],False)
            up=ci(q['coefficient_MC'],q['contribution_bound'],case['N'],r['alpha_each'],True)
            up=0. if q['lambda']>=case['density_max'] else min(up,q['lambda'])
            diff=max(diff,abs(lo-q['CI_lower']),abs(up-q['CI_upper']))
        low=max(1.,max(q['CI_lower']for q in curve))
        upper=max(1.,max(min(math.exp(q2['h']),case['density_max'])/q['lambda']*q['CI_upper']
                         for q,q2 in zip(curve,curve[1:])))
        assert low==v['full_weak_lower_CI'] and upper==v['full_weak_upper_CI']
        row['models'][model]={'lower':low,'upper':upper,
                             'best_MC':v['best_frozen_grid_MC']['coefficient_MC'],
                             'best_MC_ESS':v['best_frozen_grid_MC']['tail_ESS']}
    table.append(row)
assert math.isclose(family*r['alpha_each'],.05,rel_tol=1e-14)
out={'scope':'Saved-summary arithmetic audit only; no MC; FP inputs not interval certified.',
     'hash_checks':hashes,'CI_family_count':family,'alpha_sum':family*r['alpha_each'],
     'max_CI_recalculation_absolute_difference':diff,'cases':table,
     'confirmed_joint_weak_gain':any(z['models']['joint']['lower']>
             min(z['models'][m]['upper']for m in ['fixed0','fixed1','prior_mix'])for z in table)}
(HERE/'audit_saved_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
