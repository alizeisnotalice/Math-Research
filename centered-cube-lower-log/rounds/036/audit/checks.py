"""Exact small Dirichlet checks and saved-statistics audit only; no MC."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
cell=ROOT/'mod';theory=ROOT/'theory'
hashes={}
for script,registration,result in [('probe.py','registration.json','results.json'),('short_probe.py','short_registration.json','short_results.json')]:
    reg=json.loads((theory/registration).read_text());res=json.loads((theory/result).read_text())
    hashes['theory/'+script]=hashlib.sha256((theory/script).read_bytes()).hexdigest()==reg['script_sha256']
    hashes['theory/'+registration]=hashlib.sha256((theory/registration).read_bytes()).hexdigest()==res['registration_sha256']
reg=json.loads((cell/'registration.json').read_text());res=json.loads((cell/'RESULTS.json').read_text())
for name,key in [('probe.py','script_SHA256'),('PLAN.md','plan_SHA256'),('calibration.json','calibration_SHA256')]:
    hashes['cell/'+name]=hashlib.sha256((cell/name).read_bytes()).hexdigest()==reg[key]
hashes['cell/registration.json']=hashlib.sha256((cell/'registration.json').read_bytes()).hexdigest()==res['registration_SHA256']
assert all(hashes.values())

# Gaussian-rational roots of unity: squared Q avoids sqrt/roundoff.
count=0
roots={2:[(1,0),(-1,0)],4:[(1,0),(0,1),(-1,0),(0,-1)]}
for R,rr in roots.items():
    for N in range(2,65):
        for h in range(1,R):
            re=sum(rr[(h*j)%R][0]for j in range(N));im=sum(rr[(h*j)%R][1]for j in range(N))
            q2=F(re*re+im*im,N*N)
            c=min(F(1),F(2*N,R))
            assert q2<=1-c*c/4
            count+=1
rational=json.loads((theory/'rational_results.json').read_text())
for row in rational['rows']:
    relative=abs(F(row['joint'])/F(row['prior'])-1)
    assert relative==F(row['relative_difference'])<=F(row['universal_candidate'])

def ci(mean,cap,n,alpha,upper):
    x=mean/cap
    if upper and x>=1:return cap
    if not upper and x<=0:return 0.
    lo,hi=(x,1.)if upper else(0.,x);target=math.log(1/alpha)/n
    for _ in range(90):
        m=(lo+hi)/2
        kl=(x*math.log(x/m)if x else 0)+((1-x)*math.log((1-x)/(1-m))if x<1 else 0)
        if upper:
            if kl>target:hi=m
            else:lo=m
        else:
            if kl>target:lo=m
            else:hi=m
    return cap*(hi if upper else lo)
family=0;maxdiff=0.;rows=[]
for case in res['cases']:
    for model,value in case['models'].items():
        q=value['holdout'];family+=2
        for upper,key in [(False,'CI_lower'),(True,'CI_upper')]:
            calc=ci(q['coefficient'],q['contribution_cap'],case['N']//2,res['alpha_each'],upper)
            maxdiff=max(maxdiff,abs(calc-q[key]))
    j=case['models']['joint'];p=case['models']['prior']
    rows.append({'n':case['n'],'R':case['R'],'joint_unrestricted':j['empirical_unrestricted'],
                 'joint_ESS64':j['empirical_ESS_ge64'],
                 'joint_holdout':j['holdout'],'prior_holdout':p['holdout']})
assert family==138 and math.isclose(family*res['alpha_each'],.05,rel_tol=1e-14)
out={'scope':'No MC; exact finite Dirichlet checks plus saved summary CI recomputation.',
     'hash_checks':hashes,'exact_Dirichlet_checks':count,'rational_macro_checks':len(rational['rows']),
     'holdout_CI_family_count':family,'alpha_sum':family*res['alpha_each'],
     'max_CI_recalculation_difference':maxdiff,'cases':rows}
(HERE/'checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items()if k!='cases'},indent=2))
for x in rows:
    q=x['joint_unrestricted'];s=x['joint_ESS64'];v=x['joint_holdout']
    print(x['n'],x['R'],'unrestricted',q['coefficient'],q['ESS'],'ESS64',s['coefficient'],
          'holdout',v['coefficient'],[v['CI_lower'],v['CI_upper']])
