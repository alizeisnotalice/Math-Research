#!/usr/bin/env python3
"""New free-threshold mechanism: one fixed-tuple exact arithmetic run."""
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
import json, os, sys
P=Path(__file__).resolve().parent
reg=json.loads((P/'REGISTRATION.json').read_text())
assert all(sha256(Path(s).read_bytes()).hexdigest()==h for s,h in reg['sources'].items())
v={k:F(x) for k,x in reg['fixed_tuple'].items()}
def log24(y):
    z=(y-1)/(y+1)
    lo=2*sum((z**(2*j+1)/F(2*j+1) for j in range(24)),F())
    return lo,lo+2*z**49/(49*(1-z*z))
s,e=v['s'],v['e']; rlo=1+v['trim']; rhi=s*v['qmax']-v['trim']
qlo,qhi=log24(2*s/(s+1))
mulo=(s-1)-(s+1)*qhi; muhi=(s-1)-(s+1)*qlo
masslo,_=log24(rhi/rlo)
_,pricehi=log24(1/(1-v['sigma']))
P0=1-2*e; Ulo=P0*(1-e)/(1+e); Uhi=1/Ulo
M=2*(1+e)*Uhi; b=2*v['Ret']*(1-e)*Ulo-2*e
Gamma=v['m']*b/M
meanhi=(muhi+v['chi'])/2/(1-e)
capnu=1/s-(1+1/s)*e
eta=(1-v['qmax'])/(2*v['qmax'])
checks={
 'fixed_geometry':4<=s<=5 and 1/s<v['qmax']<1 and 1<rlo<rhi<s,
 'source_mean_positive_gap':0<mulo<muhi<v['chi'],
 'mass_slack':(1-e)*masslo>v['m']>v['gamma'],
 'retention_choice':v['gamma']/v['m']<v['Ret']<1,
 'inverse_mean_slack':meanhi<v['m_face']<v['chi'],
 'deficit_slack':(1-e)*Gamma>v['gamma'] and b>0,
 'new_primary_sigma':v['sigma']>v['chi']/(1+v['chi']),
 'new_tau_exponents':v['delta_tau']-v['chi']*(1-v['sigma'])==F(1,200) and v['sigma']-v['delta_tau']==F(1,200),
 'tau_positive_and_sublinear':0<v['delta_tau']<v['sigma']<1,
 'weighted_cap_positive':capnu>0 and eta>0,
 'exact_price_upper':pricehi<v['price_upper'],
 'complete_Htau_allowance':v['gamma']-v['price_upper']-v['c']==F(1,500)>0,
 'true_variance_bound':s/(1-e)<10,
}
assert all(checks.values()),checks
vals={'mu_lower':mulo,'mu_upper':muhi,'inverse_face_mean_upper':meanhi,
      'source_face_gap':v['chi']-v['m_face'],
      'mass_lower':(1-e)*masslo,'paid_Gamma':(1-e)*Gamma,
      'tau_to_k_rate':v['delta_tau']-v['chi']*(1-v['sigma']),
      'AB_over_tau_n_rate':v['sigma']-v['delta_tau'],
      'price_upper':pricehi,'complete_allowance_lower':v['gamma']-v['price_upper']-v['c'],
      'cap_nu':capnu,'cap_eta':eta}
assert all(sha256(Path(s).read_bytes()).hexdigest()==h for s,h in reg['sources'].items())
result={'status':'PASS','pid':os.getpid(),'python':sys.executable,'checks':checks,
        'registration_sha256':sha256((P/'REGISTRATION.json').read_bytes()).hexdigest(),
        'check_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
        'exact_values':{k:str(x) for k,x in vals.items()},
        'display_only':{k:float(x) for k,x in vals.items()},
        'source_hashes_rechecked':True,
        'scope':'One fixed tuple; no finite dimension threshold; no original/author baseline replay; Htau only, not a square-root F conclusion.'}
(P/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':'PASS','pid':result['pid'],'checks':len(checks),'new_rates':[str(vals['tau_to_k_rate']),str(vals['AB_over_tau_n_rate'])],'c':str(v['c'])}))
