from fractions import Fraction as F
from pathlib import Path
import json
import hashlib

base=Path(__file__).resolve().parent
prefix='diagonal_levy_spatial_probe_20261007'
regpath=base/(prefix+'_moment_registration.json')
reg=json.loads(regpath.read_text())
resultpath=base/(prefix+'_moment_results.json')
assert not resultpath.exists(), 'Do not overwrite.'
records=[]
for plan in reg['rounds']:
    n,L=plan['n'],F(plan['L'])
    for rstr in reg['r_values']:
        r=F(rstr)
        q=2*r*(1-r)
        rate=q/(2*(1-r))
        assert rate==r
        # Original B(v)=v/2+O(v^2), so P=sqrt(1-rv/2+O(v^2)).
        p_second_frequency_coefficient=-r*L*L/4
        marginal_variance=-2*p_second_frequency_coefficient
        assert marginal_variance==r*L*L/2
        variance_derivative=q*L*L/2+2*marginal_variance
        assert variance_derivative==r*(2-r)*L*L
        normalized_jump_variance=variance_derivative/rate
        assert normalized_jump_variance==(2-r)*L*L<=2*L*L
        records.append({'n':n,'r':rstr,'L':str(L),'rate_per_coordinate':str(rate),
            'rate_total':str(n*rate),'marginal_variance':str(marginal_variance),
            'variance_derivative':str(variance_derivative),
            'normalized_jump_variance':str(normalized_jump_variance),'passed':True})
payload={'status':'passed','exact_cases':len(records),'identity_checks':5*len(records),
    'records':records,'registration_sha256':hashlib.sha256(regpath.read_bytes()).hexdigest(),
    'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'scope':reg['scope'],'claim_limit':reg['claim_limit']}
resultpath.write_text(json.dumps(payload,indent=2)+'\n')
print(json.dumps({'status':payload['status'],'exact_cases':len(records),'identity_checks':payload['identity_checks']}))
