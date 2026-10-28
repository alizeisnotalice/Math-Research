#!/usr/bin/env python3
"""Independent bounded R088 constants; original author13 program not executed.

No finite dimension threshold or new parameter choice is evaluated.
The universal algorithmic tau and full operator chain are proved in review.md.
"""
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
import json,os,platform

HERE=Path(__file__).resolve().parent;reg=json.loads((HERE/'REGISTRATION.json').read_text())
AUTHOR=Path('/Users/zhengzhihao/Desktop/T/audit/round088/new_tau')
BASE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-7/cube-dim-7-bounds')
V={k:F(v) for k,v in reg['fixed_tuple'].items()}

def logs(y):
    z=(y-1)/(y+1)
    assert 0<z<1
    lo=2*sum((z**(2*j+1)/F(2*j+1) for j in range(20)),F(0))
    return lo,lo+2*z**41/(41*(1-z*z))

s=V['s'];chi=V['chi'];sig=V['sigma'];d=V['d_tau'];e=V['e'];m=V['m'];r0=1+V['trim'];r1=s*V['qmax']-V['trim']
Qlo,Qhi=logs(2*s/(s+1));mulo=s-1-(s+1)*Qhi;muhi=s-1-(s+1)*Qlo
masslo,_=logs(r1/r0);_,pricehi=logs(1/(1-sig))
P0=1-2*e;ulo=P0*(1-e)/(1+e);uhi=1/ulo;M=2*(1+e)*uhi
b=2*V['Ret']*(1-e)*ulo-2*e;Gamma=m*b/M
meanhi=(muhi+chi)/(2*(1-e));Cfree_lo,_=logs(s/(1+muhi))
checks={
 'independent_fixed_mass_mean_and_deficit': 0<mulo<muhi<chi
    and (1-e)*masslo>m>V['gamma'] and V['Ret']*m>V['gamma']
    and meanhi<V['m_face']<chi and (1-e)*Gamma>V['gamma']
    and 0<e<F(1,12) and M>2 and b>0,
 'new_threshold_interval_and_both_full_cost_rates': 0<chi*(1-sig)<d<sig<1
    and d-chi*(1-sig)==F(1,200) and sig-d==F(1,200),
 'complete_H_not_an_optimized_F_allowance': pricehi<V['price']
    and V['gamma']-V['price']-V['c']==F(1,500)>0,
 'strict_family_contains_point9_without_boundary': Cfree_lo>V['c']
    and 1/s<V['qmax']<1 and 1<r0<r1<s and s/(1-e)<10,
}
paths=[AUTHOR/n for n in ('PROOF.md','REGISTRATION.json','check.py','receipt.json','MANIFEST.json')]
paths += [BASE/'upper/round027/cube_cap_theory/proof.md',BASE/'revision_v3/cover/audit.md',
 BASE/'upper/round085/parameter_family/PROOF.md',BASE/'audit/round085/full/full_review.md',
 BASE/'audit/round085/full/local_uniform.md']
hashes={str(p):sha256(p.read_bytes()).hexdigest() for p in paths}
areg=json.loads((AUTHOR/'REGISTRATION.json').read_text());ar=json.loads((AUTHOR/'receipt.json').read_text())
checks['frozen_proof_and_original_source_receipt_links']=(hashes[reg['source_proof']]==reg['source_proof_sha256']
 and hashes[str(AUTHOR/'REGISTRATION.json')]==ar['registration_sha256']
 and hashes[str(AUTHOR/'check.py')]==ar['check_sha256']
 and all(sha256(Path(p).read_bytes()).hexdigest()==h for p,h in areg['sources'].items()))
assert all(checks.values()),checks
values={'mu_lower':mulo,'mu_upper':muhi,'source_mean_upper':meanhi,
 'source_mean_slack':V['m_face']-meanhi,'source_face_rate':chi-V['m_face'],
 'mass_slack':(1-e)*masslo-m,'paid_gamma_slack':(1-e)*Gamma-V['gamma'],
 'tau_to_k_rate':d-chi*(1-sig),'AQ_over_tau_n_rate':sig-d,
 'price_upper':pricehi,'full_H_allowance':V['gamma']-V['price']-V['c'],
 'Cfree_lower_at_s4p5':Cfree_lo}
out={'status':'PASS independent bounded constants; eventual full proof in review.md','pid':os.getpid(),
 'python':platform.python_version(),'method':'Fraction log20 plus exact rational rates and allowances',
 'no_author13_or_old_baseline_execution':True,'no_dimension_threshold':True,
 'checks':checks,'exact_values':{k:str(v) for k,v in values.items()},
 'display_only':{k:float(v) for k,v in values.items()},'source_sha256':hashes,
 'registration_sha256':sha256((HERE/'REGISTRATION.json').read_bytes()).hexdigest()}
(HERE/'receipt.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ('exact_values','source_sha256')},indent=2))
