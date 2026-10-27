#!/usr/bin/env python3
"""One new exact budget comparison; no Bessel, optimization or baseline replay."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
import hashlib,json,os,sys,time
HERE=Path(__file__).resolve().parent
FIVE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds')
SIX=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-upper098/cube-dim-6-bounds/upper/round128/true_beta1_joint_moment')

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(x):
    with localcontext() as c:
        c.prec=22
        return str(Decimal(x.numerator)/Decimal(x.denominator))

def main():
    started=time.monotonic()
    assert os.getpriority(os.PRIO_PROCESS,0)>=10
    if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
    sources=[FIVE/'astra_upper014/report.md',FIVE/'astra_upper014/budget.json',
      FIVE/'upper_alpha1_029/report.md',FIVE/'upper_alpha1_029/check_results.json',
      SIX/'candidate_certificate.json',SIX/'run_prepared.py']
    source_manifest=[]
    for p in sources:
        source_manifest.append({'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size})
    budget=json.loads(sources[1].read_text())
    cases=[c for c in budget['cases'] if c['alpha']=='1' and c['beta']=='1/2']
    assert len(cases)==1 and cases[0]['b']=='1/10'
    JU=F(cases[0]['J_upper_expression']['rational_upper'])
    assert JU==F(2044690556653767,10**15)
    # New proof-derived lower bound for this UPPER-BUDGET EXPRESSION.
    # H<=16/min(E_i), hence E H^-1/2>=sqrt(pi)/32.
    # B6=2^-5*(JU/(8sqrt(pi)))^5*E H^-1/2
    # >=JU^5/(2^25*pi^2)>JU^5/(2^25*(22/7)^2).
    lower=JU**5/F(2**25)/F(22,7)**2
    source128=json.loads(sources[4].read_text())
    theta0=F(8,9)**6*24576**5/65000**6
    theta_star=F(8,9)**6*24576**5/63500**6
    assert theta0==F(source128['threshold_M_power6_exact'])
    mu128=F(source128['M_power6_upper'])
    targets={'Theta_g0_1_65000':theta0,'Theta_gstar_1_63500':theta_star,'saved128_full_Mu_power6':mu128}
    separator=F(1,10**7)
    assert lower>separator
    comparisons={}
    for name,t in targets.items():
        passed=lower>t
        assert passed and t<separator
        comparisons[name]={'bridge_budget_lower_strictly_above_target':passed,
          'target_strictly_below_separator':True,'target_display_only':dec(t),
          'bridge_lower_over_target_display_only':dec(lower/t)}
    certificate={'status':'PASS_current_M5_certificate_plus_Holder_bridge_budget_insufficient',
      'authorization':'parent message explicitly authorized this one standard-library Fraction batch with three comparisons',
      'Jcap_exact':str(JU),'M5_upper_expression':'Jcap/(8*sqrt(pi))',
      'source_law_match':'V5={(1,+/-1,+/-1,+/-1,+/-1)}; iid Laplace exp(-sqrt(lambda)); A=sum S_v vv^T; S=1/(4Y), Y~Gamma(1/2,rate1)',
      'J_to_M5_identity':'J=16*Gamma(3/2)*M5=8*sqrt(pi)*M5',
      'B_definition':'E[(sum_(j=1)^16 1/E_j)^(-1/2)], E_j iid Exp(rate1)',
      'B_lower_proof':'H<=16/min(E_j); min(E_j)~Exp(rate16); E sqrt(min(E_j))/4=sqrt(pi)/32',
      'bridge_upper_budget_power6':'Jcap^5*B/(2^20*pi^(5/2))',
      'bridge_budget_lower_strict_exact':str(lower),
      'bridge_budget_lower_display_only':dec(lower),
      'strict_lower_justification':'B>=sqrt(pi)/32 and pi<22/7, so bridge budget power6>Jcap^5/(2^25*(22/7)^2)',
      'shared_rational_separator':str(separator),
      'bridge_budget_lower_strictly_above_separator':True,
      'comparisons':comparisons,'arithmetic':'Python stdlib Fraction, all decisions exact; Decimal display only',
      'scientific_batch_count':1,'comparison_count':3,'processes':1,
      'nice':os.getpriority(os.PRIO_PROCESS,0),'wall_cap':None,
      'Bessel_evaluations':0,'optimization_calls':0,'old_baseline_reruns':0,
      'actual_M6_lower_bound':False,'actual_C6_lower_bound':False,
      'all_possible_M5_certificates_excluded':False,'bridge_family_globally_excluded':False,
      'global_g0_proved':False,'global_gstar_proved':False,'new_C6_bound':False,
      'script_sha256':sha(Path(__file__)),'sources':source_manifest,
      'seconds':time.monotonic()-started}
    (HERE/'candidate_certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
    (HERE/'source_manifest.json').write_text(json.dumps(source_manifest,indent=2)+'\n')
    print(json.dumps({k:certificate[k] for k in ['status','bridge_budget_lower_display_only',
      'shared_rational_separator','comparisons','seconds']},indent=2))
if __name__=='__main__':main()
