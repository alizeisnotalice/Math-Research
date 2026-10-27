import json,hashlib,sys,os,time
from pathlib import Path
from fractions import Fraction as Q
sys.set_int_max_str_digits(0);os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
A=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-gaussian106/cube-dim-6-bounds/upper/round138/M5_schur_bridge');P=Path(__file__).resolve().parent
rd=lambda p:json.loads(Path(p).read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
sources=rd(A/'source_manifest.json');cert=rd(A/'candidate_certificate.json')
assert all(sha(r['path'])==r['sha256']and Path(r['path']).stat().st_size==r['bytes']for r in sources)
budget=rd(sources[1]['path']);oldcheck=rd(sources[3]['path']);mu=rd(sources[4]['path']);case=[c for c in budget['cases']if c['alpha']=='1'and c['beta']=='1/2'];assert len(case)==1 and case[0]['b']=='1/10'
JU=Q(case[0]['J_upper_expression']['rational_upper']);assert JU==Q(cert['Jcap_exact'])==Q(2044690556653767,10**15)
Jinterval=oldcheck['independent_512bit_checks']['J'];JL=Q(Jinterval['independent_interval']['lower_fraction']);JH=Q(Jinterval['independent_interval']['upper_fraction']);assert Q(2)<JL<=JH<JU and Jinterval['contained_in_locked_source_ball']is True
writer=Path(sources[1]['path']).parent/'verify_budget.py';assert sha(writer)==budget['source_hashes']['verify_budget.py'];assert any(r['sha256']==sha(writer)and r['path'].endswith('/verify_budget.py')for r in oldcheck['sources'])
# The source writer's float selects a candidate only; actual bracketing is proved by directed Arb comparisons. It is read, not run.
code=writer.read_text();assert 'assert ar(lo)<z and z<ar(hi)'in code and "'rational_upper':str(hi)"in code
lower=JU**5/(2**25*Q(22,7)**2);assert lower==Q(cert['bridge_budget_lower_strict_exact'])and lower>Q(1,10**7)
T0=Q(8,9)**6*24576**5/Q(65000)**6;Ts=Q(8,9)**6*24576**5/Q(63500)**6;M=Q(mu['M_power6_upper']);assert T0==Q(mu['threshold_M_power6_exact'])
comparisons={}
for name,target in [('Theta_g0_1_65000',T0),('Theta_gstar_1_63500',Ts),('saved128_full_Mu_power6',M)]:
 assert lower>target and target<Q(1,10**7)and cert['comparisons'][name]['bridge_budget_lower_strictly_above_target']is True
 comparisons[name]={'target_exact':str(target),'strict_lower_above_target':True,'strict_margin_exact':str(lower-target)}
out={'status':'PASS_C138_SCHUR_PROOF_SOURCE_QUALIFICATION_AND_THREE_EXACT_BUDGET_COMPARISONS','Jcap_exact':str(JU),'Jcap_is_certified_rational_outer_endpoint_not_display':True,'archived_J_expression_lower_exact':str(JL),'archived_J_expression_upper_exact':str(JH),'source_writer_sha256':sha(writer),'source_5D_law':'same16 unnormalized body columns, iid stable Laplace exp(-sqrt(lambda)), underlying A=sum S vv^T','matrix_law_bridge':'A_6block has law4Sigma5; not the Gaussian covariance','J_to_M5':'J=16Gamma(3/2)M5=8sqrt(pi)M5','exact_strict_budget_lower':str(lower),'comparison_count':3,'comparisons':comparisons,'new_Bessel_Gamma_quadrature_optimization_old_baseline_calls':0,'actual_M6_or_C6_lower_bound':False,'all_M5_certificates_or_bridge_family_excluded':False,'C5_operator_upper_transferred':False,'author_certificate_sha256':sha(A/'candidate_certificate.json'),'source_manifest_sha256':sha(A/'source_manifest.json')}
(P/'review_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':out['status'],'comparison_count':3,'Jcap_qualified':True,'strict_budget_lower_above1e7_separator':True,'new_special_function_or_old_baseline_calls':0}))
