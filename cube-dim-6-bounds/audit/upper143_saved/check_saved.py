#!/usr/bin/env python3
"""C143 independent finite Fraction bookkeeping, no functional replay."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
import gzip,hashlib,json,os,sys,time,traceback
sys.set_int_max_str_digits(0);sys.dont_write_bytecode=True
P=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-gaussian106/cube-dim-6-bounds/upper/round143/adaptive_second_count_budget')
S=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round130/true_second_count_boxes')
T=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round128/true_beta1_joint_moment')
A=Path(__file__).resolve().parent
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(n,d):
 with (A/n).open('x')as f:json.dump(d,f,indent=2);f.write('\n')
def display(q):
 with localcontext()as c:c.prec=30;return str(Decimal(q.numerator)/Decimal(q.denominator))
budget=read(P/'saved_budget_analysis.json');ob=read(P/'obstacle_certificate.json');old=read(S/'results.json')
files=[P/n for n in ['inspect_saved_budget.py','check_saved_obstacles.py','saved_budget_analysis.json','obstacle_certificate.json','actual_fraction_receipt.json','proof.md','report.md','source_manifest.json']]+[S/n for n in ['proof.md','run_prepared.py','inputs.json','results.json','cone_table.json']]+[T/n for n in ['proof.md','run_prepared.py','joint_core.py','upstream.json','candidate.json','candidate_certificate.json']]+[S/v['path']for v in old['leaf_files']]
pins={str(p):sha(p)for p in files}
write('finite_registration.json',{'status':'REGISTERED_BEFORE_FRACTION_RECHECK','payload_hashes':pins,'finite_scope':{'saved_C130_cones':21,'saved_C130_leaf_integral_intervals':5376,'saved_candidate_ray':1,'replacement_cap_coordinates':1,'prefix_H_certificate_scalar_records':1},'new_functional_evaluations':0,'new_Arb_log_exp_root_native_optimizer_calls':0,'old_functional_replays':0,'time_limit':None,'reviewer_script_sha256':sha(__file__)})
os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)));start=time.monotonic();cpu=time.process_time()
try:
 for row in budget['source_hashes']:assert sha(row['path'])==row['sha256'] and Path(row['path']).stat().st_size==row['bytes']
 for row in ob['source_hashes']:assert sha(row['path'])==row['sha256']
 for row in read(P/'source_manifest.json')['sources']:assert sha(row['path'])==row['sha256'] and Path(row['path']).stat().st_size==row['bytes']
 assert sha(P/'inspect_saved_budget.py')==budget['script_sha256'] and sha(P/'check_saved_obstacles.py')==ob['script_sha256']
 authorreceipt=read(P/'actual_fraction_receipt.json');assert authorreceipt['status']=='PASS' and authorreceipt['exit_code']==0 and authorreceipt['child_reaped'] and authorreceipt['script_sha256']==ob['script_sha256'] and authorreceipt['new_functional_erfc_exp_log_root_optimizer_calls']==0 and authorreceipt['time_limit']is None
 assert old['status']=='COMPLETE_FIXED_FINITE_BOX_PILOT' and old['h']=='1/128' and not old['probability_renormalized']
 cones=read(S/'cone_table.json')['cones'];spec=read(S/'inputs.json');up=read(T/'upstream.json');cert=read(T/'candidate_certificate.json');cand=read(T/'candidate.json')
 beta=F(up['beta_by_h'][5][1]);assert F(up['h_values'][5])==F(spec['h'])==F(1,128) and F(old['old_beta2_interval'][1])==F(spec['old_beta2_interval'][1])==beta>0
 muL,muU=map(F,old['moment_interval']);factorL,factorU=map(F,old['pi_factor_interval']);assert 0<factorL<=factorU and 0<muL<=muU
 lows=F(0);upps=F(0);widths=[];ranked=[];vertex_count=0;delta=F(0);cut=1-F(1,2**24)
 assert len(cones)==len(old['leaf_files'])==21
 for cone,entry in zip(cones,old['leaf_files']):
  path=S/entry['path'];assert sha(path)==entry['sha256'];data=json.load(gzip.open(path,'rt'));assert data['cone_id']==cone['id'] and data['cone_table_sha256']==sha(S/'cone_table.json') and len(data['leaves'])==256
  seen={};alpha=F(cone['alpha']);assert alpha>0;prefixes=[v['prefix']for v in data['leaves']];assert len(prefixes)==len(set(prefixes)) and sum((F(1,2**len(p))for p in prefixes),F(0))==1
  assert not any(a!=b and b.startswith(a)for a in prefixes for b in prefixes)
  for leaf in data['leaves']:
   L,U=map(F,leaf['integral_interval']);assert 0<=L<=U
   lo=list(map(F,leaf['lower_t']));hi=list(map(F,leaf['upper_t']));v=F(1)
   for a,b in zip(lo,hi):assert 0<=a<b<=cut;v*=b-a
   assert v==F(leaf['volume']);w=factorU*alpha*(U-L);widths.append(w);lows+=alpha*L;upps+=alpha*U
   for side in ('lower','upper'):
    vertex=tuple(leaf[side+'_t']);vector=leaf['kernel_at_'+side+'_t'];assert len(vector)==32
    digest=hashlib.sha256(json.dumps(vector,separators=(',',':')).encode()).hexdigest()
    if vertex in seen:assert seen[vertex]==digest
    else:
     seen[vertex]=digest
     for l,u in vector:
      l,u=F(l),F(u);assert 0<=l<=u<=F(spec['p_interval'][1]);delta=max(delta,u-l)
   ranked.append({'cone_id':cone['id'],'prefix':leaf['prefix'],'weighted_width':str(w),'integral_lower':str(L),'integral_upper':str(U),'lower_t':leaf['lower_t'],'upper_t':leaf['upper_t']})
  assert len(seen)<2048;vertex_count+=len(seen)
 assert len(widths)==5376 and vertex_count==10691 and 32*vertex_count==342112 and old['kernel_calls']==343456 and old['kernel_calls']-32*vertex_count==21*64
 uncertainty=factorU*sum((F(c['alpha'])for c in cones),F(0))*cut**6*1984*F(spec['p_interval'][1])*delta
 assert F(budget['global_endpoint_uncertainty_upper']['exact'])==uncertainty and F(budget['max_saved_kernel_interval_width']['exact'])==delta
 tailprob=1-cut**6;pU=F(spec['p_interval'][1]);tail_per_cone=496*pU**2*tailprob
 tail_upper=factorU*sum((F(c['alpha'])for c in cones),F(0))*tail_per_cone
 assert factorL*lows==muL and factorU*upps+tail_upper==muU
 width=sum(widths);assert width==F(old['weighted_box_width']) and tail_upper==F(old['tail_upper']) and (factorU-factorL)*lows==F(old['pi_factor_width'])
 assert muU-muL==width+tail_upper+F(old['pi_factor_width'])
 ranked.sort(key=lambda r:(-F(r['weighted_width']),r['cone_id'],r['prefix']))
 widths.sort(reverse=True);gap=muU-beta;assert gap>0;s15=sum(widths[:15]);s2188=sum(widths[:2188]);s2189=s2188+widths[2188]
 assert s15<=s2188<=gap<s2189 and s2188<gap
 patch=ob['patch_obstacle'];assert F(patch['old_upper_minus_beta2_exact'])==gap and F(patch['top15_recoverable_width_exact'])==s15 and F(patch['top2188_recoverable_width_exact'])==s2188 and F(patch['top2189_recoverable_width_exact'])==s2189 and patch['necessary_minimum_patched_old_leaves']==2189
 ratio=gap/width;assert F(budget['required_fraction_of_box_width_recovery']['exact'])==ratio and budget['necessary_minimum_old_leaves_to_patch_even_with_ideal_width_recovery']==2189
 # Prefix lower is the positive subprobability expectation lower; never an M lower.
 HL,HU=map(F,cert['H_interval']);assert 0<HL<=HU
 assert HL==F(int(cert['prefix_lower_at_pU_numerator']),2**cert['weighted_sum_denominator_power2'])
 assert HU==F(int(cert['prefix_upper_at_pL_numerator'])+int(cert['tail_upper_numerator']),2**cert['weighted_sum_denominator_power2'])
 HL5=HL**5;HU5=HU**5;assert HL5==F(cert['prefix_only_H_power5_lower']) and HU5==F(cert['M_power6_upper'])
 assert not cert['actual_M_lower_bound'] and not cert['joint_family_optimum_excluded'] and cert['prefix_states']==735471 and cert['both_endpoint_probability_laws_have_exact_unit_mass'] and not cert['prefix_renormalized']
 D=int(cand['denominator']);nums=cand['numerators'];assert D==2**44 and len(nums)==49 and nums[0]>0 and min(nums)>=0 and sum(nums)==D
 x=F(nums[32],D);assert x==F(2514712686865,2**44) and 0<x<1
 # Zero-based32 = 1+6*5+1: h1/128, order2; unnormalized coefficient a=x/beta.
 assert 32==1+6*5+1;coeff=[F(nums[0],D)/F(up['beta0'])]+[F(nums[1+6*j+k],D)/F(up['beta_by_h'][j][k])for j in range(8)for k in range(6)]
 caps=[F(up['beta0'])]+[F(v)for row in up['beta_by_h']for v in row]
 assert len(coeff)==49 and min(coeff)>=0 and coeff[0]>0 and sum(a*b for a,b in zip(coeff,caps))==1 and coeff[32]*beta==x
 Zlo=sum(a*(muL if i==32 else caps[i])for i,a in enumerate(coeff));assert Zlo==1-x+x*muL/beta>0
 theta=F(cert['threshold_M_power6_exact']);assert theta==F(up['threshold_M_power6_exact'])==F(8,9)**6*24576**5/65000**6
 bound=Zlo*HL5;assert bound>theta and bound/theta>F('1.0007863')
 ray=ob['fixed128_coefficient_ray_obstacle'];assert F(ray['weight_x_exact'])==x and F(ray['true_beta2_known_lower_exact'])==muL and F(ray['best_optimistic_fixed_ray_objective_power6_lower_exact'])==bound and ray['strict_full_fixed_ray_g0_obstacle'] and not ray['all_reoptimized_weights_excluded']
 reqL=beta*(1-(1-theta/HL5)/x);reqU=beta*(1-(1-theta/HU5)/x);assert muL>reqL>=reqU>0
 assert not ob['actual_M6_lower_bound'] and not ob['global_g0_proved'] and not ob['new_C6_bound'] and ob['new_functional_kernel_erfc_exp_log_root_optimizer_calls']==0
 for p,h in pins.items():assert sha(p)==h
 out={'verdict':'PASS_C143_CORRECTED_SAVED_OBSTACLES','source_obstacle_certificate_sha256':sha(P/'obstacle_certificate.json'),'registration_sha256':sha(A/'finite_registration.json'),'true_mu2_lower_reconstructed_exact':str(muL),'prefix_H_lower_identity_verified':True,'prefix_HL_power5_exact':str(HL5),'fixed_weight_x_exact':str(x),'Z_at_true_mu2_lower_exact':str(Zlo),'normalized_H_power5_lower_exact':str(bound),'threshold_theta_exact':str(theta),'lower_over_theta_exact':str(bound/theta),'lower_over_theta_display_only':display(bound/theta),'width_recovery_fraction_exact':str(ratio),'width_recovery_fraction_display_only':display(ratio),'old_leaf_count':5376,'top2188_width_strictly_below_required_recovery':True,'top2189_width_strictly_above_required_recovery':True,'minimum_necessary_old_leaf_replacements':2189,'fixed_partition_and_cost_scheme_only':True,'fixed_coefficient_ray_only':True,'actual_M_lower_bound':False,'reoptimized_family_excluded':False,'new_C6_bound':False,'new_functional_calls':0,'new_Arb_log_exp_root_native_optimizer_calls':0,'old_functional_replays':0,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
 write('verification_result.json',out);write('exit_receipt.json',{'status':out['verdict'],'exit_code':0,'result_sha256':sha(A/'verification_result.json'),'new_functional_calls':0,'new_Arb_log_exp_root_native_optimizer_calls':0,'old_functional_replays':0,'wall_seconds':out['wall_seconds'],'CPU_seconds':out['CPU_seconds'],'time_limit':None});print(json.dumps({k:out[k]for k in ['verdict','lower_over_theta_display_only','width_recovery_fraction_display_only','minimum_necessary_old_leaf_replacements','new_functional_calls','wall_seconds','CPU_seconds']}))
except BaseException:
 write('exit_receipt.json',{'status':'FAILED_STOP_NO_RETRY','failure':traceback.format_exc(),'new_functional_calls':0});raise
