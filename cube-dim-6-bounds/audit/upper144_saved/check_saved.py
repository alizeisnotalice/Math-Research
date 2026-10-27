#!/usr/bin/env python3
"""One fixed new C144 candidate diagnostic replay, no optimizer or old functional."""
from pathlib import Path
from fractions import Fraction as F
from math import comb,factorial
from itertools import combinations_with_replacement
import hashlib,importlib.util,json,os,sys,time,traceback
sys.dont_write_bytecode=True;sys.set_int_max_str_digits(0)
P=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-gaussian106/cube-dim-6-bounds/upper/round144/optimistic_joint_beta2')
G=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round144/primary_gate')
A=Path(__file__).resolve().parent
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(n,d):
 with(A/n).open('x')as f:json.dump(d,f,indent=2);f.write('\n')
files=[P/n for n in ['run_prepared.py','PROOF.md','registration.json','optimistic_inputs.json','source_manifest.json','parent_release.json','started.json','result.json','exit_receipt.json']]+[G/n for n in ['terminal_receipt.json','stdout.txt','stderr.txt']]
pins={str(p):sha(p)for p in files}
write('finite_registration.json',{'status':'REGISTERED_BEFORE_FIXED_NEW_CANDIDATE_FLOAT_REPLAY','payload_hashes':pins,'fixed_new_candidates':1,'objective_replays':1,'prefix_states':735471,'tail_rows':16,'max_reciprocal_fifth_power_values':735487,'probability_method':'direct whole multinomial factor/power products; no probability log/exp','gradient_tangent_replays':0,'optimizer_invocations':0,'new_Arb_log_exp_native_kernel_calls':0,'old_functional_replays':0,'time_limit':None,'workers':1,'nice_minimum':10,'reviewer_script_sha256':sha(__file__)})
os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)));start=time.monotonic();cpu=time.process_time()
try:
 assert sha(P/'run_prepared.py')=='df219f5b955837a7f42b59e8191619ac79b9017a31e2329bcc413279423190e3'
 spec=importlib.util.spec_from_file_location('c144_saved_identity',P/'run_prepared.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);inp,reg,source=m.preflight()
 rel=read(P/'parent_release.json')
 for k,v in m.release_requirements().items():assert rel[k]==v
 assert rel['main_approved'] and rel['independent_review_approved'] and rel['independent_review_sha256']==sha(A.parent/'upper144_preflight/REVIEW.md')
 r=read(P/'result.json');x=read(P/'exit_receipt.json');s=read(P/'started.json');terminal=read(G/'terminal_receipt.json')
 assert r['status']==x['status']=='COMPLETE_ONE_OPTIMISTIC_NUMERICAL_DIAGNOSTIC' and x['result_sha256']==sha(P/'result.json')
 assert r['run_sha256']==terminal['source_sha256']==sha(P/'run_prepared.py') and r['release_sha256']==s['release_sha256']==terminal['release_sha256']==sha(P/'parent_release.json')
 assert terminal['exit_code']==0 and terminal['restarts']==0 and terminal['time_limit']is None and s['pid']==terminal['pid'] and s['nice']>=10
 assert sha(G/'stdout.txt')==terminal['stdout_sha256'] and sha(G/'stderr.txt')==terminal['stderr_sha256'] and (G/'stderr.txt').read_bytes()==b''
 assert r['objective_gradient_calls']==x['objective_gradient_calls']==57<=150 and r['optimizer_invocations']==x['optimizer_invocations']==1 and r['restarts']==0
 assert r['outcome']['success'] and r['outcome']['status']==0 and r['outcome']['iterations']==55<=120
 assert not any(r[k]for k in ['certified_global_exclusion','new_M_upper','global_g0_proved','new_C6_bound']) and not x['new_M_upper'] and not x['new_C6_bound']
 assert r['new_kernel_erfc_heatkernel_native_MC_calls']==r['old_functional_replays']==0 and x['new_kernel_erfc_native_calls']==0
 assert r['prefix_states']==735471 and r['tail_rows']==16 and r['features']==49 and r['feature_matrix_bytes']==735487*49*8
 assert r['new_count_probability_float_exp_calls']==735471 and r['new_probability_logfactorial_calls']==33 and r['new_category_probability_float_log_calls']==9
 nums=r['candidate']['numerators'];D=int(r['candidate']['denominator']);assert D==2**44 and len(nums)==49 and all(type(v)is int and v>=0 for v in nums) and nums[0]>0 and sum(nums)==D
 assert r['optimistic_beta2_lower']==inp['optimistic_beta2_lower']==inp['beta_by_h'][5][1]
 # Exact normalization/cost identities. Only true upper caps may enter a future Holder cost.
 old=read(source['upstream.json']);true130=read(source['results.json']);bL=F(inp['optimistic_beta2_lower']);bOld=F(old['beta_by_h'][5][1]);b130U=F(true130['moment_interval'][1]);assert 0<bL<=bOld<=b130U
 weights=[F(n,D)for n in nums];caps=[F(inp['beta0'])]+[F(v)for row in inp['beta_by_h']for v in row];coeff=[w/b for w,b in zip(weights,caps)];assert sum(a*b for a,b in zip(coeff,caps))==1
 xbeta=weights[32];assert 32==1+6*5+1 and coeff[32]*bL==xbeta
 costOld=sum(a*(bOld if i==32 else caps[i])for i,a in enumerate(coeff));cost130=sum(a*(b130U if i==32 else caps[i])for i,a in enumerate(coeff))
 assert costOld==1-xbeta+xbeta*bOld/bL and cost130==1-xbeta+xbeta*b130U/bL and 1<costOld<=cost130
 # One independent floating evaluation for this new dyadic candidate, using direct probabilities.
 assert all(os.environ.get(k)=='1'for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'])
 import numpy as np
 assert np.__version__==reg['numpy_version']
 n=735471;states=np.empty((n,8),dtype=np.uint8);countK=[0]*17
 for i,state in enumerate(combinations_with_replacement(range(17),8)):states[i]=state;countK[state[-1]]+=1
 assert i+1==n and countK==[comb(k+7,7)for k in range(17)] and hashlib.sha256(states.tobytes()).hexdigest()==r['state_bytes_sha256']
 mids=[float((F(a)+F(b))/2)for a,b in inp['p_intervals']];q=[mids[0]]+[mids[j]-mids[j-1]for j in range(1,8)]+[1-mids[-1]];assert min(q)>0
 cats=np.empty((n,9),dtype=np.uint8);cats[:,0]=states[:,0];cats[:,1:8]=states[:,1:]-states[:,:-1];cats[:,8]=32-states[:,-1]
 qtable=np.array([[p**k/factorial(k)for k in range(33)]for p in q],dtype=float)
 prob=np.full(n,float(factorial(32)))
 for j in range(9):prob*=qtable[j,cats[:,j]]
 assert np.all(prob>0) and np.all(np.isfinite(prob))
 a0=coeff[0];table=[[sum((coeff[1+6*j+k-1]*int(inp['Phi_table'][state][k])for k in range(1,7)),F(0))for state in range(33)]for j in range(8)]
 wf=np.full(n,float(a0))
 for j in range(8):wf+=np.array(list(map(float,table[j])))[states[:,j]]
 assert np.all(wf>0) and np.all(np.isfinite(wf))
 hprefix=float(np.sum(prob*wf**(-.2)));massprefix=float(np.sum(prob));htail=0.;masstail=0.
 for k in range(17,33):
  p=comb(32,k)*mids[-1]**k*(1-mids[-1])**(32-k);wt=float(a0+table[7][k]);assert wt>0 and p>0;htail+=p*wt**(-.2);masstail+=p
 h=hprefix+htail;power5=h**5;mass=massprefix+masstail
 assert abs(mass-1)<1e-12 and abs(h-r['H_float_at_exact_candidate'])<=2e-12*abs(h) and abs(power5-r['H_power5_float'])<=1e-11*abs(power5)
 theta=float(F(inp['threshold_M_power6_exact']));assert (power5<theta)==r['H_power5_below_g0_float_diagnostic_only']
 for p,v in pins.items():assert sha(p)==v
 out={'verdict':'PASS_C144_SAVED_NUMERICAL_DIAGNOSTIC_ONLY','registration_sha256':sha(A/'finite_registration.json'),'saved_result_sha256':sha(P/'result.json'),'exact_candidate_simplex_verified':True,'source_actual_objective_gradient_calls':57,'source_optimizer_invocations':1,'source_terminal_exit_code':0,'source_terminal_wall_seconds':terminal['wall_seconds'],'source_terminal_CPU_seconds':terminal['CPU_seconds'],'independent_fixed_candidate_objective_replays':1,'independent_H_float':h,'independent_H_power5_float':power5,'independent_prefix_H_float':hprefix,'independent_tail_H_float':htail,'independent_mass_unrenormalized_float':mass,'objective_relative_difference':abs(h-r['H_float_at_exact_candidate'])/abs(h),'changed_beta_weight_exact':str(xbeta),'normalization_cost_formula':'C(b)=1-x+x*b/bL; b must be a certified true moment upper','true_cap_old128_upper_exact':str(bOld),'cost_using_valid_old128_upper_exact':str(costOld),'true_cap_C130_upper_exact':str(b130U),'cost_using_valid_C130_upper_exact':str(cost130),'cost_Old_times_float_H5_diagnostic_only':float(costOld)*power5,'new_optimizer_invocations':0,'new_Arb_log_exp_native_kernel_calls':0,'old_functional_replays':0,'new_M_upper':False,'certified_global_exclusion':False,'new_C6_bound':False,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
 write('verification_result.json',out);write('exit_receipt.json',{'status':out['verdict'],'exit_code':0,'result_sha256':sha(A/'verification_result.json'),'new_fixed_candidate_objective_replays':1,'new_optimizer_invocations':0,'new_Arb_log_exp_native_kernel_calls':0,'old_functional_replays':0,'wall_seconds':out['wall_seconds'],'CPU_seconds':out['CPU_seconds'],'time_limit':None});print(json.dumps({k:out[k]for k in ['verdict','independent_H_power5_float','objective_relative_difference','cost_Old_times_float_H5_diagnostic_only','new_optimizer_invocations','wall_seconds','CPU_seconds']}))
except BaseException:
 write('exit_receipt.json',{'status':'FAILED_STOP_NO_RETRY','failure':traceback.format_exc(),'new_optimizer_invocations':0,'new_native_calls':0});raise
