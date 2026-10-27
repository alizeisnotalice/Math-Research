#!/usr/bin/env python3
"""Independent C144 prepared-packet consistency; zero optimizer/scientific calls."""
from pathlib import Path
from fractions import Fraction as F
import ast,hashlib,importlib.util,json,sys,time
sys.dont_write_bytecode=True;sys.set_int_max_str_digits(0)
P=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-gaussian106/cube-dim-6-bounds/upper/round144/optimistic_joint_beta2')
A=Path(__file__).resolve().parent
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
start=time.monotonic()
assert sha(P/'run_prepared.py')=='df219f5b955837a7f42b59e8191619ac79b9017a31e2329bcc413279423190e3'
spec=importlib.util.spec_from_file_location('c144_independent_saved_preflight',P/'run_prepared.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
inp,reg,source=m.preflight();assert read(P/'release_requirements.json')==m.release_requirements()
manifest=read(P/'source_manifest.json')
for v in manifest['sources']:assert sha(v['path'])==v['sha256'] and Path(v['path']).stat().st_size==v['bytes']
old=read(source['upstream.json']);true=read(source['results.json']);warm=read(source['candidate.json']);cert=read(source['candidate_certificate.json'])
assert len(inp['h_values'])==len(inp['beta_by_h'])==len(inp['p_intervals'])==8
assert len(inp['Phi_table'])==33 and all(len(row)==7 and row[0]==1 for row in inp['Phi_table'])
assert all(type(v)is int and v>=0 for row in inp['Phi_table']for v in row)
assert all(len(row)==6 and all(F(v)>0 for v in row)for row in inp['beta_by_h'])
assert F(inp['h_values'][5])==F(true['h'])==F(1,128) and true['status']=='COMPLETE_FIXED_FINITE_BOX_PILOT'
changes=[(j,k)for j in range(8)for k in range(6)if F(inp['beta_by_h'][j][k])!=F(old['beta_by_h'][j][k])];assert changes==[(5,1)]
assert F(inp['beta_by_h'][5][1])==F(inp['optimistic_beta2_lower'])==F(true['moment_interval'][0]) and 0<F(inp['optimistic_beta2_lower'])<F(inp['original_beta2_upper'])==F(old['beta_by_h'][5][1])==F(true['old_beta2_interval'][1])
assert inp['optimistic_changed_entry']==[5,1] and inp['model']==old['model']
for k in ['h_values','beta0','Phi_table','p_intervals','threshold_M_power6_exact','prefix_Kmax','prefix_state_count','tail_K_values']:assert inp[k]==old[k]
assert inp['prefix_Kmax']==16 and inp['prefix_state_count']==735471 and inp['tail_K_values']==list(range(17,33))
for side in (0,1):
 ps=[F(v[side])for v in inp['p_intervals']];assert 0<ps[0]<ps[-1]<1 and all(a<b for a,b in zip(ps,ps[1:]))
 assert min([ps[0]]+[ps[j]-ps[j-1]for j in range(1,8)]+[1-ps[-1]])>0
assert all(F(a)<=F(b)for a,b in inp['p_intervals'])
nums=warm['numerators'];D=int(warm['denominator']);assert inp['warm_start_49']==warm and D==2**44 and len(nums)==49 and nums[0]>0 and min(nums)>=0 and sum(nums)==D
assert F(nums[0],D)>=F(1,2**24)
assert reg['max_objective_gradient_calls']==150 and reg['max_optimizer_calls_before_final_candidate']==149 and reg['max_iterations']==120 and reg['optimizer_invocations']==reg['science_invocations']==1 and reg['automatic_restarts']==0 and reg['time_or_CPU_limit']is None
assert not reg['new_M_upper'] and not reg['new_C6_bound'] and not reg['certified_global_exclusion']
# AST count is structural only; do not enumerate states, build features, call log/exp, import numpy/scipy or optimize.
tree=ast.parse((P/'run_prepared.py').read_text());mincalls=[n for n in ast.walk(tree)if isinstance(n,ast.Call)and isinstance(n.func,ast.Name)and n.func.id=='minimize'];assert len(mincalls)==1
assert not any((P/n).exists()for n in ['parent_release.json','started.json','result.json','exit_receipt.json'])
report={'verdict':'PASS_C144_PRE_EXECUTION_STATIC_REVIEW','run_sha256':sha(P/'run_prepared.py'),'proof_sha256':sha(P/'PROOF.md'),'registration_sha256':sha(P/'registration.json'),'inputs_sha256':sha(P/'optimistic_inputs.json'),'source_manifest_sha256':sha(P/'source_manifest.json'),'source_files_verified':len(manifest['sources']),'sole_changed_beta_coordinate':[5,1],'all_parameter_feature_nonnegativity_verified':True,'positive_whole_CDF_laws_verified':True,'warm_start_exact_simplex_verified':True,'single_static_minimize_call_site':True,'optimizer_call_cap':149,'final_exact_candidate_diagnostic_reserved_calls':1,'total_objective_gradient_cap':150,'new_objective_gradient_calls':0,'new_optimizer_invocations':0,'new_log_exp_root_native_kernel_calls':0,'science_started':False,'scope':'optimistic numerical selection only; no strict M upper or global exclusion','wall_seconds':time.monotonic()-start}
with (A/'static_checks.json').open('x')as f:json.dump(report,f,indent=2);f.write('\n')
print(json.dumps(report,indent=2))
