#!/usr/bin/env python3
"""Freeze one new raw candidate and true cost; zero H/root computation."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,importlib.util,json,sys
sys.set_int_max_str_digits(0);sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent
OWN=P.parents[1]/'round144/optimistic_joint_beta2'
BASE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper')
S128=BASE/'round128/true_beta1_joint_moment';S130=BASE/'round130/true_second_count_boxes'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(n,x):(P/n).write_text(json.dumps(x,indent=2)+'\n')
def main():
    diagnostic=read(OWN/'result.json');model=read(OWN/'optimistic_inputs.json');old=read(S128/'upstream.json');r130=read(S130/'results.json')
    cand=diagnostic['candidate'];x=[F(v,int(cand['denominator']))for v in cand['numerators']]
    assert len(x)==49 and sum(x)==1 and x[0]>0 and min(x)>=0
    bL=F(r130['moment_interval'][0]);b=min(F(old['beta_by_h'][5][1]),F(r130['moment_interval'][1]))
    assert F(model['beta_by_h'][5][1])==bL<b==F(old['beta_by_h'][5][1])
    raw=[F(model['beta0'])]+[F(v)for row in model['beta_by_h']for v in row]
    caps=[F(old['beta0'])]+[F(v)for row in old['beta_by_h']for v in row];caps[32]=b
    assert all(raw[i]==caps[i]for i in range(49)if i!=32)
    a=[v/q for v,q in zip(x,raw)];C=sum(v*q for v,q in zip(a,caps));assert C==1-x[32]+x[32]*b/bL
    sourcefiles=[OWN/n for n in ['result.json','optimistic_inputs.json','run_prepared.py','exit_receipt.json']]+[S128/n for n in ['joint_core.py','upstream.json','candidate_certificate.json']]+[S130/'results.json']
    write('source_manifest.json',{'sources':[{'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size}for p in sourcefiles],'joint_core_path':str(S128/'joint_core.py'),'raw_a_quantized':False,'old_functional_calls_permitted':0})
    write('inputs.json',{'round':'146','raw_model_input':model,'candidate':cand,'raw_coefficients_exact':list(map(str,a)),'valid_cap_vector_exact':list(map(str,caps)),'changed_index':32,'x_changed_exact':str(x[32]),'bL_definition_only':str(bL),'current_valid_beta2_upper':str(b),'current_cost_exact':str(C),'threshold_g0_power6_exact':old['threshold_M_power6_exact'],'saved128_certificate_path':str(S128/'candidate_certificate.json'),'expected_state_hash':diagnostic['state_bytes_sha256'],'not_a_cost_one_candidate':True})
    q,r=divmod(735471,15);jobs=[];start=0
    for j in range(15):
        end=start+q+(j<r);jobs.append({'job':j,'start':start,'end_exclusive':end,'states':end-start});start=end
    assert start==735471
    write('jobs.json',jobs)
    write('registration.json',{'round':'146','status':'PREPARED_UNRUN','science_invocations':1,'new_candidates':1,'workers':15,'prefix_states':735471,'tail_rows':16,'logical_root_requests':735487,'arithmetic_bits':192,'root_cache_per_process_cap':65536,'nice_minimum':10,'BLAS_threads':1,'time_or_CPU_limit':None,'automatic_restarts':0,'raw_a_quantized':False,'new_optimizer_kernel_erfc_native_calls':0,'run_sha256':sha(P/'driver.py'),'proof_sha256':sha(P/'PROOF.md'),'inputs_sha256':sha(P/'inputs.json'),'jobs_sha256':sha(P/'jobs.json'),'source_manifest_sha256':sha(P/'source_manifest.json')})
    spec=importlib.util.spec_from_file_location('c146_frozen_driver',P/'driver.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);write('release_requirements.json',m.release_requirements())
    print(json.dumps({'status':'PREPARED_UNRUN','candidate':1,'prefix_states':735471,'workers':15,'raw_a_quantized':False,'new_H_or_root_calls':0}))
if __name__=='__main__':main()
