#!/usr/bin/env python3
"""Build the one new optimistic model, using JSON/Fraction only."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,importlib.metadata,importlib.util,json,sys
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent
A=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round128/true_beta1_joint_moment')
B=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round130/true_second_count_boxes')
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(n,x):(P/n).write_text(json.dumps(x,indent=2)+'\n')
def main():
    old=read(A/'upstream.json');new=read(B/'results.json');candidate=read(A/'candidate.json')
    beta=old['beta_by_h'];bL=F(new['moment_interval'][0]);bU=F(beta[5][1])
    assert old['h_values'][5]==new['h']=='1/128'and 0<bL<bU==F(new['old_beta2_interval'][1])
    assert len(beta)==8 and all(len(row)==6 and all(F(v)>0 for v in row)for row in beta)
    phi=old['Phi_table'];assert len(phi)==33 and all(len(row)==7 and row[0]==1 and all(isinstance(v,int)and v>=0 for v in row)for row in phi)
    assert F(old['beta0'])>0 and all(phi[n+1][k]>=phi[n][k]for n in range(32)for k in range(7))
    assert len(candidate['numerators'])==49 and min(candidate['numerators'])>=0 and candidate['numerators'][0]>0 and sum(candidate['numerators'])==int(candidate['denominator'])==2**44
    for side in (0,1):
        ps=[F(row[side])for row in old['p_intervals']];assert 0<ps[0]and ps[-1]<1 and all(a<b for a,b in zip(ps,ps[1:]))
    assert all(F(a)<=F(b)for a,b in old['p_intervals'])
    copied=[row.copy()for row in beta];copied[5][1]=str(bL)
    inp={k:old[k]for k in ['model','h_values','beta0','p_intervals','Phi_table','threshold_M_power6_exact','prefix_Kmax','prefix_state_count','tail_K_values']}
    inp.update({'beta_by_h':copied,'warm_start_49':candidate,'optimistic_changed_entry':[5,1],'optimistic_beta2_lower':str(bL),'original_beta2_upper':str(bU),'qualification':'optimistic count model only; lower moment endpoint is not a valid moment upper cap'})
    write('optimistic_inputs.json',inp)
    sources=[A/n for n in ['upstream.json','candidate.json','candidate_certificate.json','optimize_core.py','joint_core.py']]+[B/n for n in ['results.json','inputs.json']]
    write('source_manifest.json',{'sources':[{'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size}for p in sources],'read_live129':False,'old130_functional_calls_permitted':0})
    write('registration.json',{'round':'144','status':'PREPARED_ONLY_NOT_LAUNCHED','science_invocations':1,'optimizer_invocations':1,'max_objective_gradient_calls':150,'max_optimizer_calls_before_final_candidate':149,'max_iterations':120,'prefix_states':735471,'tail_rows':16,'features':49,'precision':'float64 diagnostic only','candidate_denominator':2**44,'x0_numerical_floor':'1/16777216','workers':1,'nice_minimum':10,'time_or_CPU_limit':None,'automatic_restarts':0,'new_kernel_erfc_heatkernel_native_MC_calls':0,'source_beta2_direction':'smaller positive beta2 lowers objective for every nonnegative parameter vector; optimistic only','certified_global_exclusion':False,'new_M_upper':False,'new_C6_bound':False,'numpy_version':importlib.metadata.version('numpy'),'scipy_version':importlib.metadata.version('scipy'),'run_sha256':sha(P/'run_prepared.py'),'proof_sha256':sha(P/'PROOF.md'),'inputs_sha256':sha(P/'optimistic_inputs.json'),'source_manifest_sha256':sha(P/'source_manifest.json')})
    spec=importlib.util.spec_from_file_location('c144_prepared_runner',P/'run_prepared.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    write('release_requirements.json',mod.release_requirements())
    print(json.dumps({'status':'PREPARED_ONLY','whole_domain_nonnegative_features_checked':True,'beta2_lower_strictly_below_old_upper':True,'new_objective_calls':0,'new_kernel_erfc_calls':0}))
if __name__=='__main__':main()
