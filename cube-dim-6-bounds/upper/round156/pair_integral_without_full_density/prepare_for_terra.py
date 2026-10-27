#!/usr/bin/env python3
"""Terra-only freeze for one finite bound-family diagnostic; Sol unrun."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,json,shutil,sys
sys.set_int_max_str_digits(0)
SOURCE=Path(__file__).resolve().parent
BASE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper')
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def new(p,n,x):
    with(p/n).open('x')as f:json.dump(x,f,indent=2);f.write('\n')
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);args=ap.parse_args();P=Path(args.root).resolve()
    assert P!=SOURCE and P.parents[3]!=SOURCE.parents[3]and(P.parents[3]/'.git').exists()
    assert P.name=='pair_integral_without_full_density'and P.parent.name=='round156'and P.parents[1].name=='upper'and P.parents[2].name=='cube-dim-6-bounds'
    assert not any((P/n).exists()for n in ['inputs.json','registration.json','started.json','result.json'])
    P.mkdir(parents=True,exist_ok=True)
    for n in ['PROOF.md','prepare_for_terra.py','diagnose.py','HANDOFF_TO_TERRA.md']:
        if(P/n).exists():assert sha(P/n)==sha(SOURCE/n)
        else:shutil.copy2(SOURCE/n,P/n)
    C=BASE/'round130/true_second_count_boxes';U=BASE/'round128/true_beta1_joint_moment';R=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-gaussian106/cube-dim-6-bounds/upper/round146/raw_candidate_cost_certificate/candidate_certificate.json')
    assert sha(C/'cone_table.json')=='86f81fd6c69d8c9a6b7194d34201732c528ac4cde320aefbb19570d854097040'
    assert sha(R)=='1f529e6e2f8d414ea1bf05c7e756d3294b9b0d6cf1f72aa2b6137b110c630b2a'
    cones=read(C/'cone_table.json')['cones'];spec=read(C/'inputs.json');old=read(C/'results.json');u=read(U/'upstream.json');r=read(R)
    assert spec['h']==u['h_values'][5]=='1/128'and len(cones)==21
    assert Q(7,10**8)<Q(r['true_beta2_upper_threshold_for_strict_g0']['sufficient_exact'])and Q(8,10**8)<Q(r['true_beta2_upper_threshold_for_strict_saved128_improvement']['sufficient_exact'])
    files=[C/n for n in ['cone_table.json','inputs.json','results.json','box_core.py','proof.md']]+[U/'upstream.json',R]
    new(P,'source_manifest.json',{'sources':[{'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size}for p in files],'no149_or151_reserved_density_read':True})
    new(P,'inputs.json',{'h':'1/128','p_interval':u['p_intervals'][5],'kappa_interval':old['pi_factor_interval'],'R0':spec['R0'],'cones':cones,'saved146_legal_beta2_upper':r['current_valid_beta2_upper'],'goal7':'7/100000000','goal8':'1/12500000','source_p_is_exact_saved128_interval':True,'anchors':['0','sqrt(h)'],'max_raw_pairs':10416,'max_caps_per_pair':35,'max_cauchy_sqrt_cache_misses':156240,'global_budget_width_gate':'1/'+'1'+'0'*35})
    new(P,'registration.json',{'round':'156','status':'TERRA_PREPARED_UNRUN','science_invocations':1,'Arb_bits':512,'raw_pairs':10416,'max_caps_per_pair':35,'new_I0_closed_evaluations':1,'new_positive_k_anchor_evaluations':1,'new_exp_calls':3,'new_erfc_calls':2,'new_erf_calls':0,'max_cauchy_sqrt_cache_misses':156240,'workers':1,'nice_minimum':10,'time_or_CPU_limit':None,'automatic_restarts':0,'density_generation_MC_native_old6D_calls':0,'diagnose_sha256':sha(P/'diagnose.py'),'proof_sha256':sha(P/'PROOF.md'),'inputs_sha256':sha(P/'inputs.json'),'source_manifest_sha256':sha(P/'source_manifest.json')})
    print(json.dumps({'status':'PREPARED_UNRUN','root':str(P),'new_Arb_or_bound_family_evaluations':0,'awaiting_dual_gate':True}))
if __name__=='__main__':main()
