#!/usr/bin/env python3
"""Terra-only source freeze; Sol must not execute this preparation."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,sys
sys.set_int_max_str_digits(0)
P=Path(__file__).resolve().parent
M=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-localmax-audit/cube-dim-6-bounds/upper/round147/pair_pushforward')
C=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round130/true_second_count_boxes')
A=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-pushforward147-audit/cube-dim-6-bounds/audit')
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def new(n,x):
    with(P/n).open('x')as f:json.dump(x,f,indent=2);f.write('\n')
def main():
    assert not(P/'started.json').exists()and not(P/'certificate.json').exists()
    assert sha(M/'manifest.json')=='2500cf19538335300f9c4f2db74c914620acba279af8e182a1f43c96eb9de69f'
    assert sha(M/'prototype_result.json')=='5a492e403113d73aae9665735f0da68f8b501416e0c42548b0fae69e4b103f5b'
    assert sha(M/'next_result.json')=='1233ffbd6c9b093f3e99715affe3ef274e4d7fff374e23f9f3aa4c2c9889acd8'
    assert sha(A/'upper147_saved/manifest.json')=='2c64fc195f331064abb5c55b21165af30f3e68b26b8bf622546744d0da05ca6b'
    proto=read(M/'prototype_result.json');cones=read(C/'cone_table.json')['cones'];spec=read(C/'inputs.json')
    assert proto['cone_id']==0 and proto['body_pair_zero_based']==[16,17]and proto['nonnegative_density']is True
    cone=cones[0];b=cone['nodes_by_body'][16];c=cone['nodes_by_body'][17]
    assert b==['0','0','0','0','1/60','1/30']and c==['0','0','0','0','1/60','0']
    assert cone['alpha']==proto['alpha']=='1/171992678400'and spec['h']=='1/128'
    assert Q(1,60)*0-Q(1,30)*Q(1,60)==Q(-1,1800)and Q(1800,2*30**2)==1
    sources=[M/n for n in ['manifest.json','PROOF.md','prototype_result.json','next_result.json']]+[C/n for n in ['cone_table.json','inputs.json','box_core.py','proof.md']]+[A/'upper147_review/REVIEW.md',A/'upper147_saved/REVIEW.md',A/'upper147_saved/manifest.json']
    new('frozen_sources.json',{'sources':[{'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size}for p in sources],'h':'1/128','lambda':'30','cone_id':0,'body_pair_zero_based':[16,17],'alpha':'1/171992678400','density_coefficient':'1800','pair_factor':'900','C130_factor_numerator':'46080','raw_pair_multiplicity':1,'is_full_beta2':False,'source149_live_read':False})
    new('registration.json',{'round':'152','status':'TERRA_PREPARED_UNRUN','new_closed_transform_evaluations':1,'lambda':'30','h':'1/128','Arb_bits':512,'width_gate':'1/'+'1'+'0'*80,'new_erf_calls':1,'new_erfc_calls':1,'new_exp_calls':1,'new_sqrt_h_calls':1,'new_pi_constants':1,'pointwise_k_calls':0,'quadrature_calls':0,'old_343456_functional_replays':0,'new_density_generations':0,'MC_native_calls':0,'workers':1,'nice_minimum':10,'time_or_CPU_limit':None,'automatic_restarts':0,'case_sha256':sha(P/'single_case.py'),'proof_sha256':sha(P/'PROOF.md'),'general_interface_sha256':sha(P/'GENERAL_2D_INTERFACE.md'),'frozen_sources_sha256':sha(P/'frozen_sources.json')})
    print(json.dumps({'status':'PREPARED_UNRUN','real_cone_pair':[0,16,17],'new_numerical_evaluations':0,'scope':'one prototype contribution only'}))
if __name__=='__main__':main()
