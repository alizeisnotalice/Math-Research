#!/usr/bin/env python3
"""Terra-only explicit-root source copy/freeze; never launch science."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,json,shutil,sys
sys.set_int_max_str_digits(0)
SOURCE=Path(__file__).resolve().parent
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def new(root,n,x):
    with(root/n).open('x')as f:json.dump(x,f,indent=2);f.write('\n')
def main():
    ap=argparse.ArgumentParser(description='C152 prepare-only; no Arb, generator, preflight or science import')
    ap.add_argument('--root',required=True,help='Terra own-clone destination packet directory')
    ap.add_argument('--primary-root',required=True,help='read-only main cube-dim-6-bounds directory')
    ap.add_argument('--prototype-root',required=True,help='read-only terminal147 pair_pushforward directory')
    ap.add_argument('--audit-root',required=True,help='read-only terminal147 independent audit directory containing upper147_review and upper147_saved')
    args=ap.parse_args();root=Path(args.root).resolve();primary=Path(args.primary_root).resolve();prototype=Path(args.prototype_root).resolve();audit=Path(args.audit_root).resolve()
    assert root!=SOURCE and root.parents[3]!=SOURCE.parents[3], 'EXECUTE_ONLY_IN_TERRA_INDEPENDENT_CLONE'
    assert root.name=='laplace_kernel_first_case'and root.parent.name=='round152'and root.parents[1].name=='upper'and root.parents[2].name=='cube-dim-6-bounds'
    assert(root.parents[3]/'.git').exists(), 'EXISTING_TERRA_CLONE_REQUIRED'
    assert not any((root/n).exists()for n in ['started.json','certificate.json','frozen_sources.json','registration.json','prepare_only_receipt.json']), 'NO_EXISTING_TASK_OR_PREPARED_OUTPUT_OVERWRITE'
    manifest=read(SOURCE/'MANIFEST_STATIC.json');whitelist=read(SOURCE/'WHITELIST_STATIC.json')
    for row in manifest['files']:assert sha(SOURCE/row['path'])==row['sha256'], 'STATIC_PACKET_HASH:'+row['path']
    root.mkdir(parents=True,exist_ok=True)
    for name in whitelist['source_files']:
        target=root/name
        if target.exists():assert sha(target)==sha(SOURCE/name), 'DO_NOT_OVERWRITE_DIFFERENT_SOURCE:'+name
        else:shutil.copy2(SOURCE/name,target)
    pins=read(root/'SOURCE_PINS_STATIC.json')
    # Exact terminal paths are explicit command inputs; no149 live source.
    expected_prototype=Path(pins['sources'][0]['path']).parent
    expected_primary=Path(pins['sources'][4]['path']).parents[3]
    expected_audit=Path(pins['sources'][8]['path']).parents[1]
    assert prototype==expected_prototype and primary==expected_primary and audit==expected_audit, 'READONLY_ROOTS_MUST_MATCH_FROZEN_PROVENANCE'
    for row in pins['sources']:assert sha(row['path'])==row['sha256'], 'TERMINAL_SOURCE_HASH:'+row['path']
    C=primary/'upper/round130/true_second_count_boxes';M=prototype
    proto=read(M/'prototype_result.json');cones=read(C/'cone_table.json')['cones'];spec=read(C/'inputs.json');cone=cones[0]
    assert proto['cone_id']==0 and proto['body_pair_zero_based']==[16,17]and proto['nonnegative_density']is True
    assert cone['nodes_by_body'][16]==['0','0','0','0','1/60','1/30']and cone['nodes_by_body'][17]==['0','0','0','0','1/60','0']
    assert cone['alpha']==proto['alpha']=='1/171992678400'and spec['h']=='1/128'
    assert Q(1,60)*0-Q(1,30)*Q(1,60)==Q(-1,1800)and Q(1800,2*30**2)==1
    sources=[{'path':row['path'],'sha256':row['sha256'],'bytes':Path(row['path']).stat().st_size}for row in pins['sources']]
    sources.append({'path':str(root/'SCOPE_ADDENDUM_v2.md'),'sha256':sha(root/'SCOPE_ADDENDUM_v2.md'),'bytes':(root/'SCOPE_ADDENDUM_v2.md').stat().st_size})
    new(root,'frozen_sources.json',{'sources':sources,'h':'1/128','lambda':'30','cone_id':0,'body_pair_zero_based':[16,17],'alpha':'1/171992678400','density_coefficient':'1800','pair_factor':'900','C130_factor_numerator':'46080','raw_pair_multiplicity':1,'is_full_beta2':False,'source149_live_read':False,'contribution_scope':pins['contribution_scope'],'not_original_fixed_pair_full_R6_event':True,'scope_addendum_sha256':sha(root/'SCOPE_ADDENDUM_v2.md')})
    new(root,'registration.json',{'round':'152','status':'TERRA_PREPARED_UNRUN_V2_SCOPE_EXPLICIT','new_closed_transform_evaluations':1,'lambda':'30','h':'1/128','Arb_bits':512,'width_gate':'1/'+'1'+'0'*80,'new_erf_calls':1,'new_erfc_calls':1,'new_exp_calls':1,'new_sqrt_h_calls':1,'new_pi_constants':1,'pointwise_k_calls':0,'quadrature_calls':0,'old_343456_functional_replays':0,'new_density_generations':0,'MC_native_calls':0,'workers':1,'nice_minimum':10,'time_or_CPU_limit':None,'automatic_restarts':0,'case_sha256':sha(root/'single_case.py'),'proof_sha256':sha(root/'PROOF.md'),'general_interface_sha256':sha(root/'GENERAL_2D_INTERFACE.md'),'frozen_sources_sha256':sha(root/'frozen_sources.json'),'scope_addendum_sha256':sha(root/'SCOPE_ADDENDUM_v2.md'),'not_original_fixed_pair_full_R6_event':True})
    new(root,'prepare_only_receipt.json',{'status':'PREPARE_ONLY_COMPLETE_NO_SCIENCE_LAUNCH','root':str(root),'primary_root':str(primary),'prototype_root':str(prototype),'audit_root':str(audit),'source_manifest_sha256':sha(root/'MANIFEST_STATIC.json'),'preparer_sha256':sha(root/'prepare_for_terra_v2.py'),'registration_sha256':sha(root/'registration.json'),'frozen_sources_sha256':sha(root/'frozen_sources.json'),'Arb_or_kernel_evaluations':0,'generator_runs':0,'scientific_preflight_runs':0,'scientific_launches':0})
    print(json.dumps({'status':'PREPARE_ONLY_COMPLETE','root':str(root),'science_started':False,'awaiting_main_and_independent_exact_hash_gate':True}))
if __name__=='__main__':main()
