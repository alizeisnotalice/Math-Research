#!/usr/bin/env python3
"""Metadata/hash-only final148 index preparation; no node/tree/margin scan."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib,json

HERE=Path(__file__).resolve().parent
T=Path('/Users/zhengzhihao/Desktop/T');S=T/'cube-finite-dim-scratch'
M=T/'Math-Research-cube-dim-5/cube-dim-5-bounds'
A=T/'Math-Research-cube5-audit123-sol/cube-dim-5-bounds'
B=T/'Math-Research-cube5-auditcurvature129-sol/cube-dim-5-bounds'
U=T/'Math-Research-cube5-upper/cube-dim-5-bounds'

def read(path):return json.loads(Path(path).read_bytes())
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        while b:=f.read(8*1024*1024):h.update(b)
    return h.hexdigest()
def pairs(x):
    out=set()
    if isinstance(x,dict):
        if isinstance(x.get('path'),str) and isinstance(x.get('sha256'),str):out.add((x['path'],x['sha256']))
        for k,v in x.items():
            if isinstance(k,str) and k.startswith('/') and isinstance(v,str) and len(v)==64:out.add((k,v))
            out.update(pairs(v))
    elif isinstance(x,list):
        for v in x:out.update(pairs(v))
    return out
def find(x,want,pointer):
    if x==want:return pointer
    if isinstance(x,dict):
        for k,v in x.items():
            answer=find(v,want,pointer+[k])
            if answer is not None:return answer
    elif isinstance(x,list):
        for i,v in enumerate(x):
            answer=find(v,want,pointer+[i])
            if answer is not None:return answer

def main():
    assert not (HERE/'index.json').exists(),'Preserve the reviewed final index'
    P={
        'seed':S/'upper_continue123/seed.json','mapper':S/'upper_continue123/accepted_mapper.json',
        'state123':S/'upper_continue123/batch0002/state.json',
        'curvature128':S/'upper_curvature128/evaluation001/candidate_results.jsonl',
        'cross137':S/'upper_parallel137/run001/master_attempt0001_state.json',
        'numeric138a':S/'upper_parallel137/residual_preflight/future_batch0001/state.json',
        'numeric138':S/'upper_parallel137/residual_preflight/future_batch0002/state.json',
        'cross141':S/'upper_parallel137/next2350_preflight/run001/master_attempt0001_state.json',
        'numeric143':S/'upper_parallel137/next1679_H_preflight/future_batch0001/state.json',
        'cross145':S/'upper_parallel137/next234_preflight/run001/master_attempt0001_state.json',
        'numeric147':S/'upper_parallel137/next126_H_preflight/future_batch0001/state.json',
        'primary014':M/'round014/primary_components_acceptance.json',
        'base123lock':A/'audit_g123_sol/source_lock.json','base123manifest':A/'audit_g123_sol/manifest.json',
        'base130lock':B/'audit_curvature130/numeric_run_plan.json','base130manifest':B/'audit_curvature130/manifest.json',
        'proof075':M/'audit_frame075/report.md','proof128':M/'upper_curvature128/report.md',
        'proof029':M/'upper_alpha1_029/report.md','proof083':M/'audit_globalg083/report.md',
        'proof_primary129':M/'round014/primary_anchor_curvature_proof.md',
        'constructor105':U/'upper_g110105/run.py','constructor107':U/'upper_g110107/certificate.py',
        'constructor118':U/'upper_reuse118/run.py','constructor123':U/'upper_continue123/certificate.py',
        'constructor_reuse123':U/'upper_continue123/reuse.py','constructor_cross':U/'upper_crossanchor131/consume.py',
        'constructor_residual':U/'upper_parallel137/residual_certificate.py',
        'source_H126_run':U/'upper_parallel137/H126_run.py',
        'source_H126_summary':S/'upper_parallel137/next126_H_preflight/future_batch0001/summary.json',
    }
    D={'137':M/'audit_parallel137','138':M/'audit_residual138','141':B/'audit_next141',
       '143':A/'audit_next143','145':B/'audit_later145','147':A/'audit_final147'}
    lock_names={'137':'audit_plan_v2.json','138':'source_lock.json','141':'audit_plan.json',
                '143':'source_lock.json','145':'audit_plan_ready.json','147':'source_lock.json'}
    for stage,folder in D.items():
        P['lock'+stage]=folder/lock_names[stage]
        P['manifest'+stage]=folder/'manifest.json'
        P['final'+stage]=folder/'final_check.json'
        P['report'+stage]=folder/'report.md'
    def row(item):
        name,p=item;return name,{'path':str(p),'sha256':sha(p)}
    with ThreadPoolExecutor(max_workers=15) as pool:sources=dict(pool.map(row,P.items()))
    assert sources['manifest145']['sha256']=='4bbc7d67aab766dd9ddac947dec8346cc187e0a9d940cd82e3a0a9f874746fde'
    assert sources['manifest147']['sha256']=='141eaf7bf1883d5fabeaa22b76eb07ede520648de5cdab56504cbd98f98e6ebb'
    configs={
      'base':(['seed','mapper','state123','curvature128'],['base123lock','base130lock']),
      '137':(['cross137'],['lock137']),'138':(['numeric138a','numeric138'],['lock138']),
      '141':(['cross141'],['lock141']),'143':(['numeric143'],['lock143']),
      '145':(['cross145'],['lock145']),'147':(['numeric147'],['lock147'])}
    components={stage:{'production_sources':prod,'lock_documents':locks,'qualification_requirements':[]}
                for stage,(prod,locks) in configs.items()}
    def need(stage,source,field,value):
        pointer=field if isinstance(field,list) else [field]
        v=read(P[source])
        for key in pointer:v=v[key]
        assert v==value,(stage,source,pointer,v)
        components[stage]['qualification_requirements'].append({'source':source,'pointer':pointer,'equals':value})
    need('base','primary014','status','primary_accepts_123_values_and128_130_partial_curvature_component')
    need('base','primary014','passed_frontier',104482)
    need('base','base123manifest','status','both_original_finite_new_value_and_structure_batches_independently_qualified')
    need('base','base123manifest','source_lock_sha256',sources['base123lock']['sha256'])
    need('base','base130manifest','status','finite11330_new_curvature_component_independently_qualified_partial')
    expected={
      '137':{'status':'1153_pass_certificates_protocol_sources_fulltree_independently_qualified_partial','passed':105635,'residual':4274,'all144_input_locks_unchanged':True,'winning_integer_gzip_losslessly_decoded_rows':1153,'source_lock_plan_sha256':sources['lock137']['sha256']},
      '138':{'status':'all4274_fresh_values_Delta_gradient129q_corner_journal_and_tree_gates_passed','new_H':4274,'new_passed':3099,'new_splits':1175,'no_faults':True,'source_lock_sha256':sources['lock138']['sha256']},
      '141':{'status':'671_single_anchor_passes_all_sources_protocol_and_complete_tree_qualified_partial','new_passes':671,'actual_exit_code':0,'all4698_input_hashes_unchanged':True},
      '143':{'status':'all1679_new_values_Delta_gradient129q_corners_journal_tree_and_exit_gates_passed','new_H':1679,'new_passed':1562,'no_fault':True,'source_lock_sha256':sources['lock143']['sha256']},
      '145':{'status':'finite_next234_integer_source_journal_fulltree_and_exit_gates_passed','winning_certificates':108,'actual_audit_exitcode':0,'source_locks_actually_checked':6565,'plan_sha256':sources['lock145']['sha256']},
      '147':{'status':'last126_new_values_Delta_gradient129q_corners_journal_fulltree_component_passed','new_H':126,'new_passed':126,'new_geometry_splits':0,'no_fault':True,'source_lock_sha256':sources['lock147']['sha256']}}
    for stage,requirements in expected.items():
        for field,value in requirements.items():need(stage,'final'+stage,field,value)
    # Each verdict and its exact lock document must be a member of its frozen manifest.
    for stage in D:
        manifest=read(P['manifest'+stage])
        for source in ('lock'+stage,'final'+stage):
            pointer=find(manifest['files'],sources[source]['sha256'],['files']);assert pointer is not None,(stage,source)
            need(stage,'manifest'+stage,pointer,sources[source]['sha256'])
    for source,manifestname in [('base123lock','base123manifest'),('base130lock','base130manifest')]:
        pointer=find(read(P[manifestname])['files'],sources[source]['sha256'],['files']);assert pointer is not None
        need('base',manifestname,pointer,sources[source]['sha256'])
    lineage_receipts={}
    for stage,cfg in components.items():
        inherited=set()
        for lock in cfg['lock_documents']:inherited.update(pairs(read(P[lock])))
        for source in cfg['production_sources']:
            r=sources[source];assert (r['path'],r['sha256']) in inherited,(stage,source)
        lineage_receipts[stage]=cfg['production_sources']
    index={'task':'audit_global148_final_index','all_components_final':True,'sources':sources,'components':components,
      'scope':'metadata freeze only; final leaf/tree/min-margin integration not executed',
      'model':'n5 alpha1 body16 unnormalised inverse Fourier(2pi)^-5',
      'nominal_target':'11/100000','formal_Cn_acceptance_owner':'primary'}
    with (HERE/'index.json').open('x') as f:json.dump(index,f,ensure_ascii=False,sort_keys=True,indent=2);f.write('\n')
    with (HERE/'index_preparation_receipt.json').open('x') as f:
        json.dump({'status':'all_final_component_hashes_and_direct_lineage_locks_frozen_no_global_scan',
          'index_sha256':sha(HERE/'index.json'),'source_count':len(sources),'component_lineage_bindings':lineage_receipts,
          'global_leaf_or_tree_or_min_margin_cases_executed':0,'new_H_log_gradient_corner_calls':0},f,sort_keys=True,indent=2);f.write('\n')
    print(json.dumps({'index_sha256':sha(HERE/'index.json'),'sources':len(sources),'global_cases_executed':0},sort_keys=True))

if __name__=='__main__':main()
