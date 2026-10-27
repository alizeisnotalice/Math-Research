#!/usr/bin/env python3
"""Freeze455 exact identity/provenance inputs before endpoint-field audit."""
from pathlib import Path
import ast,gzip,hashlib,json,os
H=Path(__file__).resolve().parent
T=Path('/Users/zhengzhihao/Desktop/T');B=T/'Math-Research-cube-dim-5/cube-dim-5-bounds';P=T/'Math-Research-cube5-upper/cube-dim-5-bounds/upper_bandcurvature173'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    p=Path(p)
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return json.load(f)
def put(p,x):
    with Path(p).open('x') as f:json.dump(x,f,sort_keys=True,indent=2);f.write('\n')
def main():
    assert not (H/'plan.json').exists();mp=P/'prepared_manifest.json';assert sha(mp)=='b218200d7ae0d185509ecb39bc0f11abc52a12842219bd885bbc57dd4bf1b3dd';m=read(mp)
    locks={str(mp):sha(mp),**m['source_locks'],**m['prepared_external_files'],**{str(P/n):h for n,h in m['repository_files'].items()}}
    registry=read(P/'packet/source_registry.json.gz');assert len(registry)==455
    for cid,s in registry.items():
        md=s['source171_metadata']['new_mass_donor'];locks[md['packet_path']]=md['packet_sha256'];bd=s['donor_binding']
        if bd['kind']=='legacy092_compact':locks[bd['original_raw_path']]=bd['original_raw_sha256']
    names=['upper_global080/engine.py','upper_global080/engine_v2.py','upper_global080/execution_policy_v2.json','upper_global080/roots.json.gz','upper_global080/initial_cells.json.gz','upper_global080/geometry017.json.gz',
      'upper_adaptive095/jobs.json.gz','upper_adaptive095/children_certificates.json.gz','upper_adaptive103/jobs.json.gz','upper_adaptive103/children_certificates.json.gz',
      'upper_global080/report.md','audit_range094/report.md','audit_adaptive102/report.md','audit_adaptive102/results.json','audit_adaptive113/report.md','audit_adaptive113/results.json',
      'audit_path121/report.md','audit_reuse168/execution_report.md','audit_reuse168/execution_manifest.json','audit_reuse168/candidate_binding001.json','audit_reuse168/run001/result.json',
      'round019/final_upper_review.json','strict_g180/portable/evidence/report.md']
    for n in names:locks[str(B/n)]=sha(B/n)
    binding168=read(B/'audit_reuse168/candidate_binding001.json')
    for p,v in binding168['candidate_locks'].items():locks[p]=v['sha256'] if isinstance(v,dict) else v
    G=P.parent/'upper_global151'
    for n in ('build151.py','check151.py','metadata.json','manifest.json'):locks[str(G/n)]=sha(G/n)
    old=H.parent/'manifest.json';assert sha(old)=='05c00900b867bad231394285f74167687977f09a8e3cb777824fbd07266a1179';locks[str(old)]=sha(old)
    for p,h in locks.items():assert sha(p)==h,p
    for p in H.glob('*.py'):ast.parse(p.read_text())
    plan={'task':'174_round020_fixed455_source_identity_and_provenance_audit','producer_manifest_sha256':sha(mp),'producer_directory':str(P),'registry_path':str(P/'packet/source_registry.json.gz'),
      'finite_plan_path':str(P/'packet/finite_plan.json'),'sources':455,'source_shards':15,'source_counts':[31]*5+[30]*10,'workers':15,'source_locks':locks,
      'code_locks':{p.name:sha(p) for p in H.glob('*.py')},'endpoint_field_decoding_authorized':True,'new_q_bandproducts_corner_H_log_gradient_evaluations_authorized':False,
      'identity_checks':['same point/model/epsilon','all73timekeys','092 original full raw equality and080 evidence table','095103 registration/kernel/run/fingerprint/jobs/certificate','151 jet/raw identity and physical component mapping','original155 Delta/gamma/q byte identity','168 complete baseline and71guard records'],
      'future_math_scope_only':{'source_q':455,'weighted_bands':22295,'maximum_pairs':35872,'corners_per_pair':16},'wall_cap':None,'retry':'no automatic restart;successful prefix preserved;fault stops15processes'}
    put(H/'plan.json',plan);print(json.dumps({'pid':os.getpid(),'plan_sha256':sha(H/'plan.json'),'source_locks':len(locks),'code_files':len(plan['code_locks']),'real_q_band_corner_calls':0}))
if __name__=='__main__':main()
