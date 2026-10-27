"""Freeze exact actual173 outputs + code. AST/metadata/SHA only; no math import."""
from pathlib import Path
import ast,hashlib,json,os,time
H=Path(__file__).resolve().parent;T=Path('/Users/zhengzhihao/Desktop/T');B=H.parent/'binding020';P=T/'Math-Research-cube5-upper/cube-dim-5-bounds/upper_bandcurvature173';O=T/'cube-finite-dim-scratch/upper_bandcurvature173/run_v3_001'
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(Path(p).read_text())
def put(p,x):
    with Path(p).open('x') as f:json.dump(x,f,sort_keys=True,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def main():
    began=time.monotonic();assert not (H/'plan.json').exists();out=T/'cube-finite-dim-scratch/audit_bandcurvature174/numeric020_run001';assert not out.exists()
    for p in H.glob('*.py'):ast.parse(p.read_text(),filename=str(p))
    assert sha(B/'qualified_donor_endpoint_provenance_v2.json')=='87fa22239a24b1c2c818348648c80899b6246e1a3830eaa7512884647495834c'
    assert sha(B/'qualified172_174.json')=='9d4bf090d18541d532aa1373da4a16919a7df084fb6c65b7893f353365399625'
    assert sha(O/'math_outputs_manifest.json')=='28568080f4779a79b134c45b9e707ceb9d9bb2a78b64f17cff1b0a587995ce29';m=read(O/'math_outputs_manifest.json')
    assert m['prepared_v3_manifest_sha256']=='cceaa29e847719bc90cbbc1ef7005f8cf2c62c1924772d942003740233b773eb' and (m['actual_products'],m['actual_corner_cases'],m['actual_corner_points'])==(22295,34170,546720)
    assert sha(P/'final_local_binding173.json')=='da5fee821f1c7658cd50232b22a342cb51419ef1190276587a1a4357cd49d257'
    inputs=dict(read(B/'plan.json')['source_locks']);inputs.update(m['files']);inputs[str(O/'math_outputs_manifest.json')]=sha(O/'math_outputs_manifest.json');inputs[str(P/'final_local_binding173.json')]=sha(P/'final_local_binding173.json')
    prepared=P/'reviewed_v3/prepared_manifest.json';pm=read(prepared);assert sha(prepared)==m['prepared_v3_manifest_sha256'];inputs[str(prepared)]=sha(prepared)
    inputs.update(pm['source_locks']);inputs.update(pm['prepared_external_files']);inputs.update({str(prepared.parent/p):h for p,h in pm['repository_files'].items()})
    inputs[str(B/'qualified_donor_endpoint_provenance_v2.json')]=sha(B/'qualified_donor_endpoint_provenance_v2.json');inputs[str(B/'qualified172_174.json')]=sha(B/'qualified172_174.json')
    for p,h in inputs.items():assert sha(p)==h,p
    protected={str(B/'manifest.json'):sha(B/'manifest.json'),str(H.parent/'manifest.json'):sha(H.parent/'manifest.json')}
    for base in (B,H.parent):
        for p,h in read(base/'manifest.json')['files'].items():protected[str(base/p)]=h
    for p,h in protected.items():assert sha(p)==h,p
    finitepath=P/'packet/finite_plan.json';finite=read(finitepath);assert (finite['sources'],finite['weighted_band_terms'],finite['fixed_pairs'],finite['receivers'])==(455,22295,35872,1120)
    # Join only already accepted projection hashes; do not rerun455 source qualification.
    sourcechecks={}
    for p in sorted((B/'evidence/source_audit').glob('worker*/source_checks.jsonl')):
        for line in p.read_text().splitlines():r=json.loads(line);assert r['canonical_id'] not in sourcechecks;sourcechecks[r['canonical_id']]=r
    assert len(sourcechecks)==455
    masses=[];corners=[];seen=set()
    for sid in range(15):
        k=(7*sid)%15;ms=finite['phase1_mass_shards'][k];cs=finite['phase2_corner_shards'][k];md=O/'mass'/f'shard{k:02d}';cd=O/'corners'/f'shard{k:02d}';gp=O/f'original_g_input_shard{k:02d}.json';g=read(gp)
        for cid in read(ms['path'])['source_ids']:
            assert cid not in seen;seen.add(cid);assert g[cid]['minimal_projection_sha256']==sourcechecks[cid]['original_jet_projection_sha256'] and g[cid]['raw_record_sha256']==sourcechecks[cid]['original_raw_jet_sha256']
        common={'audit_shard':sid,'producer_shard':k}
        masses.append({**common,'input_path':ms['path'],'g_input_path':str(gp),'event_path':str(md/'events.jsonl'),'certificate_path':str(md/'source_certificates.json'),'checkpoint_path':str(md/'checkpoint.json'),'summary_path':str(md/'summary.json'),'registration_path':str(md/'registration.json')})
        corners.append({**common,'input_path':cs['input_path'],'event_path':str(cd/'events.jsonl'),'best_path':str(cd/'receiver_best.json'),'checkpoint_path':str(cd/'checkpoint.json'),'summary_path':str(cd/'summary.json'),'registration_path':str(cd/'registration.json')})
    assert seen==set(sourcechecks)
    plan={'task':'audit174_actual173v3_fixed_finite_math','producer_output':str(O),'producer_output_manifest_sha256':sha(O/'math_outputs_manifest.json'),'producer_prepared_manifest_sha256':sha(prepared),'registry_path':finite['source_registry_path'],'finite_plan_path':str(finitepath),'common_log_path':finite['old_common_log_path'],'times':finite['time_endpoints'],'mass_shards':masses,'corner_shards':corners,'output_path':str(out),'workers':15,'pilot':'one real fresh source then automatic remainder;one real corner event then remainder','source_cases':455,'weighted_band_products':22295,'all_pair_events':35872,'actual_corner_cases':34170,'actual_corner_points':546720,'source_checkpoint_every':1,'corner_checkpoint_at_most':128,'success_ACK_every':1,'normal_wall_cap':None,'automatic_restart':False,'resume_requires_primary_prefix_SHA_authorization':True,'old_source455_qualification_replayed':False,'new_H_log_gradient_calls':0,'shared_H_truth_dependency':'original017/DD/Arb accepted080/094/102/113/121 and g151/155/157/180;167 only transport','input_locks':inputs,'protected_locks':protected,'code_locks':{p.name:sha(p) for p in sorted(H.iterdir()) if p.is_file()},'actual_math_authorized':False}
    put(H/'plan.json',plan);put(H/'prepare_receipt.json',{'pid':os.getpid(),'seconds':time.monotonic()-began,'inputs':len(inputs),'protected':len(protected),'code_locks':len(plan['code_locks']),'new_math_function_calls':0,'no_live_pool_created':True,'plan_sha256':sha(H/'plan.json'),'actual_results_bound':True,'accepted455_source_projection_hash_join_only':True})
    files={p.name:sha(p) for p in sorted(H.iterdir()) if p.is_file()};put(H/'prepared_manifest.json',{'task':plan['task'],'files':files,'new_real_math_calls':0,'requires_primary_actual_authorization':True,'plan_sha256':sha(H/'plan.json'),'input_locks':len(inputs),'protected_locks':len(protected),'code_files':sum(p.endswith('.py') for p in plan['code_locks'])})
    print(json.dumps({'plan_sha256':sha(H/'plan.json'),'prepared_manifest_sha256':sha(H/'prepared_manifest.json'),'input_locks':len(inputs),'protected_locks':len(protected),'code_locks':len(plan['code_locks']),'new_real_math_calls':0},sort_keys=True))
if __name__=='__main__':main()
