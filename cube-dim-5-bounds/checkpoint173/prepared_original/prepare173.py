"""Metadata and existing coeff copies only. Does NOT open jets or H packets."""
from pathlib import Path
from collections import Counter
import ast,gzip,hashlib,json,os,time,resource
from util173 import sha,digest,read,create
H=Path(__file__).resolve().parent;D=Path('/Users/zhengzhihao/Desktop/T');MAIN=D/'Math-Research-cube-dim-5/cube-dim-5-bounds';S=D/'cube-finite-dim-scratch';W=S/'upper_bandcurvature173';B=H.parent;T=S/'upper_tailmass171';G=B/'upper_global151'
EXPECTED_MODELS={'upper_global080':'n5 alpha1 body16 unnormalized inverse Fourier (2pi)^-5','upper_adapt082':'n5 alpha1 body16 unnormalized Fourier(2pi)^-5','upper_adaptive095':'n5 alpha1 body16 unnormalized inverse Fourier(2pi)^-5','upper_adaptive103':'n5 alpha1 body16 unnormalized inverse Fourier(2pi)^-5'}
def main():
 began=time.monotonic();assert not (W/'registration.json').exists();m171=B/'upper_tailmass171/manifest.json';assert sha(m171)=='7d68ad04d62d4fcd4db0596c1eb1ddeb67f784bbbc241202a5b714f71df6bdc6';mm=read(m171);locks={str(m171):sha(m171)}
 for p,h in mm['prepared_external_files'].items():assert sha(p)==h
 for p in [MAIN/'primary/band_curvature172.md',G/'manifest.json',G/'metadata.json',B/'upper_retarget155/result_manifest.json',S/'upper_retarget155/run001/common_log_certificate.json',S/'upper_reuse158/run001/complete_domain_candidate.json',S/'upper_reuse158/run001/receiver_best_updates.json',B/'upper_retarget155/packet/leaf_margins-0001.jsonl.gz']:
  locks[str(p)]=sha(p)
 create(W/'registration.json',{'pid':os.getpid(),'code_sha256':sha(__file__),'phase':'173consumer preparation only','source171_manifest_sha256':sha(m171),'sources':455,'pairs':35872,'receivers':1120,'source_workers_future':15,'corner_workers_future':15,'max_concurrent_workers':15,'wall_deadline':None,'donor_H_files_opened':0,'original_jet_data_files_opened':0,'new_q_band_corner_H_log_gradient_calls':0})
 fixed={r['canonical_id']:r for r in map(json.loads,open(T/'fixed_sources.jsonl'))};assert len(fixed)==455;models={};times=None;table='8dd2e9895d7465ac599d43feb1b5099506f4e43d1248215eba0a77e11f4d13f9'
 for stage,expected in EXPECTED_MODELS.items():
  p=MAIN/stage/'registration.json';reg=read(p);locks[str(p)]=sha(p);assert reg['model']==expected and reg['table_sha256']==table;model={'stage':stage,'registration_path':str(p),'registration_sha256':sha(p),'model':expected,'bodytable_sha256':table,'model_equality_basis':'locked n5 alpha1 body16 unnormalised inverse Fourier(2pi)^-5;spelling variations explicitly bound by source registrations and172'}
  manifest=MAIN/stage/'manifest.json';locks[str(manifest)]=sha(manifest);model['source_manifest_path']=str(manifest);model['source_manifest_sha256']=sha(manifest)
  if stage in ('upper_adaptive095','upper_adaptive103'):
   for n,key in [('run.py','run_sha256'),('kernel.py','kernel_sha256')]:assert sha(MAIN/stage/n)==reg[key];locks[str(MAIN/stage/n)]=reg[key]
   model['expected_packet_fingerprint']=hashlib.sha256((reg['run_sha256']+reg['kernel_sha256']+sha(p)).encode()).hexdigest()
   if times is None:times=reg['time_endpoints']
   else:assert reg['time_endpoints']==times
  models[stage]=model
 assert len(times)==73 and times[24]=='1/64' and times[32]=='1/32' and times[-1]=='1';r92=MAIN/'upper_range092/nodes.json.gz';locks[str(r92)]=sha(r92);legacy={r['packet_path']:r for r in read(r92)};assert len(legacy)==56139;reg92=read(MAIN/'upper_range092/registration.json');assert reg92['time_endpoints']==times
 for p in [MAIN/'upper_range092/registration.json',MAIN/'upper_range092/manifest.json',MAIN/'upper_range092/prepare.py',MAIN/'upper_path112/manifest.json',MAIN/'upper_path112/nodes.json.gz']:locks[str(p)]=sha(p)
 meta151=read(G/'metadata.json');assert meta151['model']=='n5 alpha1 body16 unnormalised inverse Fourier(2pi)^-5' and meta151['bodytable_document_id']==table;m151=read(G/'manifest.json');coeff={}
 for part in m151['datasets']['coefficients']:
  p=G/part['path'];assert sha(p)==part['sha256'];locks[str(p)]=part['sha256']
  with gzip.open(p,'rt') as f:
   for r in map(json.loads,f):
    if r['canonical_id'] in fixed:coeff[r['canonical_id']]=r
 assert set(coeff)==set(fixed);oldq={};m155=read(B/'upper_retarget155/result_manifest.json')
 # q_updates are existing coefficient records, no H interval for these EPS64 records.
 for p in sorted((S/'upper_retarget155/run001').glob('shard*/q_updates.jsonl')):
  locks[str(p)]=sha(p)
  for index,event in enumerate(map(json.loads,open(p))):
   r=event['body'];cid=r['canonical_id']
   if cid in fixed:
    assert r['epsilon']=='1/64' and r['mass_qualification'] is None and r['canonical_id'] not in oldq;ref=fixed[cid]['original_Delta_gamma_and_oldqualified155q_reference'];assert digest(r)==ref['source155_body_sha256'] and event['hash']==ref['source155_event_hash'];oldq[cid]={'epsilon':r['epsilon'],'point':r['point'],'Delta110_lower':r['Delta110_lower'],'Delta115_lower':r['Delta115_lower'],'loggradient_directed':r['loggradient_directed'],'oldqualified155q':r['q_new'],'raw_record_sha256':r['raw_record_sha256'],'source155_q_rule':r['q_new_rule'],'source155_original_path':str(p),'source155_line_index':index,'source155_event_hash':event['hash'],'source155_body_sha256':digest(r)}
 assert set(oldq)==set(fixed);sources={};stages=Counter()
 for cid in sorted(fixed):
  r=fixed[cid];co=coeff[cid];c=co['core'];q=oldq[cid];assert digest(c)==cid and c['epsilon']==q['epsilon']=='1/64' and c['point']==q['point']==r['point'];assert c['Delta110_lower']==q['Delta110_lower'] and c['loggradient_directed']==q['loggradient_directed'];assert co['raw_record_sha256']==q['raw_record_sha256']==r['original_H64_U_reference']['raw_record_sha256'];donor=r['new_mass_donor'];path=donor['packet_path']
  if path in legacy:
   n=legacy[path];assert n['packet_sha256']==donor['packet_sha256'] and n['key'] in donor['node_keys'];rawpath=n['original_raw_path'];stage=next((st for st in ('upper_global080','upper_adapt082') if '/'+st+'/' in rawpath),None);assert stage is not None;binding={'kind':'legacy092_compact','packet_key':n['key'],'original_raw_path':rawpath,'original_raw_sha256':n['original_raw_sha256'],'model_binding':models[stage],'require_raw_endpoint_array_exact_equal':True}
  else:
   stage=next((st for st in ('upper_adaptive095','upper_adaptive103') if '/'+st+'/' in path),None);assert stage is not None;binding={'kind':'full_adaptive_evidence','packet_key':donor['node_keys'][0],'model_binding':models[stage]}
  stages[stage]+=1;sources[cid]={'canonical_id':cid,'source171_metadata':r,'old155':q,'original_g151_jet_ref':{'jet_id':co['jet_id'],'raw_record_sha256':co['raw_record_sha256'],'model':meta151['model'],'epsilon':'1/64','point':c['point'],'interval_location':'values[0]=[lower,upper];read only by authorized actual consumer'},'donor_binding':binding}
 create(W/'source_registry.json',sources);mass=[]
 for sid in range(15):
  ids=[cid for i,cid in enumerate(sorted(sources)) if i%15==sid];p=W/f'mass_shard{sid:02d}_input.json';create(p,{'shard_id':sid,'source_ids':ids,'sources':{cid:sources[cid] for cid in ids}});mass.append({'shard_id':sid,'path':str(p),'sha256':sha(p),'sources':len(ids)})
 corner=read(T/'shard_plan.json')['shards']
 for spec in corner:assert sha(spec['input_path'])==spec['sha256'];locks[spec['input_path']]=spec['sha256']
 for part in m151['datasets']['jets']:locks[str(G/part['path'])]=part['sha256'] # hash verification reads bytes, not H fields
 create(W/'source_locks.json',locks);plan={'sources':455,'bands_per_source':49,'weighted_band_terms':22295,'fixed_pairs':35872,'receivers':1120,'max_corner_points':573952,'workers_per_phase':15,'max_concurrent_workers':15,'phase1_mass_shards':mass,'phase2_corner_shards':corner,'source_registry_path':str(W/'source_registry.json'),'source_registry_sha256':sha(W/'source_registry.json'),'original_g151_jet_parts':[{'path':str(G/p['path']),'sha256':p['sha256']} for p in m151['datasets']['jets']],'g151_model':meta151['model'],'g151_model_metadata_sha256':sha(G/'metadata.json'),'time_endpoints':times,'band_first_time':'1/64','band_last_finite_time':'1','last_band_infinity_included':True,'band_denominator':'16*max(original_g64_lower,donor_H64_lower)','same_cutoff_gate':'max(gL,d64L)<=min(gU,d64U);nonempty exact overlap required;else fault no clip','split_branch':'153 L=donorH32lower, U=original_g64_upper;0<=L<U','q_rule':'172 min(oldqualified155q,Qband,153 uniform+realH32split);no new H/log/gradient/mean','old_tested_pair_guard_count':71,'source_donor_stage_counts':dict(stages),'old_common_log_path':str(S/'upper_retarget155/run001/common_log_certificate.json'),'old_common_log_sha256':sha(S/'upper_retarget155/run001/common_log_certificate.json'),'source155_margin_part':str(B/'upper_retarget155/packet/leaf_margins-0001.jsonl.gz'),'source155_margin_sha256':sha(B/'upper_retarget155/packet/leaf_margins-0001.jsonl.gz'),'base158_updates_path':str(S/'upper_reuse158/run001/receiver_best_updates.json'),'base158_updates_sha256':sha(S/'upper_reuse158/run001/receiver_best_updates.json'),'base158_candidate_path':str(S/'upper_reuse158/run001/complete_domain_candidate.json'),'base158_candidate_sha256':sha(S/'upper_reuse158/run001/complete_domain_candidate.json'),'qualification_pending':['168','172 primary and174 theory/interface','donor endpoint source/model/identity'],'dispatch_authorized':False,'wall_deadline':None};create(W/'finite_plan.json',plan)
 for p,h in locks.items():assert sha(p)==h
 create(W/'preparation_receipt.json',{'pid':os.getpid(),'sources':455,'bands_future':22295,'pairs_future':35872,'source_stages':dict(stages),'old_coefficients_copied_only':True,'g151_jet_data_opened':False,'donor_H_files_opened':0,'new_q_band_log_corner_H_gradient_calls':0,'wall_seconds':time.monotonic()-began,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss});print(json.dumps(read(W/'preparation_receipt.json')),flush=True)
if __name__=='__main__':main()
