"""Final data/AST/SHA gates only; no endpoint JSON, band/q or corner call."""
from pathlib import Path
import ast,gzip,hashlib,json,os,time
H=Path(__file__).resolve().parent;W=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_bandcurvature173')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
def read(p):return json.load(open(p))
def write(p,x):
 assert not p.exists();p.write_text(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n')
def main():
 began=time.monotonic();assert not (H/'prepared_manifest.json').exists() and not (W/'run001').exists();plan=read(W/'finite_plan.json');sources=read(W/'source_registry.json');locks=read(W/'source_locks.json');assert len(sources)==455 and plan['sources']==455 and plan['weighted_band_terms']==22295 and plan['fixed_pairs']==35872 and plan['receivers']==1120 and plan['max_concurrent_workers']==15;assert len(plan['time_endpoints'])==73 and plan['time_endpoints'][24]=='1/64' and plan['time_endpoints'][32]=='1/32' and plan['time_endpoints'][-1]=='1';massseen=set();pairseen=set();receiverseen=set();used=set();prior=0
 for spec in plan['phase1_mass_shards']:
  assert sha(spec['path'])==spec['sha256'];d=read(spec['path']);assert d['shard_id']==spec['shard_id'] and set(d['source_ids'])==set(d['sources']) and len(d['sources'])==spec['sources']
  for cid,s in d['sources'].items():assert cid not in massseen and s==sources[cid];massseen.add(cid)
 assert massseen==set(sources)
 for spec in plan['phase2_corner_shards']:
  assert sha(spec['input_path'])==spec['sha256'];d=read(spec['input_path']);assert len(d['pairs'])==spec['pairs'] and len(d['receivers'])==spec['receivers']
  for k,r in d['receivers'].items():i=int(k);assert i%15==spec['shard_id'] and i not in receiverseen;receiverseen.add(i)
  for pair in d['pairs']:
   n=pair['pair_index'];i=pair['receiver_index'];cid=pair['canonical_id'];assert n not in pairseen and str(i) in d['receivers'] and cid in sources;pairseen.add(n);used.add(cid);prior+=pair['prior_test_same_old155_q'] is not None;assert d['sources'][cid]==sources[cid]['source171_metadata'] and pair['node']==d['receivers'][str(i)]['node']
 assert pairseen==set(range(35872)) and receiverseen==set(range(1120)) and used==set(sources) and prior==71
 for p,h in locks.items():assert sha(p)==h
 for p in H.glob('*.py'):ast.parse(p.read_text(),filename=str(p))
 P=H/'packet';assert not P.exists();P.mkdir();parts=[]
 for name in ['source_registry.json']:
  src=W/name;dest=P/(name+'.gz');hh=hashlib.sha256();size=0
  with open(src,'rb') as f,open(dest,'xb') as raw:
   with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0,compresslevel=9) as z:
    while b:=f.read(8*1024*1024):hh.update(b);size+=len(b);z.write(b)
  expected=hh.hexdigest();hd=hashlib.sha256();nd=0
  with gzip.open(dest,'rb') as f:
   while b:=f.read(8*1024*1024):hd.update(b);nd+=len(b)
  assert expected==hd.hexdigest()==sha(src) and nd==size and dest.stat().st_size<100_000_000;parts.append({'path':str(dest.relative_to(H)),'sha256':sha(dest),'bytes':dest.stat().st_size,'uncompressed_sha256':expected,'uncompressed_bytes':size,'actual_decode_verified':True})
 for name in ['finite_plan.json','source_locks.json','registration.json','preparation_receipt.json']:(P/name).write_bytes((W/name).read_bytes())
 receipt={'pid':os.getpid(),'scope':'prepared only; no actual endpoint source or math execution','sources':455,'future_band_weights':22295,'future_pairs':35872,'receivers':1120,'mass_unique_source_and_corner_unique_receiver_gates':True,'prior_tested_pairs_guarded':71,'code_validation':'AST only; math functions not called','jet_compressed_bytes_SHA_checked':True,'g151_jet_interval_fields_decoded_or_accessed':False,'donor_packets_opened_or_H_fields_read':False,'old_coefficients_copied_verbatim':True,'new_band_q_log_corner_H_gradient_mean_split_calls':0,'no_live_math_pool_created':True,'gzip_bytes':sum(p['bytes'] for p in parts),'all_decode_SHA_checked':True,'wall_seconds':time.monotonic()-began,'preparation_receipt_read_boundary':'g151_jet_data_opened=false denotes no JSON/value field parsing; compressed source bytes were read for SHA'};write(W/'prepared_validation_receipt.json',receipt);(P/'prepared_validation_receipt.json').write_bytes((W/'prepared_validation_receipt.json').read_bytes());manifest={'schema':'173_band_curvature_prepared_v1','repository_files':{p.name:sha(p) for p in sorted(H.iterdir()) if p.is_file()},'source_locks':locks,'prepared_external_files':{str(p):sha(p) for p in sorted(W.iterdir()) if p.is_file()},'compact_source_metadata_parts':parts,'packet_other_files':{str(p.relative_to(H)):sha(p) for p in P.iterdir() if p.suffix!='.gz'},'sources':455,'weighted_band_terms':22295,'fixed_pairs':35872,'receivers':1120,'max_concurrent_workers':15,'dispatch_authorized':False,'new_band_q_log_corner_H_gradient_mean_split_calls':0,'qualification_pending':['168 source best acceptance','172 primary and174 review','existing donor H source qualification'],'primary_166_paper_candidate_unaffected':True};write(H/'prepared_manifest.json',manifest);print(json.dumps({'prepared_manifest_sha256':sha(H/'prepared_manifest.json'),'finite_plan_sha256':sha(W/'finite_plan.json'),'source_registry_sha256':sha(W/'source_registry.json'),'receipt':receipt},sort_keys=True),flush=True)
if __name__=='__main__':main()
