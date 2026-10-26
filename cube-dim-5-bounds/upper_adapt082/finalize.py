#!/usr/bin/env python3
"""Freeze two completed adaptive batches; no kernels or further dispatch."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import gzip,json,hashlib,time,os
import kernel as K
H=Path(__file__).resolve().parent;S=K.S;OLD=H.parent/'upper_global080';N=1<<160

def main():
 t=time.monotonic();state=json.loads(gzip.decompress((S/'state.json.gz').read_bytes()));assert len(state['batches'])==2;receipts=[];verifications=[];jobs=[];certs=[]
 for b in range(2):
  tag=f'{b:04d}';r=json.loads((H/('batch'+tag+'_receipt.json')).read_text());v=json.loads((H/('batch'+tag+'_verification.json')).read_text());assert r['status']=='complete' and v['status']=='pass' and r['workers_joined'] and v['all12_workers_joined'];receipts.append(r);verifications.append(v);jobs+=json.loads(gzip.decompress((H/('batch'+tag+'_jobs.json.gz')).read_bytes()));certs+=json.loads(gzip.decompress((H/('batch'+tag+'_cells.json.gz')).read_bytes()))
 assert len({j['key'] for j in jobs})==len(jobs)==len(certs);mapping={};payloads={}
 for j in jobs:
  ev=json.loads(gzip.decompress(K.evidencepath(j['key']).read_bytes()));ev.pop('runtime');mapping[j['key']]=hashlib.sha256(K.serial(ev).encode()).hexdigest();payloads[j['key']]=ev
 K.gzatomic(H/'mathematical_hashes.json.gz',mapping);ordered=sorted(certs,key=lambda c:c['weight']*(c['bounds']['integral_upper_i']-c['bounds']['integral_lower_i']),reverse=True);samples={jobs[0]['key'],jobs[-1]['key']}|{c['key'] for c in ordered[:8]};K.gzatomic(H/'endpoint_samples.json.gz',{k:payloads[k] for k in sorted(samples)})
 expanded={k:n for k,n in state['nodes'].items() if n['children']};updates={k:{a:state['nodes'][k][a] for a in ('selfL','selfU','effL','effU','children')} for k in state['roots'] if state['nodes'][k]['children']};K.gzatomic(H/'numeric_hierarchy.json.gz',{'units_denominator':str(N),'initial_source080_manifest_sha256':K.sha(OLD/'manifest.json'),'expanded_nodes':expanded,'root_updates':updates,'effective_global':receipts[-1]['after'],'physical_frontier_file':'batch0001_frontier.json.gz'})
 active=[]
 for r in receipts+verifications:
  for pid in r['worker_pids']:
   try:os.kill(pid,0);active.append(pid)
   except ProcessLookupError:pass
 assert not active and not list((S/'partial').glob('*.json.gz'))
 after=receipts[-1]['after'];summary={'status':'two_finite_batches_complete_verified_stop','children_per_batch':[r['completed_children'] for r in receipts],'parent_counts':[64,256],'new_children':len(certs),'new_H_count':sum(r['new_H_count'] for r in receipts),'cached_H1over64_count':sum(r['reused_H1over64_count'] for r in receipts),'new_jet_count':0,'new_cone_count':3504*len(certs),'raw_childsum_worse_parent_count':sum(len(r['raw_childsum_worse_parents']) for r in receipts),'expanded_level_counts':dict(Counter(n['level'] for n in expanded.values())),'frontier_level_counts':dict(Counter(state['nodes'][k]['level'] for k in state['frontier'])),'physical_frontier_count':len(state['frontier']),'hierarchy_node_count':len(state['nodes']),'budget':after,'production_program_walls':[r['wall_seconds'] for r in receipts],'readback_program_walls':[v['seconds'] for v in verifications],'all_workers_gone':True,'partial_checkpoints':0,'positive_band_count':sum(v['positive_bands'] for v in verifications),'endpoint_sample_count':len(samples),'all73_bands_preserved_external':True,'no085_consumed':True,'source080_manifest_sha256':K.sha(OLD/'manifest.json'),'source079_manifest_sha256':K.sha(H.parent/'upper_integral079/manifest.json'),'state_sha256':K.sha(S/'state.json.gz'),'numeric_hierarchy_sha256':K.sha(H/'numeric_hierarchy.json.gz'),'mathematical_hashes_sha256':K.sha(H/'mathematical_hashes.json.gz'),'samples_sha256':K.sha(H/'endpoint_samples.json.gz'),'seconds':time.monotonic()-t,'pid':os.getpid(),'source_sha256':K.sha(__file__)};K.atomic(H/'summary.json',summary);print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
