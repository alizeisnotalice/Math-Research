#!/usr/bin/env python3
"""C134 selection/cache preparation only; no scientific exponential or kernel call."""
import json,hashlib,sys,os,itertools
from pathlib import Path
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;CUBE=P.parents[2];S132=CUBE/'upper/round132/cross_center_band_envelopes';S133=CUBE/'upper/round133/fixed_corner_envelope_diagnostic'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def canon(x):return(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
def write(n,x):(P/n).write_bytes(canon(x))
def cachekey(k,point):return hashlib.sha256(canon({'semantic_key':k,'point':point})).hexdigest()
def main():
 os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
 c133=read(S133/'summary.json');ids=c133['targets']['1/63500']['local_subdivision_candidate_node_ids'];assert len(ids)==7 and len(set(ids))==7
 config=read(S132/'selected_boxes_and_centers.json');byid={b['node_id']:b for b in config['selected_boxes']};sources=read(S132/'selected_source_cache_entries.json')['entries'];roots=[];cache={};keys=set();inputs=[S133/'summary.json',S133/'exit_receipt.json',S133/'registration.json',S132/'selected_boxes_and_centers.json',S132/'selected_source_cache_entries.json',S132/'input_hashes.json',S132/'source_check_receipt.json',S132/'run_batch.py',S132/'PROOF.md'];old_entries=0
 for nid in sorted(ids):
  choice=byid[nid];q=S132/'box_results'/(nid+'.json');inputs.append(q);result=read(q);assert not result['passes_new_local_target']
  root={'original_node_id':nid,'box':choice['box'],'bands':[],'original_corner_certificate_path':str(q),'original_corner_certificate_sha256':sha(q)}
  for band,selected in zip(result['band_envelopes'],choice['bands']):
   ks=[x['semantic_key']for x in band['center_bounds']];assert ks==[x['semantic_key']for x in selected['candidate_keys']]
   root['bands'].append({'window':band['window'],'candidate_keys':ks});keys.update(ks)
   for center in band['center_bounds']:
    k=center['semantic_key'];assert sources[k]['expected']['center']==center['true_center'];assert len(center['corners'])==32
    for c in center['corners']:
     point=['1']+[choice['box'][i][b]for i,b in enumerate(c['corner_bits'])];ck=cachekey(k,point);value={'semantic_key':k,'point':point,'lower':c['enclosure_lower'],'upper':c['enclosure_upper'],'origin':'C132_ORIGINAL_CORNER_NO_RECOMPUTATION','source_certificate_path':str(q),'source_certificate_sha256':sha(q),'source_corner_bits':c['corner_bits'],'exact_formula_fields':{x:c[x]for x in ['displacement','gradient_dot_lower','clamped_true_value','frame32_quadratic','exponent']}}
     if ck in cache:assert cache[ck]['lower']==value['lower']and cache[ck]['upper']==value['upper']
     else:cache[ck]=value
     old_entries+=1
  assert [len(b['candidate_keys'])for b in root['bands']]==[16,16,0,0,0]
  roots.append(root)
 # Reuse every already saved C132 source-point value for this unchanged fixed source union,
 # including values originally evaluated on another of the93 boxes. This expands cache only.
 eligible_records=0
 for q in sorted((S132/'box_results').glob('*.json')):
  if q not in inputs:inputs.append(q)
  saved=read(q);qsha=sha(q)
  for band in saved['band_envelopes']:
   for center in band['center_bounds']:
    k=center['semantic_key']
    if k not in keys:continue
    for c in center['corners']:
     eligible_records+=1;point=['1']+[saved['box'][i][b]for i,b in enumerate(c['corner_bits'])];ck=cachekey(k,point);value={'semantic_key':k,'point':point,'lower':c['enclosure_lower'],'upper':c['enclosure_upper'],'origin':'C132_ORIGINAL_CORNER_NO_RECOMPUTATION','source_certificate_path':str(q),'source_certificate_sha256':qsha,'source_corner_bits':c['corner_bits'],'exact_formula_fields':{x:c[x]for x in ['displacement','gradient_dot_lower','clamped_true_value','frame32_quadratic','exponent']}}
     if ck in cache:assert cache[ck]['lower']==value['lower']and cache[ck]['upper']==value['upper']
     else:cache[ck]=value
 write('initial_point_cache.json',{'entries':cache,'unique_source_point_entries':len(cache),'C132_original_corner_records_referenced':old_entries,'all93_eligible_C132_corner_records_referenced':eligible_records,'C132_original_box_certificate_files_consulted':93,'source_keys_still_fixed_original_union':len(keys),'new_exp_calls':0})
 write('fixed_roots_and_atlas.json',{'target':'1/63500','roots':roots,'root_count':7,'source_C133_summary_sha256':sha(S133/'summary.json'),'max_depth':6,'max_leaves_per_original_box':64,'max_nodes_total':889,'max_split_attempts_total':441,'max_new_exp':225792,'midface_points_per_split':16,'active_source_jets_per_original_box':32,'split_rule':'entire closed rectangle; choose longest exact rational side, ties lowest free coordinate index; split at exact midpoint; no ordered clipping/tightening','traversal':'original node IDs sorted; breadth-first queue, lower child before upper child','early_stop':'on any new-midface-point fixed-family robust upper<target, record full point and exact ordering, stop this original box; no ordered-domain claim if point outside wedge','missing_bands_zero':True,'native_calls_permitted':0,'precision_bits':512,'workers':1,'nice_minimum':10,'time_limit':None})
 write('fixed_source_entries.json',{'entries':{k:sources[k]for k in sorted(keys)},'source_count':len(keys),'no_expanded_atlas':True})
 originalpins=read(S132/'input_hashes.json');sourcepins={}
 for k in sorted(keys):
  m=sources[k]
  for path in [m['native_path'],m['receipt_path'],m['jet_path']if'jet_path'in m else m['jet_ref']['path']]:
   expected=originalpins['selected_original_artifacts'][path];assert sha(path)==expected;sourcepins[path]=expected
 frozen={str(q):sha(q)for q in inputs}
 for path,h in originalpins['frozen_sources'].items():
  # Only original math/authentication dependencies, no live runtime reads.
  if Path(path).name in ['registration_127.json','precision_family_helpers.py','native_validator_family.py','cover_core.py','inputs.json']:assert sha(path)==h;frozen[path]=h
 write('input_hashes.json',{'frozen_inputs':frozen,'fixed_selected_source_artifacts':sourcepins,'native_calls':0,'exp_calls':0,'C132_original_corner_records_referenced':old_entries,'unique_initial_source_point_entries':len(cache),'live129_inputs':0,'C132_C133_files_modified':False})
 print({'status':'PREPARED_NO_EXP','boxes':7,'unique_source_jets':len(keys),'initial_source_point_entries':len(cache),'old_records_reused':old_entries,'new_exp_cap':225792,'nodes_cap':889})
if __name__=='__main__':main()
