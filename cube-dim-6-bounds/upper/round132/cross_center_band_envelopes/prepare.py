#!/usr/bin/env python3
"""Preparation only: exact nearest-center selection and conditional Fraction check, no exp."""
import json,hashlib,math,sys,os,ast
from pathlib import Path
from fractions import Fraction as Q
sys.dont_write_bytecode=True;sys.set_int_max_str_digits(0)
P=Path(__file__).resolve().parent;CUBE=P.parents[2];ROOT=CUBE.parent
CACHE=CUBE/'audit/round127/crash_recovery/combined_cache.json'
SEED=CUBE/'upper/round129/recovered_continuation/continuation_seed129.json'
M6=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-upper098/cube-dim-6-bounds/upper/round128/true_beta1_joint_moment/candidate_certificate.json')
THEOREM=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-theory/cube-dim-6-bounds/upper/round057/positive_mixture_curvature/theorem.md')
REG=CUBE/'upper/round127/continuation_v121/registration_127.json'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(n,x):(P/n).write_text(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n')
def main():
 os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
 cache=read(CACHE);seed=read(SEED);reg=read(REG)
 roots=seed['roots'];rid=next(k for k in sorted(roots)if not roots[k]['root_complete']);s=roots[rid];ids=sorted(set(s['queue']+s['blocked_leaves']));chosen=ids[:128]
 windows=[{'kind':'B','a':a,'T':t}for a,t in [('1/32','1/2'),('1/64','1/32'),('1/128','1/64'),('1/250','1/128'),('1/1000','1/250')]]
 candidates=[[ (k,m,[Q(x)for x in m['expected']['center'][1:]])for k,m in cache['entries'].items()if m['expected']['window']==w]for w in windows]
 boxes=[{'node_id':nid,'box':s['nodes'][nid]['box'],'original_status':s['nodes'][nid]['status'],'center':[sum(map(Q,b))/2 for b in s['nodes'][nid]['box']]}for nid in chosen]
 denominators=[q.denominator for group in candidates for _,_,c in group for q in c]+[q.denominator for box in boxes for q in box['center']];scale=math.lcm(*denominators)
 lattice=[[(k,m,[int(q*scale)for q in c])for k,m,c in group]for group in candidates];selected={};total=0
 for box in boxes:
  c=[int(q*scale)for q in box.pop('center')];box['bands']=[]
  for wi,w in enumerate(windows):
   nearest=sorted(((sum((x-y)**2 for x,y in zip(center,c)),k,m)for k,m,center in lattice[wi]),key=lambda x:(x[0],x[1]))[:16]
   rows=[]
   for dist,k,m in nearest:
    selected[k]=m;rows.append({'semantic_key':k,'distance_squared':str(Q(dist,scale*scale)),'true_center':m['expected']['center']})
   box['bands'].append({'window':w,'candidate_keys':rows,'missing_band_lower':'0'if not nearest else None});total+=len(nearest)*32
 exact=Q(read(M6)['M_power6_upper']);g=Q(1,63500);u5=Q(9,8)**6*exact/(384**5*g**6);margin=64**5-u5
 conditional={'status':'EXACT_CONDITIONAL_CHECK_ONLY','source_path':str(M6),'source_sha256':sha(M6),'source_field':'M_power6_upper','M6_exact':str(exact),'target_g':str(g),'U_power5_formula':'(9/8)^6*M6/(384^5*g^6)','U_power5_exact':str(u5),'comparison_rhs_64_power5':str(64**5),'rhs_minus_U_power5_exact':str(margin),'strict_less_than_64_power5':u5<64**5,'global_g_floor_assumed_not_proved':True,'global_g_certificate':False,'new_C6_bound':False}
 write('conditional_M6_check.json',conditional)
 config={'schema':'C132-cross-center-band-envelopes-v1','root_id':rid,'source_unresolved_leaf_count':len(ids),'selected_boxes':boxes,'selected_box_count':len(boxes),'selection_policy':'first lex unclosed original root; unresolved node IDs sorted, first at most128; exact Euclidean distance on common integer lattice, semantic key tie-break','all_root_unresolved_leaves_selected':len(ids)<=128,'source_cache_count':len(cache['entries']),'windows':windows,'nearest_centers_per_window':16,'candidate_counts_by_window':[len(g)for g in candidates],'selected_unique_jet_count':len(selected),'precision_bits':512,'target':'1/63500','reference_target':'1/65000','max_corner_exp':128*5*16*32,'actual_planned_corner_exp':total,'time_limit':None,'nice_minimum':10,'workers':1,'native_calls_permitted':0,'reuse_source_centers_and_gradients_unchanged':True,'same_band_aggregation':'max complete box lower endpoints across true centers; never add same-band centers','cross_band_aggregation':'sum five disjoint positive band lower bounds; missing bands0','scope':'local conditions on original unresolved boxes only; no global g/C6 claim'}
 write('selected_boxes_and_centers.json',config);write('selected_source_cache_entries.json',{'entries':selected})
 sourcepins={str(q):sha(q)for q in [CACHE,SEED,M6,THEOREM,REG,CUBE/'audit/round127/crash_recovery/recovery_receipt.json',CUBE/'upper/round129/recovered_continuation/unknown_geometry_receipt.json',CUBE/'upper/round127/continuation_v121/precision_family_helpers.py',CUBE/'upper/round127/continuation_v121/native_validator_family.py']}
 # Only chosen source artifacts are read; no whole original baseline validation.
 jetpins={}
 for k,m in selected.items():
  paths=[m['native_path'],m['receipt_path'],m['jet_path']if'jet_path'in m else m['jet_ref']['path']]
  for q in paths:jetpins[q]=sha(q)
  assert sha(m['native_path'])==m['native_sha256']
  assert sha(m['receipt_path'])==m['receipt_sha256']
  if'jet_path'in m:assert sha(m['jet_path'])==m['jet_sha256']
 inputs=read(reg['inputs']['117_inputs_path']);sourcepins[inputs['base_path']]=sha(inputs['base_path']);sourcepins[reg['inputs']['117_inputs_path']]=sha(reg['inputs']['117_inputs_path'])
 write('input_hashes.json',{'frozen_sources':sourcepins,'selected_original_artifacts':jetpins,'selected_unique_jets':len(selected),'prior15193_baseline_reverified':False,'native_calls':0,'exp_calls':0,'disclosure':'One initial discovery search accidentally displayed live C129 checkpoint text. It was not used for box selection, source data, proof, validation, or conclusions. Every C132 input is explicitly frozen and hash-listed here.'})
 print({'status':'PREPARED_NO_EXP','root':rid,'boxes':len(boxes),'all_unresolved_selected':len(ids)<=128,'candidates_by_band':config['candidate_counts_by_window'],'unique_jets':len(selected),'corner_exp_planned':total,'M6_conditional_strict_under64':u5<64**5})
if __name__=='__main__':main()
