#!/usr/bin/env python3
"""Prepare only: bind C141 saved rational data, semantic tables and finite diagnostic/anchor limits."""
import json,hashlib,os,sys
from pathlib import Path
from fractions import Fraction as Q
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;CUBE=P.parents[2];S141=CUBE/'upper/round141/boundary_anchor_rational_box_probe';S129=CUBE/'upper/round129/recovered_continuation';S127=CUBE/'upper/round127/continuation_v121';REC=CUBE/'audit/round127/crash_recovery';S137=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-upper098/cube-dim-6-bounds/upper/round137/extra_time_band_corner_seed')
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(n,x):(P/n).write_text(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n')
def canon(x):return(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
def semantic(center,w,model,definition,fan):
 return hashlib.sha256(canon({'center':[str(Q(x))for x in center],'a':str(Q(w['a'])),'T':str(Q(w['T'])),'formula_family_id':'C016-cauchy-frame32-outward-precision-family-v1','fan_sha256':fan,'model':model,'model_definition':definition,'quadratic':'frame32'})).hexdigest()
def main():
 os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
 ex=read(S141/'exit_receipt.json');assert ex['exit_code']==0 and sha(S141/'result.json')==ex['result_sha256']and sha(S141/'rational_point_cache.json')==ex['rational_cache_sha256']
 result=read(S141/'result.json');cache=read(S141/'rational_point_cache.json')['entries'];assert len(result['boxes'])==128 and len(cache)==5970 and result['gstar_local_box_pass_count']==result['g0_local_box_pass_count']==26 and result['unresolved_both_count']==102
 maxm=0;sum_m=0;mins={};maxs={}
 for k,row in cache.items():
  r=Q(row['exponent_exact']);a=Q(row['own_time_a']);f=Q(row['clamped_true_value']);assert f>0 and r==Q(row['gradient_dot_lower'])/f-Q(row['frame32_quadratic'])/a**2
  m=0;t=abs(r)
  while t>Q(1,2):t/=2;m+=1
  assert m<=16;maxm=max(maxm,m);sum_m+=m;w=row['window']['a'];mins[w]=min(mins.get(w,r),r);maxs[w]=max(maxs.get(w,r),r)
 cfg141=read(S141/'selected_boxes_and_existing_jets.json');r137=read(S137/'registration.json');i137=read(S137/'inputs.json')
 for meta in cfg141['jets']:assert semantic(cfg141['anchor'],meta['window'],i137['model'],i137['model_definition'],r137['fan_sha256'])==meta['semantic_key']
 prior=read(S127/'cache_crosswalk_12081.json');c127=read(REC/'recovered_attempt_index.json');c129=read(S129/'semantic_runtime_index_v129.json');statuses={k:'AVAILABLE'for k in prior['entries']}
 for index in [c127,c129]:
  for k,m in index['entries'].items():
   status=m['status'];assert status in ['AVAILABLE','UNKNOWN_NO_REISSUE']
   if k in statuses:assert statuses[k]==status
   statuses[k]=status
 for k in cfg141['source_five_semantic_keys_no_reissue']:statuses[k]='AVAILABLE_R137'
 universe=read(S141/'all_residual_leaf_metadata.json');assert universe['total']==len(universe['leaves'])==10266
 write('existing_semantic_statuses.json',{'entries':statuses,'frozen_only':True,'all_prior_and_unknown_queries_must_not_reissue':True,'future_main_must_merge_other_rounds_and_exclusive_owner_registry':True})
 cfg={'round':145,'status':'PREPARED_UNRUN','C141_result_path':str(S141/'result.json'),'C141_rational_cache_path':str(S141/'rational_point_cache.json'),'C129_residual_metadata_path':str(S141/'all_residual_leaf_metadata.json'),'source_boxes':128,'source_lower_unresolved_boxes':102,'source_existing_pass_boxes':26,'source_point_count':5970,'max_point_rational_enclosures':5970,'small_negative_argument_interval':['-1/2','0'],'initial_Taylor_lower_degree':15,'initial_Taylor_upper_degree':16,'dyadic_outward_fraction_bits':256,'max_range_reduction_squarings':16,'actual_saved_data_max_squarings':maxm,'actual_saved_data_sum_squarings':sum_m,'max_Taylor_updates':5970*16,'max_dyadic_squaring_steps':5970*16,'max_box_corner_cases':128*32,'targets':['1/63500','1/65000'],'source_windows':cfg141['bands'],'max_projected_boundary_candidate_pool':64,'max_shared_anchor_proposals':8,'max_geometry_box_pairs':64*10266,'geometry_proxy':'K_first=max_corner ||x-c||²/(32*(1/32)^2)<=1/4; equivalently K_second<=1; geometric proposal proxy only, not proof of new source effectiveness','candidate_source':'only ordered corners with reliable original fixed-anchor point-family upper<gstar, projected to c2=1 and optionally c3=1; no box-midpoint anchors','candidate_new_query_status':'proposal only; no native owner authorization, no execution path','max_future_semantic_keys_in_proposals':40,'future_native_attempt_cap_proposal_only':120,'future_parallel_lanes_suggestion_only':15,'future_main_exclusive_owner_required':True,'source_model':i137['model'],'source_model_definition':i137['model_definition'],'fan_sha256':r137['fan_sha256'],'max_new_native':0,'max_new_exp':0,'max_checker_calls':0,'max_new_subdivisions':0,'workers':1,'nice_minimum':10,'time_or_CPU_ceiling':None,'saved_exponent_ranges_by_lower_a':{a:{'min':str(mins[a]),'max':str(maxs[a])}for a in mins},'classification_safety':'true-T upper is upper of robust certificate surrogate only, never actual G; lower failure alone unresolved; ordered clipping distinction preserved'}
 write('inputs.json',cfg)
 paths=[S141/n for n in ['run_batch.py','PROOF.md','result.json','rational_point_cache.json','exit_receipt.json','selected_boxes_and_existing_jets.json','registration.json','all_residual_leaf_metadata.json']]+[S129/'semantic_runtime_index_v129.json',S127/'cache_crosswalk_12081.json',REC/'recovered_attempt_index.json',S137/'registration.json',S137/'inputs.json',S137/'results.json']
 audit141=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-upper098/cube-dim-6-bounds/audit/upper141_saved');paths.extend([audit141/'scientific_seal.json',audit141/'verified_local_passes.json'])
 write('input_hashes.json',{'frozen_inputs':{str(p):sha(p)for p in paths},'C141_exit_result_cache_hash_chain_verified':True,'source_129_141_137_not_modified':True,'old_raw_jet_or_native_recomputations':0,'new_rational_exp_enclosures_computed_in_preparation':0,'native_calls':0,'transcendental_exp_calls':0,'checker_calls':0})
 print({'status':'PREPARED_UNRUN','saved_source_points':5970,'existing_pass':26,'unresolved':102,'observed_max_squarings':maxm,'sum_saved_squarings':sum_m,'frozen_semantic_keys':len(statuses),'future_anchor_proposal_cap':8,'new_native':0,'new_exp':0})
if __name__=='__main__':main()
