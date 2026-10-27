#!/usr/bin/env python3
"""C129 terminal structure analysis + finite C141 input selection. No new box certificate, exp, or native call."""
import os,sys,json,hashlib,collections,subprocess
from pathlib import Path
from fractions import Fraction as Q
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;CUBE=P.parents[2];S129=CUBE/'upper/round129/recovered_continuation';S137=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-upper098/cube-dim-6-bounds/upper/round137/extra_time_band_corner_seed')
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(n,x):(P/n).write_text(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n')
def main():
 os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
 terminal=read(S129/'final_runtime_receipt_v129.json');counter=read(S129/'controller_counters_v129.json');detached=read(S129/'detached_exit_receipt.json');rr=read(S129/'run_receipt_v129.json');index=read(S129/'semantic_runtime_index_v129.json');forest=read(S129/'forest_v129.json');checkpoint=read(S129/'checkpoint_v129.json');initial=read(S129/'continuation_seed129.json')
 assert detached['exit_code']==0 and terminal['status']==forest['status']=='FINITE_FRONTIER_RUN_COMPLETE' and terminal['live_owned_processes']=={}
 assert counter['started_attempts']==counter['reaped_attempts']==18257 and counter['active_children']==0 and counter['reserved_attempts']==21384
 assert sha(S129/'semantic_runtime_index_v129.json')==terminal['semantic_runtime_index_sha256']==rr['semantic_runtime_index_sha256']
 assert sha(S129/'reservations_v129.jsonl')==terminal['reservation_journal_sha256']==index['reservation_journal_sha256']
 assert sha(S129/'controller_counters_v129.json')==index['controller_counters_sha256']
 assert len(forest['roots'])==325 and len(forest['continued_parent_roots'])==55
 # Validate all R137 declared original prepared bytes before consuming the five saved jets.
 r137reg=read(S137/'registration.json')
 for name,h in r137reg['prepared_pins'].items():assert sha(S137/name)==h
 assert sha(S137/'run_prepared.py')==r137reg['script_sha256']
 r137=read(S137/'results.json');jets=[];pinpaths=[S137/n for n in ['registration.json','inputs.json','results.json','run_receipt.json','source_manifest.json','parent_release.json','proof.md','run_prepared.py']]
 for row in r137['accepted_windows']:
  assert row['status']=='AVAILABLE'and row['accepted_precision_bits']==128
  for pathkey,hashkey in [('raw_path','raw_sha256'),('receipt_path','receipt_sha256'),('jet_path','jet_sha256')]:assert sha(row[pathkey])==row[hashkey];pinpaths.append(Path(row[pathkey]))
  ad=Path(row['attempt_directory']);nr=read(row['receipt_path']);assert nr['exit_code']==0 and nr['child_reaped']and nr['width_acceptable']
  for n,hkey in [('intent.json','intent_sha256'),('original_process.json','process_sha256'),('stderr.txt','stderr_sha256')]:assert sha(ad/n)==nr[hkey];pinpaths.append(ad/n)
  assert sha(row['raw_path'])==nr['stdout_sha256'];intent=read(ad/'intent.json');jet=read(row['jet_path']);assert jet['center']==r137['point']and jet['window']==row['window']and intent['semantic_key']==row['semantic_key'];jets.append({'semantic_key':row['semantic_key'],'window':row['window'],'jet_path':row['jet_path'],'jet_sha256':row['jet_sha256'],'raw_path':row['raw_path'],'raw_sha256':row['raw_sha256'],'receipt_path':row['receipt_path'],'receipt_sha256':row['receipt_sha256'],'intent_path':str(ad/'intent.json'),'intent_sha256':nr['intent_sha256']})
 assert len(jets)==5 and len(set(x['semantic_key']for x in jets))==5
 anchor=list(map(Q,r137['point']));g0=Q(1,65000);gs=Q(1,63500);rootrows=[];residual=[];leafcounts=collections.Counter();splitprefix=collections.Counter();frontierprefix=collections.Counter();terminalnodes=0;splitnodes=0;mechanism=collections.Counter();saved_crosschecks=0
 for rid,s in forest['roots'].items():
  nodes=s['nodes'];terminalnodes+=len(nodes);groups=[s[k]for k in ['passed_leaves','empty_leaves','queue','blocked_leaves']];flat=sum(groups,[]);assert len(flat)==len(set(flat));leaves={k for k,n in nodes.items()if n['status']!='split'};assert set(flat)==leaves
  assert set(s['passed_leaves'])=={c['node']for c in s['closed_coverage']}and len(s['passed_leaves'])==len(s['closed_coverage'])
  assert s['root_complete']==(not s['queue']and not s['blocked_leaves']);assert set(s['residual_frontier'])==set(s['queue']+s['blocked_leaves']);assert s['inflight']is None
  for node in nodes.values():
   if node['status']=='split':
    splitnodes+=1;assert len(node['children'])==2
    for cid in node['children']:assert cid in nodes and nodes[cid]['parent']==node['id']
    pre=node.get('partial_prefix')or{};plen=len(pre.get('artifact_keys',pre.get('keys',[])));splitprefix[str(plen)]+=1
    if pre.get('F_interval')and node.get('whole_box_interval'):
     flo,fhi=map(Q,pre['F_interval']);blo,bhi=map(Q,node['whole_box_interval']);saved_crosschecks+=1
     if flo>g0 and blo<g0:mechanism['center_above_g0_but_whole_box_failed_then_split']+=1;mechanism['this_split_prefix_length_'+str(plen)]+=1
     if flo<gs:mechanism['split_center_lower_not_yet_gstar']+=1
  for name,ids in zip(['passed','empty','queued','blocked'],groups):leafcounts[name]+=len(ids)
  if rid in forest['continued_parent_roots']:
   assert checkpoint['roots'][rid]['status']==s['status']and checkpoint['roots'][rid]['queue_count']==len(s['queue'])and checkpoint['roots'][rid]['passed_count']==len(s['passed_leaves'])and checkpoint['roots'][rid]['blocked_count']==len(s['blocked_leaves'])
   rootrows.append({'root_id':rid,'status':s['status'],'root_complete':s['root_complete'],'passed_leaves':len(s['passed_leaves']),'empty_leaves':len(s['empty_leaves']),'queue':len(s['queue']),'blocked':len(s['blocked_leaves']),'residual_frontier':len(s['residual_frontier']),'cumulative_C127_C129_native_used':s['native_reserved_129'],'cumulative_C127_C129_nodes_processed':s['nodes_processed_129'],'node_cap':1024,'native_cap':512})
  if not s['root_complete']:
   assert s['status']=='NATIVE_BUDGET_EXHAUSTED'and s['native_reserved_129']==512
   for nid in s['residual_frontier']:
    n=nodes[nid];pre=n.get('partial_prefix')or{};plen=len(pre.get('artifact_keys',pre.get('keys',[])));frontierprefix[str(plen)]+=1
    qmax=sum(max((Q(b[0])-c)**2,(Q(b[1])-c)**2)for b,c in zip(n['box'],anchor[1:]))/32
    residual.append({'root_id':rid,'node_id':nid,'box':n['box'],'status':n['status'],'stop_reason':n.get('stop_reason'),'saved_prefix_length':plen,'anchor_worst_frame32_quadratic':str(qmax),'source_complete_root':False})
 assert sum(s['root_complete']for s in initial['roots'].values())==276 and sum(initial['roots'][r]['root_complete']for r in forest['continued_parent_roots'])==6
 assert sum(s['root_complete']for s in forest['roots'].values())==290 and sum(r['root_complete']for r in rootrows)==20
 assert sum(r['cumulative_C127_C129_native_used']for r in rootrows)==21384
 residual.sort(key=lambda r:(Q(r['anchor_worst_frame32_quadratic']),r['root_id'],r['node_id']));selected=residual[:128]
 # Register the full finite candidate universe and every selected original leaf. Ranking is geometry only.
 write('terminal_structure_analysis.json',{'status':'TERMINAL_METADATA_AND_TREE_PARTITION_CHECKS_PASS_NOT_CERTIFICATE_REVERIFICATION','roots_total':325,'root_complete_claimed_by_forest':290,'unclosed_roots':35,'continued_roots':55,'continued_closed':20,'continued_preclosed_before129':6,'new_root_closures_in129':14,'whole_forest_closed_before129':276,'whole_forest_closed_after129':290,'newly_closed_root_ids':[r for r in forest['continued_parent_roots']if forest['roots'][r]['root_complete']and not initial['roots'][r]['root_complete']],'continued_unclosed':35,'root_status_counts':dict(collections.Counter(s['status']for s in forest['roots'].values())),'nodes_total':terminalnodes,'split_nodes':splitnodes,'terminal_leaf_counts':dict(leafcounts),'residual_frontier_total':len(residual),'residual_saved_prefix_length_counts':dict(frontierprefix),'split_saved_prefix_length_counts':dict(splitprefix),'split_mechanism_counts':dict(mechanism),'saved_split_prefix_box_interval_comparisons':saved_crosschecks,'inherited_passed_count_from_forest':forest['inherited_passed_count'],'new117_passed_leaf_list_count':len(forest['new117_passed_leaves']),'counter':counter,'native_129_started_reaped':18257,'C127_C129_reservations':21384,'unknown_C127_no_reissue_retained_in_index':sum(v['status']=='UNKNOWN_NO_REISSUE'for v in index['entries'].values()),'global_complete_coverage':False,'global_gstar_coverage':False,'global_g_certificate':False,'new_C6_bound':False,'mathematical_certificate_reverification':'Assigned to another Sol; this file checks terminal structure/accounting and saved interval logic only; no old kernel/certificate recomputation.'})
 write('continued55_root_summary.json',rootrows);write('all_residual_leaf_metadata.json',{'leaves':residual,'total':len(residual),'selection_universe_is_only_terminal_unresolved_original_leaves':True})
 write('selected_boxes_and_existing_jets.json',{'round':141,'status':'PREPARED_UNRUN','boxes':selected,'selected_box_count':len(selected),'source_residual_leaf_count':len(residual),'source_residual_root_count':35,'selection_rule':'at most128 terminal unresolved original leaves, ranked by exact max frame32 quadratic distance from fixed R137 boundary center; root ID/node ID tie-break; no certificate evaluated during selection','anchor':r137['point'],'jets':jets,'target_g0':'1/65000','target_gstar':'1/63500','odd_degrees':[1,3,5,7,9,11,13,15],'bands':r137reg['fixed_windows'],'max_corner_band_cases':len(selected)*32*5,'max_shared_Taylor_recurrence_degree':15,'max_Taylor_recurrence_updates':len(selected)*32*5*15,'new_native_calls_permitted':0,'new_exp_calls_permitted':0,'new_subdivisions_permitted':0,'time_or_CPU_ceiling':None,'workers':1,'nice_minimum':10,'g0_and_gstar_accounted_separately':True,'source_five_semantic_keys_no_reissue':[x['semantic_key']for x in jets],'failure_meaning':'algebraic lower below target is unresolved, not a true G/g counterexample or a robust-family impossibility witness'})
 terminalpaths=[S129/n for n in ['detached_exit_receipt.json','controller_counters_v129.json','run_receipt_v129.json','final_runtime_receipt_v129.json','forest_v129.json','checkpoint_v129.json','semantic_runtime_index_v129.json','run_prepared_129.py','registration_129.json','continuation_seed129.json']]
 inp=read(S137/'inputs.json');pinpaths.extend([Path(inp[k])for k in ['native_family_registration_path','native_family_validator_path','native_family_helpers_path','base_path']]);pins={str(q):sha(q)for q in terminalpaths+pinpaths};write('input_hashes.json',{'frozen_inputs':pins,'R137_prepared_pins_verified':True,'R137_all5_raw_receipt_jet_intent_process_stderr_hashes_verified':True,'C129_terminal_index_journal_counter_links_verified':True,'old_native_recomputations':0,'new_exp_calls':0,'C129_source_or_terminal_modified':False,'scope':'only terminal C129 metadata/tree and frozen R137 five saved source jets; no live129 or old kernel run'})
 firstlo=Q(r137['prefixes'][0]['F_interval'][0]);fivelo=Q(r137['prefixes'][-1]['F_interval'][0]);write('R137_reuse_context.json',{'five_semantic_keys':[x['semantic_key']for x in jets],'source_point':r137['point'],'first_band_center_lower_exact':str(firstlo),'first_band_lower_minus_gstar_exact':str(firstlo-gs),'five_band_center_lower_exact':str(fivelo),'point_value_is_not_a_whole_box_proof':True,'no_new_or_repeated_native_queries':True})
 print({'status':'PREPARED_UNRUN','roots_closed':290,'roots_residual':35,'continued_closed':20,'terminal_nodes':terminalnodes,'residual_leaves':len(residual),'residual_prefix_counts':dict(frontierprefix),'split_prefix_counts':dict(splitprefix),'mechanism':dict(mechanism),'selected_boxes':len(selected),'selected_roots':len(set(b['root_id']for b in selected)),'corner_band_budget':len(selected)*32*5,'new_exp':0,'new_native':0})
if __name__=='__main__':main()
