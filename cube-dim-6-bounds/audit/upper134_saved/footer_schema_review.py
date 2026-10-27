"""Zero-exp correction of the final cache-counter schema assertion; no oracle rerun."""
import json,sys,os,hashlib,itertools
from pathlib import Path
from fractions import Fraction as Q
sys.set_int_max_str_digits(0);os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
P=Path(__file__).resolve().parent;rd=lambda p:json.loads(Path(p).read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
registration=rd(P/'registration.json');A=Path(registration['author_directory']);receipt=rd(P/'actual_tool_receipt.json');checkpoint=rd(P/'verification_checkpoint.json');stderr=(P/'checker_stderr.txt').read_text()
assert receipt['exit_code']==1 and receipt['child_reaped']and receipt['source_pins_unchanged_after']
assert sha(P/'check_saved_768.py')==registration['checker_sha256'] and sha(P/'checker_stderr.txt')==receipt['stderr_sha256']
assert 'line 104' in stderr and 'AssertionError: all consolidated saved counts/results' in stderr
assert checkpoint['phase']=='SAVED_BFS_PARTITION_AND_COUNT_REVIEW' and checkpoint['roots_complete']==7 and checkpoint['new768_exp']==21126
for row in rd(P/'source_pins.json')['files']:assert sha(row['path'])==row['sha256'] and Path(row['path']).stat().st_size==row['bytes']
cfg=rd(A/'fixed_roots_and_atlas.json');result=rd(A/'result.json');pc=rd(A/'point_cache.json');target=Q(cfg['target']);roots=[rd(A/'roots'/(r['original_node_id']+'.json'))for r in cfg['roots']]
summary=[{k:v for k,v in r.items()if k!='nodes'}for r in roots]
globalcounts={'nodes_evaluated':sum(r['counts']['nodes_evaluated']for r in roots),'nodes_created':sum(len(r['nodes'])for r in roots),'split_attempts':sum(r['counts']['split_attempts']for r in roots),'geometry_splits':sum(r['counts']['geometry_splits']for r in roots)}
familycounts={k:sum(r[k]for r in roots)for k in ['new_exp','cache_hits_initial','cache_hits_new']};familycounts['old_corner_exp_recomputed']=0
# Correct schema: point_cache stores family counters only; result adds global node counters.
assert result['original_box_results']==summary and result['actual_counts']=={**globalcounts,**familycounts} and pc['actual_counts']==familycounts
assert globalcounts=={'nodes_evaluated':139,'nodes_created':139,'split_attempts':66,'geometry_splits':66}
assert familycounts=={'new_exp':21126,'cache_hits_initial':28062,'cache_hits_new':126940,'old_corner_exp_recomputed':0}
assert pc['initial_entry_count']==16929 and pc['total_entries']==len(pc['entries'])==38055 and pc['old_corner_exp_recomputed']==0
assert sum(familycounts[k]for k in ['new_exp','cache_hits_initial','cache_hits_new'])==139*32*32+66*16*32 and familycounts['new_exp']<=225792
assert result['native_calls']==0 and result['global_g_certificate']is False and result['new_C6_bound']is False
counts={'ordered':0,'outside_ordered':0};ranges={tag:{k:None for k in ['min_lower','max_lower','min_upper','max_upper']}for tag in counts};unique=set();lo_pass=unresolved=bad=0;core=[];allleaves=paidleaves=depthleaves=0
for root in roots:
 nodes=root['nodes'];frontier=root['closed_partition_leaf_ids'];paid=sum(nodes[n]['status']=='LOCAL_LEAF_PASS'for n in frontier);allpass=paid==len(frontier)
 assert root['local_pass']==allpass and root['obstruction']is None
 assert all(nodes[n]['status']in ['LOCAL_LEAF_PASS','DEPTH6_UNRESOLVED']for n in frontier)
 margin=min(Q(nodes[n]['whole_box_certificate']['summed_complete_box_lower'])for n in frontier)-target
 assert not allpass or margin>0
 core.append({'original_node_id':root['original_node_id'],'local_pass':allpass,'leaves':len(frontier),'paid_leaves':paid,'unresolved_leaves':len(frontier)-paid,'minimum_leaf_lower_minus_target_exact':str(margin),'counts':root['counts'],'new_exp':root['new_exp']})
 allleaves+=len(frontier);paidleaves+=paid;depthleaves+=len(frontier)-paid
 for node in nodes.values():
  if 'split_attempt'not in node:continue
  mids=node['split_attempt']['midface_point_certificates'];assert len(mids)==16
  for m in mids:
   p=list(map(Q,m['point']));isordered=p[0]==1 and p[-1]>=0 and all(a>=b for a,b in zip(p,p[1:]));assert m['point_is_ordered']==isordered
   tag='ordered'if isordered else'outside_ordered';counts[tag]+=1;unique.add(tuple(m['point']));lo=Q(m['fixed_family_lower']);hi=Q(m['fixed_family_upper']);assert lo<=hi
   lo_pass+=lo>=target;unresolved+=lo<target<=hi;bad+=hi<target
   for name,value in [('min_lower',lo),('max_lower',lo),('min_upper',hi),('max_upper',hi)]:ranges[tag][name]=value if ranges[tag][name]is None else min(ranges[tag][name],value)if name.startswith('min')else max(ranges[tag][name],value)
assert sum(r['local_pass']for r in roots)==result['locally_passed_original_boxes']==3
assert allleaves==73==7+66 and sum(counts.values())==66*16==1056 and bad==0
out={'status':'PASS_ONE_C134_SAVED_MATH_REVIEW_WITH_ZERO_EXP_FOOTER_SCHEMA_CORRECTION','new_exp_768_verified':21126,'old16929_exp_reevaluated':0,'completed_math_prefix_assertions_before_footer':checkpoint['assertions'],'math_run_exit_code':1,'math_run_failure_scope':'only final assertion compared family-only point_cache counters to combined result counters; all new768 endpoints and all7 BFS/closed partitions already passed','zero_exp_footer_schema_pass':True,'actual_counts':{**globalcounts,**familycounts},'original_rectangle_passes':3,'total_partition_leaves':allleaves,'paid_leaves':paidleaves,'depth6_unresolved_leaves':depthleaves,'midface_point_requests_ordered':counts['ordered'],'midface_point_requests_outside_ordered':counts['outside_ordered'],'midface_fixed_family_enclosure_ranges':{tag:{name:str(v)for name,v in stats.items()}for tag,stats in ranges.items()},'midface_lower_ge_target':lo_pass,'midface_interval_unresolved':unresolved,'ordered_obstructions':0,'outside_ordered_full_rectangle_obstructions':0,'unique_midface_physical_points':len(unique),'roots':core,'native_producer_old_baseline_evaluator_calls':0,'global_g_or_C6_proved':False,'precision_bits':768,'fresh_exp_reruns_after_prefix':0}
(P/'review_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k not in ['roots','midface_fixed_family_enclosure_ranges']},indent=2))
