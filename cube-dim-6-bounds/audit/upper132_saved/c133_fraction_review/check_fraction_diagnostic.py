import json,hashlib,sys,itertools,time,os
from pathlib import Path
from fractions import Fraction as Q
sys.set_int_max_str_digits(0);sys.dont_write_bytecode=True
os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
P=Path(__file__).resolve().parent;A=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-numerics/cube-dim-6-bounds/upper/round133/fixed_corner_envelope_diagnostic')
rd=lambda p:json.loads(Path(p).read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
reg=rd(A/'registration.json');diag=rd(A/'diagnostic.json');receipt=rd(A/'exit_receipt.json');cfg=rd(reg['selected_config_path']);source=rd(reg['source_result_path']);checks=0
def ck(x,label):
 global checks
 if not x:raise AssertionError(label)
 checks+=1
ck(sha(A/'diagnose.py')==reg['program_sha256'],'registered author program hash')
for row in reg['inputs']:ck(sha(row['path'])==row['sha256'],'original C132 input hash')
ck(receipt['status']=='SOLE_FIXED_BATCH_COMPLETE'and receipt['exit_code']==0 and receipt['diagnostic_sha256']==sha(A/'diagnostic.json'),'author terminal diagnostic binding')
ck(diag['registration_sha256']==sha(A/'registration.json')and diag['source_result_sha256']==sha(reg['source_result_path']),'diagnostic source/registration binding')
ck(diag['boxes']==len(diag['results'])==len(cfg['selected_boxes'])==93 and diag['corners']==2976,'complete93x32')
ck(reg['targets']==['1/63500','1/65000'],'two registered exact targets')
bitslist=list(itertools.product((0,1),repeat=5));classes=['FIXED_CERTIFICATE_FAMILY_OBSTRUCTED_AT_CORNER','WHOLEBOX_ALREADY_PASSED','ALL_CORNERS_PAID_LOCAL_SUBDIVISION_CANDIDATE_ONLY','INTERVAL_UNRESOLVED'];reconstructed=[];feasible_counts={t:0 for t in reg['targets']};feasible_boxes={t:0 for t in reg['targets']}
for choice,out in zip(cfg['selected_boxes'],diag['results']):
 src=rd(Path(reg['box_result_directory'])/(choice['node_id']+'.json'));ck(out['node_id']==src['node_id']==choice['node_id']and out['box']==src['box']==choice['box'],'same registered local box')
 corners=[]
 for ni,bits in enumerate(bitslist):
  L=U=Q(0);bands=[]
  for sb,cb in zip(src['band_envelopes'],choice['bands']):
   ck(sb['window']==cb['window'],'same disjoint band')
   rows=sb['center_bounds'];pairs=[]
   for center in rows:
    corner=center['corners'][ni];ck(corner['corner_bits']==list(bits),'unchanged C132 corner order');pairs.append((Q(corner['enclosure_lower']),Q(corner['enclosure_upper'])))
   l=max((z[0]for z in pairs),default=Q(0));u=max((z[1]for z in pairs),default=Q(0));L+=l;U+=u
   bands.append({'window':sb['window'],'max_center_robust_lower':str(l),'max_center_robust_upper':str(u),'lower_max_keys':[c['semantic_key']for c,pair in zip(rows,pairs)if pair[0]==l],'upper_max_keys':[c['semantic_key']for c,pair in zip(rows,pairs)if pair[1]==u],'missing_band_zero':not rows})
  corner={'corner_bits':list(bits),'face_point':[choice['box'][j][bit]for j,bit in enumerate(bits)],'fixed_robust_family_lower':str(L),'fixed_robust_family_upper':str(U),'band_envelopes':bands};ck(corner==out['corners'][ni],'independent exact pointwise max/sum enclosure');corners.append(corner)
 ck(len(out['corners'])==32,'all32 corners')
 targets={}
 for t in reg['targets']:
  T=Q(t);bad=[c for c in corners if Q(c['fixed_robust_family_upper'])<T];alllower=all(Q(c['fixed_robust_family_lower'])>=T for c in corners);prior=Q(src['summed_disjoint_band_lower'])>=T;ck(not prior or alllower,'wholebox implies paid corners')
  name=classes[0]if bad else classes[1]if prior else classes[2]if alllower else classes[3]
  worst=min(corners,key=lambda c:(Q(c['fixed_robust_family_upper']),c['corner_bits']));lowest=min(corners,key=lambda c:(Q(c['fixed_robust_family_lower']),c['corner_bits']))
  target={'classification':name,'original_wholebox_pass':prior,'obstruction_corner_count':len(bad),'all32_corner_lowers_ge_target':alllower,'interval_unresolved_corner_count':sum(Q(c['fixed_robust_family_lower'])<T<=Q(c['fixed_robust_family_upper'])for c in corners),'worst_upper_corner':{'corner_bits':worst['corner_bits'],'face_point':worst['face_point'],'upper_exact':worst['fixed_robust_family_upper'],'upper_minus_target_exact':str(Q(worst['fixed_robust_family_upper'])-T),'target_minus_upper_exact':str(T-Q(worst['fixed_robust_family_upper']))},'lowest_lower_corner':{'corner_bits':lowest['corner_bits'],'lower_exact':lowest['fixed_robust_family_lower'],'lower_minus_target_exact':str(Q(lowest['fixed_robust_family_lower'])-T)}};ck(out['targets'][t]==target,'all exact threshold classes/residuals');targets[t]=target
  feasible=[]
  for c in bad:
   x=[Q(1)]+list(map(Q,c['face_point']))+[Q(0)]
   if all(a>=b for a,b in zip(x,x[1:])):feasible.append(c)
  feasible_counts[t]+=len(feasible);feasible_boxes[t]+=bool(feasible)
 reconstructed.append({'node_id':choice['node_id'],'targets':targets})
aggregates={}
for t in reg['targets']:
 counts={name:sum(x['targets'][t]['classification']==name for x in reconstructed)for name in classes};worst=min(reconstructed,key=lambda x:(Q(x['targets'][t]['worst_upper_corner']['upper_exact']),x['node_id']))
 aggregates[t]={'counts':counts,'obstruction_corners_total':sum(x['targets'][t]['obstruction_corner_count']for x in reconstructed),'local_subdivision_candidate_node_ids':[x['node_id']for x in reconstructed if x['targets'][t]['classification']==classes[2]],'global_worst_saved_family_upper_corner':{'node_id':worst['node_id'],**worst['targets'][t]['worst_upper_corner']}};ck(diag['targets'][t]==aggregates[t],'exact93-box aggregates')
ck(diag['native_calls']==diag['exp_calls']==diag['checker_calls']==diag['new_subdivisions']==0,'zero new native/exp/subdivision author calls')
ck(diag['global_g_certificate']is False and diag['new_C6_bound']is False and diag['true_g_counterexample']is False,'fixed family only')
result={'status':'PASS_C133_EXACT_SAVED_DIAGNOSTIC_AND_FIXED_FAMILY_PROOF_REVIEW','assertions':checks,'boxes':93,'corners':2976,'aggregates':{t:{'counts':aggregates[t]['counts'],'obstruction_corners_total':aggregates[t]['obstruction_corners_total'],'ordered_face_obstruction_corners':feasible_counts[t],'boxes_with_ordered_face_obstruction':feasible_boxes[t]}for t in reg['targets']},'new_exp_native_producer_old_baseline_calls':0,'scope':'same per-original-box fixed selected jets and independent endpoint rectangle formula, complete closed axis-box coverage retaining exhibited corner; not G upper or true floor obstruction','diagnostic_sha256':sha(A/'diagnostic.json'),'program_sha256':sha(A/'diagnose.py'),'registration_sha256':sha(A/'registration.json')}
(P/'review_results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
