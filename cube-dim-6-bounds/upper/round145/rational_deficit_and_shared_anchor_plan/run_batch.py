#!/usr/bin/env python3
"""Gated finite exact-rational diagnosis and geometry-only shared-anchor proposals; no native/exp calls."""
import os,sys,json,hashlib,itertools,math,time,traceback,datetime
from pathlib import Path
from fractions import Fraction as Q
sys.dont_write_bytecode=True;sys.set_int_max_str_digits(0)
P=Path(__file__).resolve().parent

def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def canon(x):return(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def atomic(n,x):
 p=P/n;t=p.with_name(p.name+'.pending')
 with t.open('wb')as f:f.write(canon(x));f.flush();os.fsync(f.fileno())
 os.replace(t,p)
def write_new(n,x):
 with(P/n).open('xb')as f:f.write(canon(x));f.flush();os.fsync(f.fileno())
def ordered(point):
 x=list(map(Q,point));return len(x)==6 and x[0]==1 and all(x[i]>=x[i+1]for i in range(5))and x[-1]>=0

def release_requirements():
 return {'round':'145','authorized':True,'scope':'SAVED_C141_RATIONAL_TRUE_T_DIAGNOSIS_AND_NONEXECUTABLE_ANCHOR_PROPOSALS_ONLY','native_calls_permitted':0,'transcendental_exp_calls_permitted':0,'checker_calls_permitted':0,'subdivisions_permitted':0,'workers':1,'nice_minimum':10,'time_or_CPU_ceiling':None,'point_enclosure_cap':5970,'Taylor_updates_cap':95520,'squaring_steps_cap':95520,'dyadic_fraction_bits':256,'max_squaring_steps_per_point':16,'max_geometry_box_pairs':657024,'max_anchor_proposals':8,'run_batch_sha256':sha(P/'run_batch.py'),'inputs_sha256':sha(P/'inputs.json'),'input_hashes_sha256':sha(P/'input_hashes.json'),'semantic_statuses_sha256':sha(P/'existing_semantic_statuses.json'),'proof_sha256':sha(P/'PROOF.md'),'prepared_packet_sha256':sha(P/'prepared_packet.json')}
def release_gate():
 p=P/'parent_release.json';require(p.is_file(),'MAIN_AND_INDEPENDENT_REVIEW_RELEASE_REQUIRED');r=read(p)
 for k,v in release_requirements().items():require(r.get(k)==v,'RELEASE_MISMATCH:'+k)
 require(r.get('main_approved')is True and r.get('independent_review_approved')is True,'DUAL_GATE_REQUIRED');return p

def source_preflight():
 for p,h in read(P/'input_hashes.json')['frozen_inputs'].items():require(sha(p)==h,'FROZEN_INPUT:'+p)
 cfg=read(P/'inputs.json');r=read(cfg['C141_result_path']);cache=read(cfg['C141_rational_cache_path'])['entries'];universe=read(cfg['C129_residual_metadata_path'])['leaves'];require(len(r['boxes'])==128 and len(cache)==5970 and len(universe)==10266,'FINITE_SOURCE_COUNTS');require(r['gstar_local_box_pass_count']==r['g0_local_box_pass_count']==26 and r['unresolved_both_count']==102,'C141_ACTUAL_CLASSIFICATION');require(cfg['dyadic_outward_fraction_bits']==256 and cfg['max_range_reduction_squarings']==16,'FIXED_RATIONAL_PARAMETERS')
 for row in cache.values():require(Q(row['clamped_true_value'])>0 and Q(row['exponent_exact'])==Q(row['gradient_dot_lower'])/Q(row['clamped_true_value'])-Q(row['frame32_quadratic'])/Q(row['own_time_a'])**2,'SAVED_TRUE_T_FORMULA_IDENTITY')
 return cfg,r,cache,universe

def exp_interval_rational(r,counters):
 """P15(-|t|)<=exp(-|t|)<=P16(-|t|), reciprocal for positive r, then outward dyadic squaring."""
 t=abs(r);m=0
 while t>Q(1,2):t/=2;m+=1
 require(m<=16,'FIXED_RANGE_REDUCTION_CAP');u=-t;term=Q(1);partial=Q(1);lo=None
 for k in range(1,17):
  term=term*u/k;partial+=term;counters['Taylor_updates']+=1;require(counters['Taylor_updates']<=95520,'TAYLOR_UPDATE_CAP')
  if k==15:lo=partial
 hi=partial;require(Q(0)<lo<=hi,'POSITIVE_SMALL_ARGUMENT_ENCLOSURE')
 if r>0:lo,hi=1/hi,1/lo
 S=1<<256
 def down(q):return Q((q.numerator*S)//q.denominator,S)
 def up(q):return Q((q.numerator*S+q.denominator-1)//q.denominator,S)
 lo,hi=down(lo),up(hi)
 for _ in range(m):
  lo,hi=down(lo*lo),up(hi*hi);counters['dyadic_squaring_steps']+=1;require(counters['dyadic_squaring_steps']<=95520,'SQUARING_STEP_CAP')
 require(Q(0)<=lo<=hi,'OUTWARD_DYADIC_ORDER')
 return lo,hi,m

def semantic(center,w,cfg):
 return hashlib.sha256(canon({'center':[str(Q(x))for x in center],'a':str(Q(w['a'])),'T':str(Q(w['T'])),'formula_family_id':'C016-cauchy-frame32-outward-precision-family-v1','fan_sha256':cfg['fan_sha256'],'model':cfg['source_model'],'model_definition':cfg['source_model_definition'],'quadratic':'frame32'})).hexdigest()

def make_anchor_proposals(cfg,boxes,universe,statuses):
 # Only actual ordered point-family deficits justify proposals; a weak polynomial failure alone does not.
 paid_for_planning={b['node_id']for b in boxes if b['targets']['1/63500']['classification']=='RATIONAL_TRUE_T_BOX_PAID'}
 universe_lookup={b['node_id']:b for b in universe}
 for box in boxes:
  if box['node_id']in paid_for_planning:require(box['node_id']in universe_lookup and box['box']==universe_lookup[box['node_id']]['box'],'EXACT_LOCAL_PASS_IDENTITY_FOR_PLANNING_ONLY')
 pool={};current_anchor=('1','1','1','3/8','1/4','1/8')
 for box in boxes:
  bad=[c for c in box['corners']if c['ordered']and Q(c['fixed_point_family_upper_exact'])<Q(1,63500)]
  if not bad:continue
  witness=min(bad,key=lambda c:(Q(c['fixed_point_family_upper_exact']),tuple(map(Q,c['point']))));gap=Q(1,63500)-Q(witness['fixed_point_family_upper_exact']);base=witness['point']
  for point in [['1','1',*base[2:]],['1','1','1',*base[3:]]]:
   point=list(map(lambda x:str(Q(x)),point));key=tuple(point)
   if key==current_anchor:continue
   require(ordered(point),'PROJECTED_BOUNDARY_ANCHOR_ORDER')
   entry=pool.setdefault(key,{'point':point,'justifying_original_boxes':[],'known_fixed_anchor_point_deficit_weight':Q(0)})
   if box['node_id']not in entry['justifying_original_boxes']:entry['justifying_original_boxes'].append(box['node_id']);entry['known_fixed_anchor_point_deficit_weight']+=gap
 candidates=sorted(pool.values(),key=lambda c:(-c['known_fixed_anchor_point_deficit_weight'],tuple(map(Q,c['point']))))[:64];eligible=[];unknown_excluded=[]
 for c in candidates:
  queries=[]
  for w in cfg['source_windows']:
   k=semantic(c['point'],w,cfg);queries.append({'window':w,'semantic_key':k,'frozen_status':statuses.get(k,'NEW_IN_FROZEN_SNAPSHOT_ONLY'),'reissue_allowed':False if k in statuses else None})
  if any(q['frozen_status']=='UNKNOWN_NO_REISSUE'for q in queries):unknown_excluded.append({'point':c['point'],'reason':'PRIOR_UNKNOWN_NO_REISSUE','queries':queries});continue
  c['queries']=queries;c['new_semantics_in_frozen_snapshot']=sum(q['frozen_status']=='NEW_IN_FROZEN_SNAPSHOT_ONLY'for q in queries);eligible.append(c)
 denominators=[Q(v).denominator for b in universe for pair in b['box']for v in pair]+[Q(v).denominator for c in eligible for v in c['point'][1:]];D=math.lcm(*denominators);integerboxes=[[(int(Q(lo)*D),int(Q(hi)*D))for lo,hi in b['box']]for b in universe];pairs=0
 known_deficits={b['node_id']:max(Q(0),Q(1,63500)-Q(b['targets']['1/63500']['true_T_box_upper_exact']))for b in boxes if b['targets']['1/63500']['classification']=='FIXED_ANCHOR_WHOLEBOX_TRUE_T_DEFICIT'}
 for c in eligible:
  center=[int(Q(x)*D)for x in c['point'][1:]];cover=set()
  for i,box in enumerate(integerboxes):
   if universe[i]['node_id']in paid_for_planning:continue
   pairs+=1;require(pairs<=657024,'FINITE_GEOMETRY_PAIR_CAP');dist=sum(max((lo-x)**2,(hi-x)**2)for(lo,hi),x in zip(box,center))
   if 128*dist<=D*D:cover.add(i)
  c['_cover']=cover
 picked=[];covered=set();remaining=eligible.copy()
 while remaining and len(picked)<8:
  def ranking(c):
   new=c['_cover']-covered;weight=sum(known_deficits.get(universe[i]['node_id'],Q(0))for i in new);return(-len(new),-weight,c['new_semantics_in_frozen_snapshot'],tuple(map(Q,c['point'])))
  best=min(remaining,key=ranking);remaining.remove(best);new=best['_cover']-covered
  if not new:break
  covered|=best['_cover'];item={k:v for k,v in best.items()if k!='_cover'};item['known_fixed_anchor_point_deficit_weight']=str(item['known_fixed_anchor_point_deficit_weight']);item['geometry_proxy_newly_covered_boxes']=len(new);item['geometry_proxy_all_covered_boxes']=len(best['_cover']);item['geometry_proxy_root_count']=len({universe[i]['root_id']for i in best['_cover']});item['geometry_proxy_covered_node_ids']=[universe[i]['node_id']for i in sorted(best['_cover'])];item['scientific_query_authorized']=False;item['unknown_new_source_values_and_gradients']=True;picked.append(item)
 return {'status':'GEOMETRY_ONLY_SHARED_BOUNDARY_ANCHOR_PROPOSALS_NOT_AUTHORIZED_QUERIES','original_terminal_residual_universe':len(universe),'locally_paid_node_ids_excluded_for_planning_only':sorted(paid_for_planning),'original_terminal_forest_not_modified_or_claimed_merged':True,'candidate_pool_before_cap':len(pool),'candidate_pool_capped':len(candidates),'unknown_anchor_candidates_excluded':unknown_excluded,'actual_geometry_box_pairs':pairs,'proposed_anchor_count':len(picked),'proposals':picked,'new_semantics_in_frozen_snapshot':sum(x['new_semantics_in_frozen_snapshot']for x in picked),'max_future_native_attempts_proposal_only':3*sum(x['new_semantics_in_frozen_snapshot']for x in picked),'main_must_merge_current_global_query_registry_and_exclusive_owner_before_any_future_issue':True,'proxy_is_not_a_certificate_or_claim_of_source_effectiveness':True,'future_same_band_combine_rule':'max whole-box lower across sources; never sum centers of one window','future_cross_band_combine_rule':'sum nonoverlapping positive bands','native_calls':0,'exp_calls':0}

def run():
 rp=release_gate();require(not(P/'started.json').exists(),'SOLE_BATCH_NO_RESTART');os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)));cfg,source,cache,universe=source_preflight();statuses=read(P/'existing_semantic_statuses.json')['entries'];write_new('started.json',{'pid':os.getpid(),'release_sha256':sha(rp),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'native_calls':0,'exp_calls':0});counters={'point_enclosures':0,'Taylor_updates':0,'dyadic_squaring_steps':0};bounds={};rows=[];t0=time.monotonic()
 try:
  for key,row in cache.items():
   lo,hi,m=exp_interval_rational(Q(row['exponent_exact']),counters);f=Q(row['clamped_true_value']);lower=max(f*lo,Q(row['band_corner_lower_exact']));upper=f*hi;require(Q(0)<=lower<=upper,'OLD_POLYNOMIAL_LOWER_MUST_FIT_TRUE_T_ENCLOSURE');bounds[key]={'semantic_key':row['semantic_key'],'point':row['point'],'range_squarings':m,'exp_lower_exact':str(lo),'exp_upper_exact':str(hi),'true_T_lower_exact':str(lower),'true_T_upper_exact':str(upper),'upper_scope':'robust tangent surrogate only, not actual band/G upper'};counters['point_enclosures']+=1;require(counters['point_enclosures']<=5970,'POINT_CAP')
  for original in source['boxes']:
   cornerrows={tuple(c['corner_bits']):{'point':cache[c['cache_key']]['point'],'L':Q(0),'U':Q(0)}for c in original['per_band'][0]['corner_references']};boxlo=Q(0);boxhi=Q(0);bandrows=[]
   for band in original['per_band']:
    lows=[];highs=[]
    for ref in band['corner_references']:
     v=bounds[ref['cache_key']];bits=tuple(ref['corner_bits']);require(v['point']==cornerrows[bits]['point']and v['semantic_key']==band['semantic_key'],'SAME_PHYSICAL_CORNER_BAND');lo,hi=Q(v['true_T_lower_exact']),Q(v['true_T_upper_exact']);lows.append(lo);highs.append(hi);cornerrows[bits]['L']+=lo;cornerrows[bits]['U']+=hi
    l=min(lows);u=min(highs);require(l<=u,'MIN_ENCLOSURE');boxlo+=l;boxhi+=u;bandrows.append({'window':band['window'],'true_T_wholebox_lower_exact':str(l),'true_T_wholebox_upper_exact':str(u),'upper_scope':'scalar robust complete-box certificate value only'})
   corners=[{'corner_bits':list(bits),'point':v['point'],'ordered':ordered(v['point']),'fixed_point_family_lower_exact':str(v['L']),'fixed_point_family_upper_exact':str(v['U'])}for bits,v in sorted(cornerrows.items())];targets={}
   for t in cfg['targets']:
    T=Q(t);old=original['passes_gstar']if t=='1/63500'else original['passes_g0'];paid=boxlo>=T;deficit=boxhi<T;classification='RATIONAL_TRUE_T_BOX_PAID'if paid else'FIXED_ANCHOR_WHOLEBOX_TRUE_T_DEFICIT'if deficit else'RATIONAL_ENCLOSURE_UNRESOLVED';bad=[c for c in corners if Q(c['fixed_point_family_upper_exact'])<T];orderedbad=[c for c in bad if c['ordered']];require(not old or paid,'PRESERVE_OLD_CERTIFIED_PASS');targets[t]={'classification':classification,'old_odd15_box_pass':old,'newly_paid_after_stronger_rational_enclosure':paid and not old,'true_T_box_lower_exact':str(boxlo),'true_T_box_upper_exact':str(boxhi),'box_upper_minus_target_exact':str(boxhi-T),'fixed_point_family_bad_corners':len(bad),'ordered_fixed_point_family_bad_corners':len(orderedbad),'only_full_rectangle_bad_corner_scope':bool(bad)and not orderedbad,'all_point_corner_lowers_paid':all(Q(c['fixed_point_family_lower_exact'])>=T for c in corners),'true_g_counterexample':False,'actual_G_upper_bound':False}
   rows.append({'node_id':original['node_id'],'root_id':original['root_id'],'box':original['box'],'targets':targets,'bands':bandrows,'corners':corners})
  require(counters['point_enclosures']==5970 and len(rows)==128,'COMPLETE_FIXED_SOURCE_SIZE');summary={}
  for t in cfg['targets']:
   summary[t]={'true_T_box_paid':sum(b['targets'][t]['classification']=='RATIONAL_TRUE_T_BOX_PAID'for b in rows),'old_odd15_paid':sum(b['targets'][t]['old_odd15_box_pass']for b in rows),'degree15_loss_newly_paid':sum(b['targets'][t]['newly_paid_after_stronger_rational_enclosure']for b in rows),'true_fixed_anchor_wholebox_deficit':sum(b['targets'][t]['classification']=='FIXED_ANCHOR_WHOLEBOX_TRUE_T_DEFICIT'for b in rows),'rational_enclosure_unresolved':sum(b['targets'][t]['classification']=='RATIONAL_ENCLOSURE_UNRESOLVED'for b in rows),'boxes_with_ordered_point_family_deficit':sum(b['targets'][t]['ordered_fixed_point_family_bad_corners']>0 for b in rows),'boxes_with_only_rectangle_point_deficit':sum(b['targets'][t]['only_full_rectangle_bad_corner_scope']for b in rows)}
  proposals=make_anchor_proposals(cfg,rows,universe,statuses);write_new('rational_true_T_point_enclosures.json',{'entries':bounds,'count':len(bounds),'actual_counts':counters,'native_calls':0,'transcendental_exp_calls':0});write_new('shared_anchor_proposals.json',proposals);result={'status':'SOLE_SAVED_RATIONAL_DIAGNOSIS_COMPLETE','targets':summary,'box_count':128,'original_odd15_unresolved':102,'boxes':rows,'actual_counts':counters,'anchor_proposal_count':proposals['proposed_anchor_count'],'native_calls':0,'transcendental_exp_calls':0,'checker_calls':0,'subdivisions':0,'wall_seconds':time.monotonic()-t0,'global_g_certificate':False,'new_C6_bound':False,'upper_scope':'fixed robust tangent/certificate surrogate only; never actual G upper or true g counterexample'};write_new('result.json',result);write_new('exit_receipt.json',{'status':result['status'],'exit_code':0,'result_sha256':sha(P/'result.json'),'enclosures_sha256':sha(P/'rational_true_T_point_enclosures.json'),'proposals_sha256':sha(P/'shared_anchor_proposals.json'),'actual_counts':counters,'native_calls':0,'exp_calls':0});print(json.dumps({k:v for k,v in result.items()if k!='boxes'},sort_keys=True))
 except BaseException:atomic('failure_receipt.json',{'status':'FAILED_STOP_NO_RETRY','traceback':traceback.format_exc(),'actual_counts':counters,'native_calls':0,'exp_calls':0});raise
if __name__=='__main__':
 if sys.argv[1:]==['--source-check-only']:
  cfg,result,cache,universe=source_preflight();print(json.dumps({'status':'SOURCE_ONLY_PASS','points':len(cache),'boxes':len(result['boxes']),'residual_universe':len(universe),'new_true_T_enclosures':0,'native_calls':0,'exp_calls':0}))
 elif not sys.argv[1:]:run()
 else:raise SystemExit('usage: run_batch.py [--source-check-only]')
