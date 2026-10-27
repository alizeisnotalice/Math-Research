#!/usr/bin/env python3
"""Independent finite C145 saved-only Fraction review. No producer imports or exp."""
import hashlib, itertools, json, math, os, sys, time
from pathlib import Path
from fractions import Fraction as F
sys.set_int_max_str_digits(0)
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
SRC=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-numerics/cube-dim-6-bounds/upper/round145/rational_deficit_and_shared_anchor_plan')
N=0
def ck(x,m):
 global N
 N+=1
 if not x: raise AssertionError(m)
def read(p): return json.loads(Path(p).read_text())
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ordered(p):
 p=list(map(F,p)); return p[0]==1 and p[-1]>=0 and all(a>=b for a,b in zip(p,p[1:]))
def canon(x):return (json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
def rounded(x,up=False):
 s=1<<256; a=x.numerator*s; b=x.denominator
 return F(-((-a)//b) if up else a//b,s)
def polynomial(u,n):
 # Independent Horner evaluation of factorial coefficients, not producer recurrence.
 p=F(1,math.factorial(n))
 for k in range(n-1,-1,-1):p=p*u+F(1,math.factorial(k))
 return p
def interval(r):
 m=0
 while abs(r)>F(1,2):r/=2;m+=1
 ck(m<=16,'range cap');u=-abs(r)
 lo,hi=polynomial(u,15),polynomial(u,16)
 ck(0<lo<=hi,'Taylor sign/positivity')
 if r>0:lo,hi=1/hi,1/lo
 lo,hi=rounded(lo),rounded(hi,True)
 for k in range(m):lo,hi=rounded(lo**2),rounded(hi**2,True)
 ck(0<=lo<=hi,'dyadic order')
 return lo,hi,m
def main():
 t0=time.monotonic();cpu=time.process_time()
 os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
 ck(not (HERE/'receipt.json').exists(),'one saved review only')
 localnames=['PROOF.md','run_batch.py','prepare.py','inputs.json','input_hashes.json','registration.json','release_requirements.json','parent_release.json','prepared_packet.json','existing_semantic_statuses.json','exit_receipt.json','result.json','rational_true_T_point_enclosures.json','shared_anchor_proposals.json']
 pins={str(SRC/n):sha(SRC/n) for n in localnames}
 frozen=read(SRC/'input_hashes.json')['frozen_inputs']
 for p,h in frozen.items():ck(sha(p)==h,'frozen input hash:'+p)
 pins.update(frozen)
 cfg=read(SRC/'inputs.json');old=read(cfg['C141_result_path']);cache=read(cfg['C141_rational_cache_path'])['entries'];universe=read(cfg['C129_residual_metadata_path'])['leaves']
 result=read(SRC/'result.json');bounds=read(SRC/'rational_true_T_point_enclosures.json');proposals=read(SRC/'shared_anchor_proposals.json');receipt=read(SRC/'exit_receipt.json')
 for name,key in [('result.json','result_sha256'),('rational_true_T_point_enclosures.json','enclosures_sha256'),('shared_anchor_proposals.json','proposals_sha256')]:ck(sha(SRC/name)==receipt[key],'exit binding')
 ck(receipt['exit_code']==0,'source exit0');reg=read(SRC/'registration.json')
 for n,k in [('run_batch.py','program_sha256'),('PROOF.md','proof_sha256'),('inputs.json','inputs_sha256'),('input_hashes.json','input_hashes_sha256'),('existing_semantic_statuses.json','semantic_statuses_sha256')]:ck(sha(SRC/n)==reg[k],'registration pin')
 ck(len(cache)==len(bounds['entries'])==5970 and len(old['boxes'])==len(result['boxes'])==128 and len(universe)==10266,'finite cardinalities')
 ck(set(cache)==set(bounds['entries']),'all saved point IDs')
 windows=[(F(w['a']),F(w['T'])) for w in cfg['source_windows']]
 ck(len(set(windows))==5 and all(0<a<b for a,b in windows),'positive five distinct bands')
 ck(all(windows[i][0]==windows[i+1][1] for i in range(4)),'disjoint adjacent windows')
 B={};steps=0;maxsteps=0
 for key,row in cache.items():
  f=F(row['clamped_true_value']);r=F(row['exponent_exact']);a=F(row['own_time_a']);q=F(row['frame32_quadratic']);b=F(row['gradient_dot_lower'])
  ck(f>0 and r==b/f-q/a**2,'source exact formula')
  lo,hi,m=interval(r);steps+=m;maxsteps=max(m,maxsteps);v=bounds['entries'][key]
  lower=max(f*lo,F(row['band_corner_lower_exact']));upper=f*hi
  ck(lower<=upper,'existing old lower contained')
  ck((lo,hi,lower,upper,m)==(F(v['exp_lower_exact']),F(v['exp_upper_exact']),F(v['true_T_lower_exact']),F(v['true_T_upper_exact']),v['range_squarings']),'all new rational endpoints')
  ck(v['point']==row['point'] and v['semantic_key']==row['semantic_key'],'source point provenance')
  B[key]=(lower,upper)
 counts={'Taylor_updates':5970*16,'dyadic_squaring_steps':steps,'point_enclosures':5970}
 ck(counts==result['actual_counts']==bounds['actual_counts']==receipt['actual_counts'],'actual budget counts')
 ck(steps==23846 and maxsteps==11 and steps<=95520,'squaring totals')
 targets=['1/63500','1/65000'];computed=[];bitsset=set(itertools.product((0,1),repeat=5))
 ck(len({b['node_id'] for b in result['boxes']})==128,'unique box IDs')
 for source,box in zip(old['boxes'],result['boxes']):
  ck(all(box[k]==source[k] for k in ['node_id','root_id','box']),'original box identity')
  ck(len(source['per_band'])==len(box['bands'])==5,'exact band count')
  points={bits:['1']+[pair[bit] for pair,bit in zip(box['box'],bits)] for bits in bitsset}
  pointlo={bits:F(0) for bits in bitsset};pointhi=pointlo.copy();L=F(0);U=F(0)
  for j,(band,saved) in enumerate(zip(source['per_band'],box['bands'])):
   ck(band['window']==cfg['source_windows'][j]==saved['window'],'all exact windows')
   refs=band['corner_references'];ck(len(refs)==32 and {tuple(r['corner_bits']) for r in refs}==bitsset,'complete 32 corner partition')
   lows=[];highs=[]
   for ref in refs:
    bits=tuple(ref['corner_bits']);row=cache[ref['cache_key']];ck(row['point']==points[bits] and row['semantic_key']==band['semantic_key'] and row['window']==band['window'],'actual physical corner identity')
    l,u=B[ref['cache_key']];lows.append(l);highs.append(u);pointlo[bits]+=l;pointhi[bits]+=u
   l,u=min(lows),min(highs);L+=l;U+=u
   ck((l,u)==(F(saved['true_T_wholebox_lower_exact']),F(saved['true_T_wholebox_upper_exact'])),'bandwise min then sum')
  ck(len(box['corners'])==32 and {tuple(c['corner_bits']) for c in box['corners']}==bitsset,'complete output corners')
  for c in box['corners']:
   bits=tuple(c['corner_bits']);ck((c['point'],c['ordered'],F(c['fixed_point_family_lower_exact']),F(c['fixed_point_family_upper_exact']))==(points[bits],ordered(points[bits]),pointlo[bits],pointhi[bits]),'point same-location sum and ordering')
  for target in targets:
   T=F(target);v=box['targets'][target];paid=L>=T;deficit=U<T;oldpaid=source['passes_gstar' if target==targets[0] else 'passes_g0']
   classification='RATIONAL_TRUE_T_BOX_PAID' if paid else 'FIXED_ANCHOR_WHOLEBOX_TRUE_T_DEFICIT' if deficit else 'RATIONAL_ENCLOSURE_UNRESOLVED'
   bad=[bits for bits in bitsset if pointhi[bits]<T];obad=[bits for bits in bad if ordered(points[bits])]
   expected={'classification':classification,'old_odd15_box_pass':oldpaid,'newly_paid_after_stronger_rational_enclosure':paid and not oldpaid,'true_T_box_lower_exact':str(L),'true_T_box_upper_exact':str(U),'box_upper_minus_target_exact':str(U-T),'fixed_point_family_bad_corners':len(bad),'ordered_fixed_point_family_bad_corners':len(obad),'only_full_rectangle_bad_corner_scope':bool(bad) and not obad,'all_point_corner_lowers_paid':all(x>=T for x in pointlo.values()),'true_g_counterexample':False,'actual_G_upper_bound':False}
   ck(v==expected,'all target classification fields');ck(not oldpaid or paid,'old paid preserved')
  computed.append(box)
 for target in targets:
  rows=[x['targets'][target] for x in computed]
  expected={'true_T_box_paid':sum(x['classification']=='RATIONAL_TRUE_T_BOX_PAID' for x in rows),'old_odd15_paid':sum(x['old_odd15_box_pass'] for x in rows),'degree15_loss_newly_paid':sum(x['newly_paid_after_stronger_rational_enclosure'] for x in rows),'true_fixed_anchor_wholebox_deficit':sum(x['classification']=='FIXED_ANCHOR_WHOLEBOX_TRUE_T_DEFICIT' for x in rows),'rational_enclosure_unresolved':sum(x['classification']=='RATIONAL_ENCLOSURE_UNRESOLVED' for x in rows),'boxes_with_ordered_point_family_deficit':sum(x['ordered_fixed_point_family_bad_corners']>0 for x in rows),'boxes_with_only_rectangle_point_deficit':sum(x['only_full_rectangle_bad_corner_scope'] for x in rows)}
  ck(expected==result['targets'][target],'target aggregate')
 # Independently reconstruct the prior frozen semantic registry, without native artifacts.
 statuses={}
 for path in frozen:
  if path.endswith('cache_crosswalk_12081.json'):statuses.update({k:'AVAILABLE' for k in read(path)['entries']})
 for filename in ['recovered_attempt_index.json','semantic_runtime_index_v129.json']:
  path=next(p for p in frozen if p.endswith(filename))
  for k,v in read(path)['entries'].items():
   ck(v['status'] in ['AVAILABLE','UNKNOWN_NO_REISSUE'],'old semantic status');ck(k not in statuses or statuses[k]==v['status'],'registry consistency');statuses[k]=v['status']
 selection=read(Path(cfg['C141_result_path']).parent/'selected_boxes_and_existing_jets.json')
 for k in selection['source_five_semantic_keys_no_reissue']:statuses[k]='AVAILABLE_R137'
 ck(statuses==read(SRC/'existing_semantic_statuses.json')['entries'] and len(statuses)==33470,'whole frozen registry')
 def semantic(point,w):
  return hashlib.sha256(canon({'center':list(map(lambda x:str(F(x)),point)),'a':str(F(w['a'])),'T':str(F(w['T'])),'formula_family_id':'C016-cauchy-frame32-outward-precision-family-v1','fan_sha256':cfg['fan_sha256'],'model':cfg['source_model'],'model_definition':cfg['source_model_definition'],'quadratic':'frame32'})).hexdigest()
 paid={b['node_id'] for b in computed if b['targets'][targets[0]]['classification']=='RATIONAL_TRUE_T_BOX_PAID'}
 ulookup={b['node_id']:b for b in universe};ck(len(ulookup)==10266,'unique residual IDs')
 for box in computed:
  ck(box['node_id'] in ulookup and box['box']==ulookup[box['node_id']]['box'],'local boxes exact residual inclusion')
 ck(sorted(paid)==proposals['locally_paid_node_ids_excluded_for_planning_only'],'paid planning exclusion')
 pool={};anchor=['1','1','1','3/8','1/4','1/8']
 for box in computed:
  bad=[c for c in box['corners'] if c['ordered'] and F(c['fixed_point_family_upper_exact'])<F(1,63500)]
  if not bad:continue
  witness=min(bad,key=lambda c:(F(c['fixed_point_family_upper_exact']),tuple(map(F,c['point']))));gap=F(1,63500)-F(witness['fixed_point_family_upper_exact']);p=witness['point']
  for point in [['1','1']+p[2:],['1','1','1']+p[3:]]:
   point=list(map(lambda x:str(F(x)),point));key=tuple(point)
   if point==anchor:continue
   ck(ordered(point),'proposal ordered projection')
   entry=pool.setdefault(key,{'point':point,'justifying_original_boxes':[],'weight':F(0)})
   if box['node_id'] not in entry['justifying_original_boxes']:entry['justifying_original_boxes'].append(box['node_id']);entry['weight']+=gap
 candidates=sorted(pool.values(),key=lambda c:(-c['weight'],tuple(map(F,c['point']))))[:64]
 ck(len(pool)==proposals['candidate_pool_before_cap'] and len(candidates)==proposals['candidate_pool_capped'],'complete ranked candidate pool')
 eligible=[];excluded=[]
 for c in candidates:
  queries=[]
  for w in cfg['source_windows']:
   k=semantic(c['point'],w);queries.append({'window':w,'semantic_key':k,'frozen_status':statuses.get(k,'NEW_IN_FROZEN_SNAPSHOT_ONLY'),'reissue_allowed':False if k in statuses else None})
  if any(q['frozen_status']=='UNKNOWN_NO_REISSUE' for q in queries):excluded.append({'point':c['point'],'reason':'PRIOR_UNKNOWN_NO_REISSUE','queries':queries});continue
  c['queries']=queries;c['new_semantics_in_frozen_snapshot']=sum(q['frozen_status']=='NEW_IN_FROZEN_SNAPSHOT_ONLY' for q in queries);eligible.append(c)
 ck(excluded==proposals['unknown_anchor_candidates_excluded'],'unknown quarantine')
 # Exact integer scaling independently evaluates the max-corner curvature proxy.
 D=math.lcm(*(F(v).denominator for b in universe for pair in b['box'] for v in pair),*(F(v).denominator for c in eligible for v in c['point'][1:]))
 iboxes=[[(int(F(l)*D),int(F(h)*D)) for l,h in b['box']] for b in universe];pairs=0
 for c in eligible:
  center=[int(F(x)*D) for x in c['point'][1:]];covered=set()
  for i,box in enumerate(iboxes):
   if universe[i]['node_id'] in paid:continue
   pairs+=1
   if sum(max(abs(l-x),abs(h-x))**2 for (l,h),x in zip(box,center))*128<=D**2:covered.add(i)
  c['cover']=covered
 ck(pairs==proposals['actual_geometry_box_pairs']==655360 and pairs<=657024,'geometry cap')
 deficits={b['node_id']:max(F(0),F(1,63500)-F(b['targets'][targets[0]]['true_T_box_upper_exact'])) for b in computed if b['targets'][targets[0]]['classification']=='FIXED_ANCHOR_WHOLEBOX_TRUE_T_DEFICIT'}
 remaining=eligible.copy();covered=set();selected=[]
 while remaining and len(selected)<8:
  def rank(c):
   new=c['cover']-covered;return(-len(new),-sum(deficits.get(universe[i]['node_id'],F(0)) for i in new),c['new_semantics_in_frozen_snapshot'],tuple(map(F,c['point'])))
  c=min(remaining,key=rank);remaining.remove(c);new=c['cover']-covered
  if not new:break
  covered|=c['cover'];v={'point':c['point'],'justifying_original_boxes':c['justifying_original_boxes'],'known_fixed_anchor_point_deficit_weight':str(c['weight']),'queries':c['queries'],'new_semantics_in_frozen_snapshot':c['new_semantics_in_frozen_snapshot'],'geometry_proxy_newly_covered_boxes':len(new),'geometry_proxy_all_covered_boxes':len(c['cover']),'geometry_proxy_root_count':len({universe[i]['root_id'] for i in c['cover']}),'geometry_proxy_covered_node_ids':[universe[i]['node_id'] for i in sorted(c['cover'])],'scientific_query_authorized':False,'unknown_new_source_values_and_gradients':True};selected.append(v)
 ck(selected==proposals['proposals'] and len(selected)==proposals['proposed_anchor_count']==result['anchor_proposal_count']==8,'all greedy selected fields')
 new=sum(c['new_semantics_in_frozen_snapshot'] for c in selected)
 ck(new==40==proposals['new_semantics_in_frozen_snapshot'] and 3*new==120==proposals['max_future_native_attempts_proposal_only'],'future proposal cap only')
 ck(all(result[k]==0 for k in ['native_calls','transcendental_exp_calls','checker_calls','subdivisions']) and result['global_g_certificate'] is False and result['new_C6_bound'] is False,'no global or oracle claim')
 for p,h in pins.items():ck(sha(p)==h,'after review frozen hash:'+p)
 out={'status':'PASS','assertions':N,'independent_method':'Fraction Horner P15/P16, exact outward dyadic squaring; band-min then sum; full ordered witness and greedy geometry reconstruction','source_path':str(SRC),'point_enclosures':5970,'max_squarings':maxsteps,'actual_counts':counts,'boxes':128,'physical_corners':4096,'targets':result['targets'],'frozen_semantic_keys':len(statuses),'prior_unknown_count':sum(v=='UNKNOWN_NO_REISSUE' for v in statuses.values()),'candidate_pool':len(pool),'geometry_pairs':pairs,'anchor_proposals':8,'shared_geometry_union_boxes':len(covered),'native_calls':0,'transcendental_exp_calls':0,'producer_calls':0,'old_baseline_recomputations':0,'wall_seconds':time.monotonic()-t0,'CPU_seconds':time.process_time()-cpu,'nice':os.getpriority(os.PRIO_PROCESS,0)}
 (HERE/'source_pins.json').write_text(json.dumps(pins,indent=2,sort_keys=True)+'\n')
 (HERE/'receipt.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
