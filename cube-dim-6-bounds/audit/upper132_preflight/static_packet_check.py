import sys,json,hashlib,math,heapq,ast
from pathlib import Path
from fractions import Fraction as Q
sys.set_int_max_str_digits(0);sys.dont_write_bytecode=True
P=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-numerics/cube-dim-6-bounds/upper/round132/cross_center_band_envelopes');O=Path(__file__).resolve().parent
rd=lambda p:json.loads(Path(p).read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
canon=lambda x:(json.dumps(x,sort_keys=True,separators=(',',':'),ensure_ascii=False)+'\n').encode()
H=rd(P/'input_hashes.json');cfg=rd(P/'selected_boxes_and_centers.json');selected=rd(P/'selected_source_cache_entries.json')['entries'];condition=rd(P/'conditional_M6_check.json')
checks={}
checks['all_frozen_sources_hash']=all(sha(p)==h for p,h in H['frozen_sources'].items())
checks['all_selected_original_artifact_hash']=all(sha(p)==h for p,h in H['selected_original_artifacts'].items())
cachepath=next(Path(p)for p in H['frozen_sources']if p.endswith('/combined_cache.json'));seedpath=next(Path(p)for p in H['frozen_sources']if p.endswith('/continuation_seed129.json'));cache=rd(cachepath);seed=rd(seedpath)
checks['frozen_cache15193_static_seed_only']=len(cache['entries'])==15193 and '/round129/'in str(seedpath) and seedpath.name=='continuation_seed129.json'
rid=next(k for k in sorted(seed['roots'])if not seed['roots'][k]['root_complete']);s=seed['roots'][rid];ids=sorted(set(s['queue']+s['blocked_leaves']))
checks['first_lex_unclosed_root_all93_boxes']=rid==cfg['root_id']=='C121.r11011101101' and len(ids)==cfg['source_unresolved_leaf_count']==cfg['selected_box_count']==93 and [x['node_id']for x in cfg['selected_boxes']]==ids and cfg['all_root_unresolved_leaves_selected']is True
windows=[{'kind':'B','a':a,'T':t}for a,t in [('1/32','1/2'),('1/64','1/32'),('1/128','1/64'),('1/250','1/128'),('1/1000','1/250')]]
checks['five_positive_disjoint_windows']=cfg['windows']==windows and len({json.dumps(w,sort_keys=True)for w in windows})==5 and all(0<Q(w['a'])<Q(w['T'])for w in windows)and all(Q(windows[j+1]['T'])==Q(windows[j]['a'])for j in range(4))
checks['every_box_five_windows_exactly_once']=all([b['window']for b in x['bands']]==windows for x in cfg['selected_boxes'])
checks['all_boxes_static_original_unresolved']=all(x['box']==s['nodes'][x['node_id']]['box'] and x['original_status']==s['nodes'][x['node_id']]['status'] and x['node_id']not in s['passed_leaves']for x in cfg['selected_boxes'])
groups=[[(k,m,[Q(z)for z in m['expected']['center'][1:]])for k,m in cache['entries'].items()if m['expected']['window']==w]for w in windows];counts=list(map(len,groups))
checks['frozen_candidate_counts']=counts==cfg['candidate_counts_by_window']==[14670,523,0,0,0]
centers=[[sum(map(Q,pair))/2 for pair in x['box']]for x in cfg['selected_boxes']];scale=math.lcm(*(z.denominator for group in groups for _,_,c in group for z in c),*(z.denominator for c in centers for z in c));lattice=[[(k,m,[int(z*scale)for z in c])for k,m,c in group]for group in groups]
nearest_ok=True;refs_ok=True;total=0;used=set()
for x,c in zip(cfg['selected_boxes'],centers):
 c=[int(z*scale)for z in c]
 for wi,band in enumerate(x['bands']):
  rows=band['candidate_keys'];keys=[r['semantic_key']for r in rows];refs_ok &= len(keys)==len(set(keys))<=16 and all(k in selected and selected[k]==cache['entries'][k] and selected[k]['expected']['window']==windows[wi]for k in keys)
  nearest=heapq.nsmallest(16,((sum((a-b)**2 for a,b in zip(v,c)),k,m)for k,m,v in lattice[wi]),key=lambda y:(y[0],y[1]))
  expected=[{'semantic_key':k,'distance_squared':str(Q(d,scale*scale)),'true_center':m['expected']['center']}for d,k,m in nearest]
  nearest_ok &= rows==expected;refs_ok &= band['missing_band_lower']==('0'if not rows else None);used.update(keys);total+=32*len(rows)
checks['every_band_unique_le16_same_window_and_missing0']=refs_ok
checks['all_exact_nearest16_and_distances']=nearest_ok
checks['325_sources_exact_frozen_cache_subset']=used==set(selected) and len(used)==cfg['selected_unique_jet_count']==325
checks['actual95232_absolute327680_exp_counts']=total==cfg['actual_planned_corner_exp']==95232 and cfg['max_corner_exp']==128*5*16*32==327680 and total<=327680
# Decode selected originals as exact dyadic metadata only: no evaluator/native/exp call.
def dyadic(v):return[Q(int(v[key]))*Q(2)**int(v['common_binary_exponent'])for key in ['lower_mantissa','upper_mantissa']]
raw_ok=all_center_ok=receipt_ok=True
for k,m in selected.items():
 e=m['expected'];jet=rd(m['jet_path']);raw=rd(m['native_path']);receipt=rd(m['receipt_path']);pr=jet['provenance']
 all_center_ok &= list(map(Q,e['center']))==list(map(Q,jet['center'])) and Q(e['center'][0])==1 and jet['window']==e['window'] and jet['model']==e['model'] and jet['model_definition']==e['model_definition']
 raw_ok &= list(map(Q,jet['F_interval']))==dyadic(raw['band']) and [list(map(Q,v))for v in jet['face_gradient_intervals']]==[dyadic(v)for v in raw['gradient'][1:]] and len(raw['gradient'])==6 and len(jet['face_gradient_intervals'])==5
 raw_ok &= dyadic(raw['band'])[0]>0 and raw['status']=='all_real_nodes_supported' and raw['precision_bits']==m['precision_bits']==pr['precision_bits'] and raw['full_group_terms']==967680 and raw['computed_sign_pairs']==483840 and pr['native_sha256']==m['native_sha256']==sha(m['native_path']) and m['jet_sha256']==sha(m['jet_path'])
 receipt_ok &= receipt['exit_code']==0 and m['receipt_sha256']==sha(m['receipt_path'])
checks['all325_original_true_center_face_and_window_preserved']=all_center_ok
checks['all325_native_value_full5grad_dyadic_match']=raw_ok
checks['all325_auth_original_receipt_exit0_hash']=receipt_ok
m6=Q(rd(condition['source_path'])['M_power6_upper']);g=Q(1,63500);u5=Q(9,8)**6*m6/(384**5*g**6)
checks['complete128_exactM6_conditional_formula']=Q(condition['M6_exact'])==m6 and condition['source_sha256']==sha(condition['source_path']) and Q(condition['U_power5_exact'])==u5 and Q(condition['rhs_minus_U_power5_exact'])==64**5-u5 and u5<64**5 and condition['strict_less_than_64_power5']is True
checks['conditional_not_global_scope']=condition['global_g_floor_assumed_not_proved']is True and condition['global_g_certificate']is False and condition['new_C6_bound']is False
checks['AST_all_python']=all(ast.parse(x.read_text())is not None for x in P.glob('*.py'))
result={'status':'PASS_STATIC_METADATA_SOURCE_AND_CONDITIONAL_FRACTION'if all(checks.values())else'FAIL','checks':checks,'selected_sources':325,'selected_original_artifact_hash_count':len(H['selected_original_artifacts']),'exp_native_old_baseline_evaluator_calls':0,'scope':'frozen metadata/hash/dyadic/Fraction only; no full15193 baseline authentication replay','input_hashes_sha256':sha(P/'input_hashes.json'),'selected_boxes_sha256':sha(P/'selected_boxes_and_centers.json'),'selected_sources_sha256':sha(P/'selected_source_cache_entries.json')}
(O/'static_packet_review.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
