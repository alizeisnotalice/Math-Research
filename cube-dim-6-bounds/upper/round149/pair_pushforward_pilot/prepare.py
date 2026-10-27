#!/usr/bin/env python3
"""Read-only Fraction selection and source qualification; no density generation."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict,Counter
import hashlib,json,os,sys,time
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent
TERRA=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-brascamp108/cube-dim-6-bounds/coordination/terra_20261027/pair_pushforward_patterns/patterns.json')
C130=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round130/true_second_count_boxes')
C147=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-localmax-audit/cube-dim-6-bounds/upper/round147/pair_pushforward')
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(n,x):
 with (P/n).open('x')as f:json.dump(x,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def canonical(cols):
 cols=[tuple(map(F,v))for v in cols];normal=tuple(sorted(cols));swapped=tuple(sorted((b,a)for a,b in cols));return min(normal,swapped)
def serial(cols):return [list(map(str,v))for v in cols]
def key(cols):return hashlib.sha256((json.dumps(serial(cols),separators=(',',':'))+'\n').encode()).hexdigest()
def features(cols):
 active=[v for v in cols if any(v)];directions=Counter(tuple(x/sum(v)for x in v)for v in active);multiplicity=Counter(active)
 return {'active_columns':len(active),'distinct_directions':len(directions),'zero_columns':6-len(active),'nonzero_repeated_column_excess':len(active)-len(multiplicity),'maximum_identical_nonzero_column_multiplicity':max(multiplicity.values()),'maximum_direction_multiplicity':max(directions.values()),'has_axis_direction':any(0 in v for v in active)}
def main():
 start=time.monotonic();cpu=time.process_time()
 assert not any((P/n).exists()for n in ['case_packet.json','source_manifest.json','registration.json','prepare_receipt.json','parent_release.json','started.json','result.json'])
 assert sha(TERRA)=='01a05b2633dfdc7994227b24618847e4b1e7da98b0cde0704acf72ec9753876a'
 for frozen_path,frozen_hash in {'/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round130/true_second_count_boxes/cone_table.json': '86f81fd6c69d8c9a6b7194d34201732c528ac4cde320aefbb19570d854097040', '/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round130/true_second_count_boxes/box_core.py': 'b1e16d942bcec857562c900bda301e817cee741090b282bf8a6fdfbeedc70275', '/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-localmax-audit/cube-dim-6-bounds/upper/round147/pair_pushforward/PROOF.md': '2ecb770367e181fa153c7cfc1d2bcb7c67f7babc17b6ed7b24b7f4c78adeb22f', '/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-localmax-audit/cube-dim-6-bounds/upper/round147/pair_pushforward/manifest.json': '2500cf19538335300f9c4f2db74c914620acba279af8e182a1f43c96eb9de69f', '/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-localmax-audit/cube-dim-6-bounds/upper/round147/pair_pushforward/prototype_result.json': '5a492e403113d73aae9665735f0da68f8b501416e0c42548b0fae69e4b103f5b', '/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-localmax-audit/cube-dim-6-bounds/upper/round147/pair_pushforward/exit_receipt.json': 'dd071ca60efd9d88a612368e06c089e15760ee803accd23ad0723450187bdb50', '/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-localmax-audit/cube-dim-6-bounds/upper/round147/pair_pushforward/next_result.json': '1233ffbd6c9b093f3e99715affe3ef274e4d7fff374e23f9f3aa4c2c9889acd8', '/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-localmax-audit/cube-dim-6-bounds/upper/round147/pair_pushforward/next_exit_receipt.json': 'c5d44067b73d8cddfac3420aa6290b0406db7b257443cfe2ec5456a798b5063f'}.items():assert sha(frozen_path)==frozen_hash
 patterns=read(TERRA)['patterns'];cones=read(C130/'cone_table.json')['cones'];assert len(patterns)==9961
 records={};all_sources=set();weight=F(0)
 for i,row in enumerate(patterns):
  cols=canonical(row['columns']);k=key(cols);assert k not in records and row['rank']==2
  w=F(0)
  for cid,a,b in row['sources']:
   assert 0<=cid<21 and 0<=a<b<32 and (cid,a,b)not in all_sources;all_sources.add((cid,a,b))
   original=list(zip(cones[cid]['nodes_by_body'][a],cones[cid]['nodes_by_body'][b]));assert canonical(original)==cols;w+=F(cones[cid]['alpha'])
  assert w==F(row['weight']) and len(row['sources'])==row['source_count'];weight+=w
  f=features(cols);assert 2<=f['active_columns']<=6 and 2<=f['distinct_directions']<=f['active_columns']
  records[k]={'canonical_key':k,'terra_pattern_index':i,'columns':serial(cols),'pattern_weight':str(w),'source_count':row['source_count'],'sources':row['sources'],'geometry':f}
 assert len(all_sources)==10416 and weight==496*sum((F(c['alpha'])for c in cones),F(0))
 old=[('prototype',0,16,17,C147/'prototype_result.json'),('m3',0,9,16,C147/'next_result.json'),('m6',0,0,1,C147/'next_result.json')];excluded=[]
 toy=read(C147/'prototype_result.json');nxt=read(C147/'next_result.json');assert toy['status']=='PASS_ONE_REAL_PAIR_EXACT_POSITIVE_2D_PUSHFORWARD' and nxt['status']=='COMPLETE_TWO_FIXED_EXACT_DENSITY_GENERATIONS'
 assert read(C147/'exit_receipt.json')['result_sha256']==sha(C147/'prototype_result.json') and read(C147/'next_exit_receipt.json')['result_sha256']==sha(C147/'next_result.json')
 for label,cid,a,b,path in old:
  raw=[tuple(map(F,v))for v in zip(cones[cid]['nodes_by_body'][a],cones[cid]['nodes_by_body'][b])];cols=canonical(raw);k=key(cols);record=records[k]
  swapped=tuple(sorted((y,x)for x,y in raw))==cols and tuple(sorted(raw))!=cols
  excluded.append({'passed_case':label,'canonical_key':k,'source_reference_path':str(path),'source_reference_sha256':sha(path),'original_cone_id':cid,'original_body_pair_zero_based':[a,b],'canonical_requires_output_coordinate_swap':swapped,'column_permutation_only':True,'kernel_product_is_symmetric_under_output_swap':True,'pattern_weight':record['pattern_weight'],'source_count':record['source_count'],'generation_repeated':False})
 excludedkeys={v['canonical_key']for v in excluded};assert len(excludedkeys)==3
 groups=defaultdict(list)
 for k,row in records.items():
  if k in excludedkeys:continue
  f=row['geometry'];groups[(f['active_columns'],f['distinct_directions'],f['nonzero_repeated_column_excess']>0)].append(k)
 for g in groups:groups[g].sort()
 selected=[];seen=set();positions={g:0 for g in groups}
 # One of every available (m,d,nonzero-repeat) stratum first, then deterministic round robin.
 order=sorted(groups,key=lambda g:(-g[1],-g[0],-int(g[2])))
 while len(selected)<100:
  progress=False
  for g in order:
   while positions[g]<len(groups[g])and groups[g][positions[g]]in seen:positions[g]+=1
   if positions[g]==len(groups[g]):continue
   k=groups[g][positions[g]];positions[g]+=1;seen.add(k);selected.append(k);progress=True
   if len(selected)==100:break
  assert progress
 assert not seen&excludedkeys and len(seen)==100
 cases=[]
 for i,k in enumerate(selected):cases.append({'case_id':'case_'+str(i).zfill(3),**records[k]})
 witness=next(v['canonical_key']for v in cases if v['geometry']['active_columns']==3 and v['geometry']['distinct_directions']==2 and v['geometry']['nonzero_repeated_column_excess']==1)
 assert {v['geometry']['active_columns']for v in cases}==set(range(2,7)) and max(v['geometry']['distinct_directions']for v in cases)==6
 caps={'density_stage_checks_per_case':5,'nonzero_terms_per_stage':2000,'literal_endpoint_substitutions_per_case':2200,'polynomial_monomials_per_stage':30000,'Laplace_monomial_summands_per_stage':30000,'expression_nodes_per_checked_expression':1000000,'checked_expression_node_visits_per_case':20000000,'rational_coefficient_bits':8192}
 packet={'status':'PREPARED_UNRUN_ONLY_NEW_CANONICAL_DENSITIES','terra_pattern_path':str(TERRA),'terra_pattern_sha256':sha(TERRA),'cone_table_path':str(C130/'cone_table.json'),'cone_table_sha256':sha(C130/'cone_table.json'),'canonicalization':'sort exact rational six columns; choose smaller of direct vs simultaneous coordinate swap; zero columns retained in key','all_canonical_keys':9961,'all_source_pairs':10416,'all_pattern_weight':str(weight),'selected_cases':cases,'excluded_passed_references':excluded,'coverage_strata':{str(g):sum(v['geometry']['active_columns']==g[0]and v['geometry']['distinct_directions']==g[1]and(v['geometry']['nonzero_repeated_column_excess']>0)==g[2]for v in cases)for g in order},'resonance_witness_key':witness,'selected_pattern_weight':str(sum((F(v['pattern_weight'])for v in cases),F(0))),'selected_source_pair_count':sum(v['source_count']for v in cases),'max_workers':15,'caps':caps,'new_k_integral_evaluations':0,'new_beta2_upper':False,'new_C6_bound':False}
 write('case_packet.json',packet)
 paths=[TERRA,C130/'cone_table.json',C130/'box_core.py',C147/'PROOF.md',C147/'manifest.json',C147/'prototype_result.json',C147/'exit_receipt.json',C147/'next_result.json',C147/'next_exit_receipt.json']
 write('source_manifest.json',{'sources':[{'path':str(v),'sha256':sha(v),'bytes':v.stat().st_size}for v in paths],'read_live129':False,'old_6D_grid_calls_permitted':0})
 reg={'status':'PREPARED_UNRUN_PENDING_DUAL_GATE','new_canonical_cases':100,'max_workers':15,'workers_set_by_main_release':True,'max_stages_total':sum(v['geometry']['active_columns']-1 for v in cases),'max_literal_endpoint_substitutions_total':220000,'caps_per_case':caps,'counter_semantics':'literal endpoints = 2*(resonant_branches + nonresonant_branches); separate branch counters','sympy_version':'1.14.0','time_limit':None,'nice_minimum':10,'automatic_restarts':0,'k_integral_native_special_function_calls':0,'case_packet_sha256':sha(P/'case_packet.json'),'source_manifest_sha256':sha(P/'source_manifest.json'),'run_sha256':sha(P/'run_prepared.py'),'core_sha256':sha(P/'density_core.py'),'prepare_sha256':sha(__file__),'proof_sha256':sha(P/'PROOF.md')}
 write('registration.json',reg)
 write('prepare_receipt.json',{'status':'PREPARED_ONLY_NO_DENSITY_GENERATIONS','new_cases':100,'excluded_passed_keys':3,'new_density_generations':0,'new_k_integral_special_function_native_kernel_calls':0,'old_6D_grid_replays':0,'script_sha256':sha(__file__),'case_packet_sha256':sha(P/'case_packet.json'),'registration_sha256':sha(P/'registration.json'),'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu})
 print(json.dumps({'status':'PREPARED_UNRUN','new_cases':100,'excluded_passed_keys':3,'source_pairs_in_new_cases':packet['selected_source_pair_count'],'stage_cap_total':reg['max_stages_total'],'coverage_strata':packet['coverage_strata'],'new_density_generations':0}))
if __name__=='__main__':main()
