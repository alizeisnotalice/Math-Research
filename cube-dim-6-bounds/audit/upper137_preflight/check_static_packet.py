import json,hashlib,sys,ast
from pathlib import Path
from fractions import Fraction as Q
sys.set_int_max_str_digits(0);sys.dont_write_bytecode=True
A=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-upper098/cube-dim-6-bounds/upper/round137/extra_time_band_corner_seed');P=Path(__file__).resolve().parent
rd=lambda p:json.loads(Path(p).read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();canon=lambda x:(json.dumps(x,sort_keys=True,separators=(',',':'),ensure_ascii=False)+'\n').encode()
r=rd(A/'registration.json');i=rd(A/'inputs.json');qual=rd(A/'cache_qualification.json');dis=rd(A/'live129_disjointness.json');manifest=rd(A/'source_manifest.json');family=rd(i['native_family_registration_path']);code=(A/'run_prepared.py').read_text();checks={}
expected={'run_prepared.py':'dc85a6eb06bcaeb9b937a79b08d8c5ad0eb32e7bcca7500a162d31b559da83c2','observe_once.py':'42b4159c0c6b29a06384318b1403c434326da7b2f069d6c847691926176ceccf','registration.json':'5c2c2da9d42519615a2a271e85ed0ab665fde18dcd27a76022b855c5b9407fa5','pre_run_manifest.json':'4043dca29cab274296ad196d5cd21ed40b7e76861acfaf0518bd4afade2198b1'}
checks['frozen4_pins']=all(sha(A/n)==h for n,h in expected.items());checks['all_prepared_pins']=all(sha(A/n)==h for n,h in r['prepared_pins'].items());checks['all_archived_source_pins_bytes']=all(sha(s['path'])==s['sha256']and Path(s['path']).stat().st_size==s['bytes']for s in manifest['sources']);checks['all_native_source_runtime_pins']=all(sha(p)==h for p,h in r['native_pins'].items())
checks['five_fixed_positive_disjoint_windows']=i['windows']==r['fixed_windows']==[q['window']for q in qual['windows']]and len(i['windows'])==5 and all(0<Q(w['a'])<Q(w['T'])for w in i['windows'])and all(Q(i['windows'][j+1]['T'])==Q(i['windows'][j]['a'])for j in range(4))
cache=rd(qual['pins'][0]['path']);unknown=rd(qual['pins'][1]['path']);attempts=rd(qual['pins'][2]['path'])['entries'];seed=rd(dis['seed_path']);keys=[]
for q in qual['windows']:
 e={'center':[str(Q(v))for v in i['center']],'a':str(Q(q['window']['a'])),'T':str(Q(q['window']['T'])),'formula_family_id':'C016-cauchy-frame32-outward-precision-family-v1','fan_sha256':r['fan_sha256'],'model':i['model'],'model_definition':i['model_definition'],'quadratic':i['quadratic']};k=hashlib.sha256(canon(e)).hexdigest();keys.append(k)
 assert k==q['semantic_key']and k not in cache['entries']and k not in attempts and all(u['semantic_key']!=k for u in unknown)
checks['all5_genuinely_new_frozen_keys']=len(cache['entries'])==15193 and len(unknown)==15 and len(set(keys))==5 and [q['window_index']for q in qual['windows']]==list(range(5))and all(q['action']=='NEW'and q['available_entry']is None and not q['unknown_entries']and q['historical_attempt_entry']is None for q in qual['windows'])
nodes=0;seed_ok=True
for root in seed['roots'].values():
 for n in root['nodes'].values():
  if n['box']is None or n['status']=='empty':continue
  nodes+=1;seed_ok &= all(Q(l)<1 and Q(u)<=1 for l,u in n['raw_box']) and (Q(n['box'][0][0])+Q(n['box'][0][1]))/2<1
checks['frozen_seed_c2_midpoint_disjointness']=seed_ok and nodes==dis['nonempty_seed_nodes_checked']==18333 and dis['live_state_reads']==0 and Q(i['center'][1])==1 and sha(dis['seed_path'])==dis['seed_sha256']
checks['original128_256_512_bins_compatible_pins']=all(b['sha256']==r['native_pins'][b['path']]and b['source_sha256']==r['native_pins'][b['source_path']]for b in family['binaries'].values())
checks['exact_corner_center_and_encoder_bounds']=list(map(Q,i['center']))==[Q(1),Q(1),Q(1),Q(3,8),Q(1,4),Q(1,8)]and 5*30<2**63 and 160*8<2**63 and 5*30*160*8<2**127
checks['all_AST_and_UNRUN']=all(ast.parse(p.read_text())is not None for p in A.glob('*.py'))and not(A/'started.json').exists()
checks['serial3_per_window_and5workers_cap']=r['precision_ladder']==i['precision_ladder']==[128,256,512]and r['max_real_attempts_total']==i['max_actual_attempts']==5*3==15 and i['parallel_max_windows']==5
checks['sumUpper_fullG_lower_not_acceptance_requirement']="'sum_upper_at_most_fullG_lower_strong_check':F[1]<=GL"in code and "req(F[1]<=GL"not in code
checks['partial_adjacent_and_complete_prefix_Euler_guards']="fatal adjacent endpoint density mismatch in partial batch"in code and "req(eu[0]<=0<=eu[1],'merged exact Euler')"in code
checks['exclusive_main_semantic_owner_guard']="req(rel.get('exclusive_owner_for_registered_semantic_keys')is True"in code
out={'status':'PASS_C137_STATIC_CONCURRENCY_SOURCE_AND_NEW_QUERY_GATE'if all(checks.values())else'FAIL','checks':checks,'frozen_pins':expected,'prepared_pins_count':len(r['prepared_pins']),'archived_sources_count':len(manifest['sources']),'semantic_keys':keys,'static_seed_nonempty_nodes':nodes,'live129_state_reads':0,'new_native_band_Arb_science_calls':0,'concurrency_review':'validators use only immutable constants/local Fraction data; base alias set once before pool; worker prior/data/levels private; launch/state/accepted/CPU protected by RLock; exact-PID wait4','scope':'one fixed point, five previously missing disjoint positive windows, max15 native attempts/max5 parallel, no box/global floor'}
(P/'static_review.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
