import sys,json,hashlib,copy,importlib.util
from pathlib import Path
from fractions import Fraction as Q
sys.dont_write_bytecode=True
p=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-numerics/cube-dim-6-bounds/upper/round129/recovered_continuation');sys.path.insert(0,str(p))
sp=importlib.util.spec_from_file_location('review129_controller',p/'run_prepared_129.py');R=importlib.util.module_from_spec(sp);sp.loader.exec_module(R)
U=R.U;H,B=R.model_modules();f=U.read(R.REC_FOREST);g=U.read(p/'continuation_seed129.json');u=U.read(R.REC_UNKNOWN);b=U.read(R.REC_BUDGET);c=U.read(p/'remaining_budgets129.json');r=U.read(p/'unknown_geometry_receipt.json');cache=U.read(R.REC_CACHE)
checks={};checks['325_roots_preserved']=len(f['roots'])==len(g['roots'])==325 and set(f['roots'])==set(g['roots'])
unknown_nodes={(x['owner_root'],x['owner_node_id']) for x in u};ukeys={x['semantic_key'] for x in u}
checks['15_distinct_unknowns']=len(u)==len(unknown_nodes)==len(ukeys)==15 and not (ukeys & cache['entries'].keys())
checks['15193_cache_retained']=len(cache['entries'])==15193 and cache['baseline_accepted_count']==12081 and cache['recovered_accepted_count']==3112
checks['276_closed_roots_preserved']=sum(s['root_complete'] for s in f['roots'].values())==sum(s['root_complete'] for s in g['roots'].values())==276
checks['closed_roots_byte_equivalent']=all(s==g['roots'][rid] for rid,s in f['roots'].items() if s['root_complete'])
checks['all_original_nodes_preserved']=all(n==g['roots'][rid]['nodes'][nid] for rid,s in f['roots'].items() for nid,n in s['nodes'].items() if (rid,nid)not in unknown_nodes)
checks['15_exact_original_split_make_node']=True;checks['child_semantic_keys_not_unknown']=True
for row in r['records']:
 rid,nid=row['root_id'],row['unknown_parent_node'];old=f['roots'][rid]['nodes'][nid];new=g['roots'][rid]['nodes'][nid];split=B.split_box(old['box']);kids=[B.make_node(nid+suffix,nid,old['depth']+1,raw,0) for suffix,raw in [('0',split['lower_raw_box']),('1',split['upper_raw_box'])]]
 checks['15_exact_original_split_make_node'] &= split==row['split'] and kids==row['children'] and all(k==g['roots'][rid]['nodes'][k['id']] for k in kids) and new['children']==[k['id'] for k in kids] and new['unknown_attempt_preserved'] in u
 for kid in kids:
  if kid['status']!='empty':
   for w in R.WINDOWS:
    e={**new['unknown_attempt_preserved']['expected'],'center':B.midpoint(kid['box']),'face_box':kid['box'],'window':w}
    checks['child_semantic_keys_not_unknown'] &= U.semantic_key(e,R.FAN_SHA)not in ukeys
checks['55_original_root_budget_no_reset']=len(b)==len(c)==55 and set(b)==set(c) and all(c[rid]['native_remaining']==512-x['native_used_127'] and c[rid]['node_remaining']==1024-x['processed_127']-c[rid]['recovery_geometry_debit'] for rid,x in b.items())
checks['3127_including15_charged']=sum(x['native_used_127'] for x in b.values())==3127
checks['remaining25033_and53198']=sum(x['native_remaining'] for x in c.values())==25033 and sum(x['node_remaining'] for x in c.values())==53198 and sum(x['recovery_geometry_debit'] for x in c.values())==15
checks['geometry_receipt_pins']=r['source_recovered_forest_sha256']==U.sha(R.REC_FOREST) and r['source_unknown_sha256']==U.sha(R.REC_UNKNOWN) and r['seed_sha256']==U.sha(p/'continuation_seed129.json')
checks['geometry_no_native_no_exp']=r['native_calls']==r['exp_calls']==0
out={'status':'PASS_SAVED_ONLY_RECOVERY_GEOMETRY' if all(checks.values()) else 'FAIL','checks':checks,'new_native_calls':0,'exponential_evaluations':0,'old_baseline_replays':0,'scope':'saved JSON/hash/graph/rational split verification only; no native or full saved certificate replay','pins':{n:U.sha(p/n)for n in ['continuation_seed129.json','unknown_geometry_receipt.json','remaining_budgets129.json']}}
Path(__file__).with_name('round129_saved_geometry_review.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
