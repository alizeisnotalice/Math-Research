"""True2350 frontier: new single-source pairs and backup jet points, no cert/H."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import gzip,importlib.util,json,os,time,resource
import numpy as np
from scipy.spatial import cKDTree
from common import sha,digest,write
import residual_run
H=Path(__file__).resolve().parent;S=residual_run.S;W=S/'upper_parallel137/next2350_preflight';R=S/'upper_parallel137/residual_preflight';OLD=S/'upper_crossanchor131'
def pk(point):return ','.join(map(str,map(F,point)))
def emit(f,x):f.write(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n')
def main():
 began=time.monotonic();assert not W.exists();W.mkdir();manifest=H/'residual_results_manifest.json';assert sha(manifest)=='86a877cab4efcd4f8b7c3d5a22196d6593bb81c01c0d4179129d0077810e7ad9';frozen=json.load(open(manifest));package=json.load(open(H/'residual_prepared_manifest.json'))
 for p,h in frozen['source_locks'].items():assert sha(p)==h
 statepath=R/'future_batch0002/state.json';assert sha(statepath)=='52e6fe15b99fe228501a2b2de41d7c23c4c41f67d1790d088fa3c1ddd2fe9385';state=json.load(open(statepath));assert state['status']=='finite_original_prefix_complete' and not state['faults'] and not state['inflight'] and not state['obstructions_H128'] and len(state['queue'])==2350
 write(W/'registration.json',{'pid':os.getpid(),'code_sha256':sha(__file__),'source_manifest':str(manifest),'source_manifest_sha256':sha(manifest),'final_state':str(statepath),'final_state_sha256':sha(statepath),'target':'11/100000','epsilon_for_new_points':'1/128','children':2350,'maximum_pairs':75200,'selection':'direct parent first +nearest15 real anchors +nearest16 parabolic centers;whole-source exactdedup;max32','source_pool_union':'old131 immutable pool plus only4274 new137 jet sources;no full oldpool copy','future_worker_count':15,'actual_preparation_workers':1,'H_gradient_kernel_log_exp_corner_calls':0,'split_calls':0,'dispatch_authorized':False})
 sources=list(map(json.loads,open(OLD/'source_pool.jsonl')));assert len(sources)==111331;sourceby={s['source_id']:s for s in sources};proofs={r['certificate']['point_key']:(r['node'],r['certificate']) for r in state['new_passed']}
 for row in state['new_splits']:assert row['failed_certificate']['point_key'] not in proofs;proofs[row['failed_certificate']['point_key']]=(row['node'],row['failed_certificate'])
 assert len(proofs)==len(state['new_cache_by_epsilon']['1/128'])==4274;new=[];parent_source={}
 for pointkey,record in state['new_cache_by_epsilon']['1/128'].items():
  node,c=proofs[pointkey];assert record['epsilon']=='1/128' and record['qualification_target']==c['target']=='11/100000' and record['pair_cache_key']==['1/128',pointkey] and F(record['values'][0][0])>0
  core={'epsilon':'1/128','point':list(map(str,map(F,record['point']))),'target':'11/100000','Delta_lower':str(F(c['Delta_lower'])),'loggradient_directed':[[str(F(a)),str(F(b))] for a,b in c['loggradient_directed']],'penalty_coefficient':str(F(c['penalty_coefficient'])),'curvature_rule':'129_anchor_Jensen_g105_U128'};sid=digest(core);assert sid not in sourceby;assert F(core['penalty_coefficient'])==1024-768*F(21,200000)/F(record['values'][0][1])
  point=list(map(F,core['point']));gm=[(F(a)+F(b))/2 for a,b in core['loggradient_directed']];q=F(core['penalty_coefficient']);source={'source_id':sid,'core':core,'H_interval':record['values'][0],'raw_record_sha256':digest(record),'aliases':[{'origin':'137_new4274','source_node':node,'source_certificate_sha256':digest(c),'source_record_sha256':digest(record),'source_original_penalty':c['penalty_coefficient']}],'H64_g105_ancestry':{'node':c['H64_lower_source_node'],'H64_component_epsilon':'1/64','H64_lower':'21/200000','point_inside_ancestor':True,'H128_upper':c['H128_upper'],'proof':c['H64_g105_ancestor_proof'],'source_final137_state_sha256':sha(statepath)},'anchor_tail_float':list(map(float,point[1:])),'parabolic_center_float':[float(point[j+1]+gm[j]/(2*q)) for j in range(4)],'negative_Delta':F(core['Delta_lower'])<0,'H_lower_positive':True,'qualification':'production candidate;new4274 PF/Delta/corner independent audit pending'};new.append(source);parent_source[tuple(node)]=sid
 new.sort(key=lambda s:s['source_id'])
 with open(W/'new_source_pool.jsonl','x') as f:
  for i,source in enumerate(new,111331):source['source_index']=i;sourceby[source['source_id']]=source;emit(f,source)
 sources.extend(new);assert len(sources)==len(sourceby)==115605
 parents={tuple(r['node']):r for r in state['new_splits']};assert len(parents)==1175;rows=[]
 for node in state['queue']:
  parent=parents[(node[0],node[1][:-1])];assert node in parent['children'];lo=[F(x,32) for x in package.get('unused',[])];rows.append({'node':node,'direct_parent_node':parent['node'],'parent_source_id':parent_source[tuple(parent['node'])],'parent_failed_margin':parent['failed_certificate']['minimum_log_margin'],'physical_volume':str(F(1,32**4*2**len(node[1])))})
 # Root tuples are deterministic ordered roots, avoiding a complete cache load.
 roots=[[a,b,c,d] for a in range(32) for b in range(a+1) for c in range(b+1) for d in range(c+1)]
 def box(node):
  lo=[F(x,32) for x in roots[node[0]]];hi=[a+F(1,32) for a in lo]
  for j,b in enumerate(node[1]):
   k=j%4;m=(lo[k]+hi[k])/2
   if b=='0':hi[k]=m
   else:assert b=='1';lo[k]=m
  return lo,hi
 rows.sort(key=lambda r:(F(r['parent_failed_margin']),-F(r['physical_volume']),*r['node']));ta=cKDTree(np.asarray([s['anchor_tail_float'] for s in sources]));tc=cKDTree(np.asarray([s['parabolic_center_float'] for s in sources]));known={tuple(r) for r in json.load(open(OLD/'known_failure_objects.json'))['source_child_curvature_failures']}
 for p in (S/'upper_parallel137/run001').glob('shard*/results.jsonl'):
  for r in map(json.loads,open(p)):
   if r['kind']=='evaluated' and not r['certificate']['passed']:known.add((r['source_id'],*r['certificate']['child_node'],r['certificate']['curvature_rule']))
 for node,c in proofs.values():
  if not c['passed']:known.add((parent_source[tuple(node)],*node,'129_anchor_Jensen_g105_U128'))
 pair_count=excluded=dedup=0;selected_ids=set();children=[]
 with open(W/'pairs.jsonl','x') as f:
  for child_index,row in enumerate(rows):
   lo,hi=box(row['node']);mid=np.asarray([float((a+b)/2) for a,b in zip(lo,hi)]);ia=ta.query(mid,k=15,workers=1)[1];ic=tc.query(mid,k=16,workers=1)[1];raw=[row['parent_source_id']]+[sources[int(i)]['source_id'] for i in ia]+[sources[int(i)]['source_id'] for i in ic];chosen=[];seen=set()
   for sid in raw:
    source=sourceby[sid];obj=(sid,*row['node'],source['core']['curvature_rule'])
    if obj in known:excluded+=1;continue
    if sid in seen:dedup+=1;continue
    seen.add(sid);chosen.append(sid);selected_ids.add(sid);emit(f,{'pair_index':pair_count,'child_index':child_index,'child_node':row['node'],'source_id':sid,'source_index':source['source_index'],'source_pool':'old131' if source['source_index']<111331 else 'new1374274','curvature_rule':source['core']['curvature_rule']});pair_count+=1
   assert chosen and len(chosen)<=32 and chosen[0]==row['parent_source_id'];children.append({**row,'child_index':child_index,'box_lower':list(map(str,lo)),'box_upper':list(map(str,hi)),'selected_sources':chosen,'pair_count':len(chosen)})
 write(W/'children.json',children);write(W/'source_pool_refs.json',{'old_pool':{'path':str(OLD/'source_pool.jsonl'),'sha256':sha(OLD/'source_pool.jsonl'),'objects':111331},'new_pool':{'path':str(W/'new_source_pool.jsonl'),'sha256':sha(W/'new_source_pool.jsonl'),'objects':4274},'whole_library':115605,'no_old_pool_copy':True,'new_source_independent_numeric_qualification':'pending audit140/primary;preparation only'})
 # Concrete backup new-jet points, not a dispatch plan and not scheduled with pairs.
 pre=json.load(open(R/'registration.json'));planned={}
 for p in pre['source_locks']:
  if p.endswith('jobs.json.gz'):
   planned[p]={pk(r['point']) for r in json.loads(gzip.decompress(Path(p).read_bytes()))}
 for p in [S/'upper_g110107/planned_requests.json',S/'upper_continue123/planned_tasks.json',R/'tasks.json']:planned[str(p)]={pk(t['point']) for t in json.load(open(p))}
 full={pk(s['core']['point']):s for s in sources if s['core']['epsilon']=='1/128'};assert len(full)==14801;spec=importlib.util.spec_from_file_location('readonly128_next_geometry',package['source_locks'] and H.parent/'upper_mixed098/epsilon_kernel.py');api=importlib.util.module_from_spec(spec);spec.loader.exec_module(api);jet=api.FirstJetKernel(json.loads(gzip.decompress(Path(package['table_path']).read_bytes())),F(1,128));points=[];used=set();blocked=[];counts=Counter()
 for row in children:
  lo=list(map(F,row['box_lower']));hi=list(map(F,row['box_upper']));mid=[(a+b)/2 for a,b in zip(lo,hi)];eta=min(b-a for a,b in zip(lo,hi))/4096;candidates=[[F(1),*mid]]+[[F(1),*[x+eta*beta**j for j,x in enumerate(mid)]] for beta in (F(1,8),F(1,16),F(1,32))];tests=[];choice=None
  for j,point in enumerate(candidates):
   p=pk(point);g=jet.eligible(point);cache=p in full;collision=[n for n,keys in planned.items() if p in keys];failed=cache and (full[p]['source_id'],*row['node'],full[p]['core']['curvature_rule']) in known;assert all(a<x<b for a,x,b in zip(lo,point[1:],hi));reason='generic280' if not g['eligible'] else 'already_failed_same_region_point' if failed else 'duplicate_newpoint' if p in used else 'scalar_or_planned_object' if collision and not cache else None;tests.append({'candidate':j,'point':list(map(str,point)),'geometry':g,'completed_fulljet':cache,'registered_collisions':collision,'rejection':reason})
   if reason is None:choice=(j,point,p,cache);break
  if choice is None:blocked.append({'node':row['node'],'tests':tests});continue
  j,point,p,cache=choice;used.add(p);counts['cache_hit' if cache else 'fresh_H128_if_later_authorized']+=1;counts['candidate'+str(j)]+=1;points.append({'node':row['node'],'child_index':row['child_index'],'epsilon':'1/128','target':'11/100000','point':list(map(str,point)),'pair_cache_key':['1/128',p],'candidate':j,'candidate_tests':tests,'cached_source_id':full[p]['source_id'] if cache else None,'same_point_H64_guard':None,'H64_g105_ancestor':sourceby[row['parent_source_id']]['H64_g105_ancestry'],'actual_request':False})
 assert len(points)+len(blocked)==2350;write(W/'backup_newpoint_candidates.json',points);write(W/'blocked_points.json',blocked)
 summary={'status':'true2350_new_pairs_and_backup_points_prepared_NOT_executed','pid':os.getpid(),'children':2350,'pairs':pair_count,'maximum_pairs':75200,'source_library':115605,'new_137_sources':4274,'selected_distinct_sources':len(selected_ids),'same_source_child_rule_excluded':excluded,'duplicate_candidate_hits_removed':dedup,'pair_count_histogram':dict(Counter(r['pair_count'] for r in children)),'backup_point_counts':dict(counts),'blocked_points':len(blocked),'actual_H_gradient_kernel_log_exp_corner_calls':0,'new_splits':0,'future_workers':15,'production_or_corner_authorized':False,'wall_seconds':time.monotonic()-began,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};write(W/'summary.json',summary);print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
