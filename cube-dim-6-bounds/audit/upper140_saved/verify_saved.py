#!/usr/bin/env python3
"""Finite independent C140 saved-candidate Arb768 verification; no exp/native."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,itertools,os,sys,time,traceback,importlib.util
sys.dont_write_bytecode=True
P=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-gaussian106/cube-dim-6-bounds/upper/round140/common_entropy_band_certificate')
A=Path(__file__).resolve().parent
G=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round140/primary_gate')
DEN=2**24
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(n,d):
 with (A/n).open('x')as f:json.dump(d,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def key(k,p):return hashlib.sha256((json.dumps({'semantic_key':k,'point':p},sort_keys=True,separators=(',',':'))+'\n').encode()).hexdigest()
def pts(box):return [['1']+[box[i][b]for i,b in enumerate(bits)]for bits in itertools.product((0,1),repeat=5)]
def vol(box):
 v=F(1)
 for a,b in box:assert F(a)<F(b);v*=F(b)-F(a)
 return v
r=read(P/'result.json');x=read(P/'exit_receipt.json');s=read(P/'started.json');release=read(P/'parent_release.json');terminal=read(G/'terminal_receipt.json')
assert sha(P/'run_batch.py')=='81588269d309deacd7cfa44fb2829a18ee1eb5db5083811fedfdfc17a547c91c'
assert r['status']==x['status']=='FINITE_COMMON_ENTROPY_BATCH_COMPLETE' and x['result_sha256']==sha(P/'result.json')
assert r['release_sha256']==s['release_sha256']==terminal['release_sha256']==sha(P/'parent_release.json')
assert terminal['exit_code']==0 and terminal['restarts']==0 and terminal['time_limit']is None and s['nice']>=10 and s['pid']==terminal['pid']
assert sha(G/'stdout.txt')==terminal['stdout_sha256'] and sha(G/'stderr.txt')==terminal['stderr_sha256']
spec=importlib.util.spec_from_file_location('c140_saved_runner',P/'run_batch.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);m.release_gate();sel,cache=m.preflight()
assert r['source_input_hashes_sha256']==sha(P/'input_hashes.json')
assert release['independent_review_report_sha256']==sha(A.parent/'upper140_preflight/REVIEW.md')
assert r['actual_counts']==x['actual_counts']
assert not r['global_g_certificate'] and not r['new_C6_bound']
leaves=r['leaf_results'];assert len(leaves)==11 and len({(l['root_id'],l['leaf_id'])for l in leaves})==11
assert [(l['root_id'],l['leaf_id'],l['box'])for l in leaves]==[(l['root_id'],l['leaf_id'],l['box'])for l in sel['selected_leaves']]
used=set();weights=set();paircount=0
for leaf,meta in zip(leaves,sel['selected_leaves']):
 rows=leaf['fixed_pair_results'];assert len(rows)==256 and {(p['band0_source_index'],p['band1_source_index'])for p in rows}==set(itertools.product(range(16),repeat=2))
 for row in rows:
  assert type(row['weight_numerator'])is int and 0<=row['weight_numerator']<=DEN and row['weight_denominator']==DEN
  assert F(row['score_min_lower'])<=F(row['score_min_upper']);assert row['all_corner_score_lowers_nonnegative']==(F(row['score_min_lower'])>=0)
  weights.add(row['weight_numerator']);paircount+=1
 b=leaf['selected_best_certificate'];best=max(rows,key=lambda p:F(p['score_min_lower']))
 assert all(b[k]==v for k,v in best.items()) and leaf['whole_leaf_pass']==(F(b['score_min_lower'])>=0)
 assert b['band0_semantic_key']==meta['source_keys'][0][b['band0_source_index']] and b['band1_semantic_key']==meta['source_keys'][1][b['band1_source_index']]
 assert len(b['corner_score_bounds'])==32 and min(F(v[0])for v in b['corner_score_bounds'])==F(b['score_min_lower']) and min(F(v[1])for v in b['corner_score_bounds'])==F(b['score_min_upper'])
 for keys in meta['source_keys']:
  for k in keys:
   for point in pts(meta['box']):
    ck=key(k,point);v=cache[ck];assert v['point']==point and v['semantic_key']==k and 0<F(v['lower'])<=F(v['upper']);used.add(ck)
assert paircount==2816 and len(used)==8544
newlogcap=len(used)+2*len(weights-{0,DEN});assert newlogcap==9190
assert r['actual_counts']=={'Arb_source_logs':8544,'Arb_entropy_logs':646,'Arb_corner_score_combinations':90112,'float_source_logs':8544,'float_entropy_logs':382976,'float_objective_calls':197120,'pair_certificates':2816,'native_calls':0,'heat_kernel_calls':0,'source_point_exp_calls':0,'final_exp_calls':0,'old_paid_leaf_retries':0}
# Independent exact recursive partition bridging; old paid status is referenced only.
S=Path(read(P/'input_hashes.json')['source_directory_allowlist']);cfg=read(S/'fixed_roots_and_atlas.json');bridge=[];oldpaid=0
for root in cfg['roots']:
 rid=root['original_node_id'];tree=read(S/'roots'/(rid+'.json'));nodes=tree['nodes'];seen=set();terminalids=[]
 def visit(nid,parent,depth):
  assert nid not in seen;seen.add(nid);n=nodes[nid];assert n['id']==nid and n['parent']==parent and n['depth']==depth
  if n['status']=='SPLIT_CLOSED_RECTANGLE':
   axis=max(range(5),key=lambda i:(F(n['box'][i][1])-F(n['box'][i][0]),-i));mid=sum(map(F,n['box'][axis]))/2
   expected=[[v.copy()for v in n['box']]for _ in range(2)];expected[0][axis][1]=str(mid);expected[1][axis][0]=str(mid);assert len(n['children'])==2
   for j,ch in enumerate(n['children']):assert nodes[ch]['box']==expected[j];visit(ch,nid,depth+1)
  else:assert n['status']in ('LOCAL_LEAF_PASS','DEPTH6_UNRESOLVED');terminalids.append(nid)
 visit('r',None,0)
 assert seen==set(nodes) and len(terminalids)==len(set(tree['closed_partition_leaf_ids'])) and set(terminalids)==set(tree['closed_partition_leaf_ids'])
 assert nodes['r']['box']==root['box']==tree['original_box'] and sum(vol(nodes[n]['box'])for n in terminalids)==vol(root['box'])
 old=[n for n in terminalids if nodes[n]['status']=='LOCAL_LEAF_PASS'];pending=[n for n in terminalids if nodes[n]['status']=='DEPTH6_UNRESOLVED'];oldpaid+=len(old)
 matches=[l for l in leaves if l['root_id']==rid];assert {l['leaf_id']for l in matches}==set(pending)
 passed=[l['leaf_id']for l in matches if l['whole_leaf_pass']];remaining=[n for n in pending if n not in passed]
 bridge.append({'root_id':rid,'old_paid':len(old),'new_passes':passed,'remaining_unresolved':remaining,'local_complete':not remaining})
assert oldpaid==62 and sum(b['local_complete']for b in bridge)==4 and sum(len(b['remaining_unresolved'])for b in bridge)==9
assert r['original_seven_local_box_completion']==[{'root_id':b['root_id'],'all_original_partition_leaves_paid':b['local_complete'],'old_paid_leaves_referenced_only':True}for b in bridge]
# Finite registration is written before importing Arb or doing new science.
files=[P/n for n in ['run_batch.py','PROOF.md','selected_metadata.json','input_hashes.json','prepared_packet.json','parent_release.json','started.json','result.json','exit_receipt.json']]+[G/'terminal_receipt.json']
registration={'status':'REGISTERED_BEFORE_NEW_LOGS','payload_hashes':{str(p):sha(p)for p in files},'precision_bits':768,'source_logs_cap':8544,'entropy_logs_cap':646,'total_logs_cap':9190,'corner_score_combinations_cap':90112,'fixed_pair_certificates':2816,'float_selector_repeated':False,'exp_calls_permitted':0,'native_calls_permitted':0,'heat_kernel_calls_permitted':0,'time_limit':None,'workers':1,'nice_minimum':10,'reviewer_script_sha256':sha(__file__),'scope':'saved C140 dyadic candidates only; negative results do not exclude other weights/sources'}
write('finite_registration.json',registration)
os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
assert all(os.environ.get(k)=='1'for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'])
import flint
from flint import arb,ctx
assert flint.__version__=='0.9.0' and flint.__FLINT_VERSION__=='3.6.0'
ctx.prec=768;logs={};entropies={};counts={'source_logs':0,'entropy_logs':0,'corner_score_combinations':0,'pair_certificates':0,'exp':0,'native':0,'heat_kernel':0};start=time.monotonic();cpu=time.process_time()
def ar(q):q=F(q);return arb(q.numerator)/q.denominator
def endpoints(v):assert v.is_finite();return F(str(v.lower().fmpq())),F(str(v.upper().fmpq()))
def lg(k,p):
 ck=key(k,p)
 if ck not in logs:
  counts['source_logs']+=1;assert counts['source_logs']<=8544;logs[ck]=ar(F(cache[ck]['lower'])*63500).log()
 return logs[ck]
def H(n):
 if n not in entropies:
  if n in (0,DEN):entropies[n]=arb(0)
  else:
   counts['entropy_logs']+=2;assert counts['entropy_logs']<=646;w=ar(F(n,DEN));z=ar(F(DEN-n,DEN));entropies[n]=-w*w.log()-z*z.log()
 return entropies[n]
write('started.json',{'registration_sha256':sha(A/'finite_registration.json'),'precision_bits':768,'nice':os.getpriority(os.PRIO_PROCESS,0),'pid':os.getpid(),'time_limit':None})
try:
 audits=[]
 for leaf,meta in zip(leaves,sel['selected_leaves']):
  vertices=pts(meta['box']);arrays=[[[lg(k,p)for p in vertices]for k in keys]for keys in meta['source_keys']];selected=leaf['selected_best_certificate'];corner_saved=selected['corner_score_bounds'];passpairs=0;negativepairs=0;best_corners=None
  for row in leaf['fixed_pair_results']:
   i=row['band0_source_index'];j=row['band1_source_index'];n=row['weight_numerator'];w=ar(F(n,DEN));z=ar(F(DEN-n,DEN));h=H(n)
   scores=[endpoints(h+w*a+z*b)for a,b in zip(arrays[0][i],arrays[1][j])];lo=min(v[0]for v in scores);hi=min(v[1]for v in scores)
   assert F(row['score_min_lower'])<=lo<=hi<=F(row['score_min_upper'])
   assert (lo>=0)==row['all_corner_score_lowers_nonnegative']
   if lo>=0:passpairs+=1
   if hi<0:negativepairs+=1
   counts['corner_score_combinations']+=32;counts['pair_certificates']+=1;assert counts['corner_score_combinations']<=90112 and counts['pair_certificates']<=2816
   if (i,j)==(selected['band0_source_index'],selected['band1_source_index']):
    assert n==selected['weight_numerator']
    for (l,u),(sl,su)in zip(scores,corner_saved):assert F(sl)<=l<=u<=F(su)
    assert (lo>=0)==leaf['whole_leaf_pass'];best_corners=[[str(l),str(u)]for l,u in scores]
  assert best_corners is not None
  audits.append({'root_id':leaf['root_id'],'leaf_id':leaf['leaf_id'],'whole_leaf_pass':leaf['whole_leaf_pass'],'saved_selected_fixed_weight':str(F(selected['weight_numerator'],DEN)),'selected_band0_key':selected['band0_semantic_key'],'selected_band1_key':selected['band1_semantic_key'],'selected_all_32_corner_bounds_Arb768':best_corners,'selected_candidate_pass_pair_count':passpairs,'selected_candidate_negative_upper_pair_count':negativepairs,'scope':'only saved dyadic candidates; no continuum/source-family exclusion'})
 assert counts['source_logs']==8544 and counts['entropy_logs']==646 and counts['corner_score_combinations']==90112 and counts['pair_certificates']==2816
 assert sum(l['whole_leaf_pass']for l in audits)==r['passed_leaf_count']==2
 for p in files:assert sha(p)==registration['payload_hashes'][str(p)]
 out={'verdict':'PASS_SAVED_RESULT_INDEPENDENT_ARB768','fixed_runner_sha256':sha(P/'run_batch.py'),'saved_result_sha256':sha(P/'result.json'),'registration_sha256':sha(A/'finite_registration.json'),'precision_bits':768,'actual_counts':counts,'leaf_audits':audits,'partition_bridging':bridge,'old_complete_boxes':3,'now_complete_boxes':4,'new_leaf_passes':2,'remaining_unresolved_boxes':3,'remaining_unresolved_leaves':9,'global_g_certificate':False,'new_C6_bound':False,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
 write('verification_result.json',out);write('exit_receipt.json',{'status':out['verdict'],'exit_code':0,'result_sha256':sha(A/'verification_result.json'),'actual_counts':counts,'wall_seconds':out['wall_seconds'],'CPU_seconds':out['CPU_seconds'],'time_limit':None});print(json.dumps({k:out[k]for k in ['verdict','actual_counts','old_complete_boxes','now_complete_boxes','new_leaf_passes','remaining_unresolved_boxes','remaining_unresolved_leaves','wall_seconds','CPU_seconds']}))
except BaseException:
 write('exit_receipt.json',{'status':'FAILED_STOP_NO_RETRY','failure':traceback.format_exc(),'actual_counts':counts,'time_limit':None});raise
