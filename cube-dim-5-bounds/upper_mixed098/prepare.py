#!/usr/bin/env python3
"""Prepare a mixed-epsilon G state and finite194 H128 input; no H calls."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from datetime import datetime,timezone
import argparse,hashlib,json,os,time,resource
import epsilon_kernel as kernel
from cache_model import lookup,pair_key
ROOT=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/cube-dim-5-bounds/upper_mixed098')
BASE=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_g105086/final/state.json')
BASE_SHA='8647a4543ef6a869e747ec06d93f5137cf11feb7b3a64a002995d25a681c4f9a'
STOP=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_g105089/batch0009/state.json')
STOP_SHA='ba9ba3ac4f8e7b6788ed75005643dae86ed4e2f139c3a11bf4f12a986b491e74'
TERM=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/cube-dim-5-bounds/upper_g105089/terminal_manifest.json')
TERM_SHA='c8b639b5ff46b36a602c1f38f2fc291d76a7b87a9b76218a6c808b9f875012c7'
TABLE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds/astra_vertex017/boundary_triangulation.json.gz')
TABLE_SHA='8dd2e9895d7465ac599d43feb1b5099506f4e43d1248215eba0a77e11f4d13f9'
ENC=json.JSONEncoder(sort_keys=True,separators=(',',':'))
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
def digest(obj):
 h=hashlib.sha256()
 for x in ENC.iterencode(obj):h.update(x.encode())
 return h.hexdigest()
def write(p,obj):
 temp=p.with_suffix(p.suffix+'.tmp')
 with open(temp,'w') as f:
  for x in ENC.iterencode(obj):f.write(x)
  f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(temp,p)
def box(node,roots):
 lo=[F(x,32) for x in roots[node[0]]];hi=[x+F(1,32) for x in lo]
 for bit in node[1]:
  j=max(range(4),key=lambda j:hi[j]-lo[j]);m=(lo[j]+hi[j])/2
  if bit=='0':hi[j]=m
  else:assert bit=='1';lo[j]=m
 return lo,hi
def main():
 import gzip
 began=time.monotonic();ap=argparse.ArgumentParser();ap.add_argument('--work-dir',type=Path,required=True);args=ap.parse_args();out=args.work_dir.resolve();out.mkdir(parents=True,exist_ok=True)
 assert sha(BASE)==BASE_SHA and sha(STOP)==STOP_SHA and sha(TERM)==TERM_SHA and sha(TABLE)==TABLE_SHA
 reg={'task':'upper_mixed098_preparation','pid':os.getpid(),'registered_utc':datetime.now(timezone.utc).isoformat(),
  'source086':str(BASE),'source086_sha256':BASE_SHA,'source089_stop':str(STOP),'source089_stop_sha256':STOP_SHA,
  'source089_terminal_manifest_sha256':TERM_SHA,'preparation_code_sha256':sha(Path(__file__)),'new_H128_kernel_sha256':sha(ROOT/'epsilon_kernel.py'),
  'cache_model_code_sha256':sha(ROOT/'cache_model.py'),'kernel_diff_sha256':sha(ROOT/'kernel.diff'),
  'object':'same full G; each closed leaf G>=H_epsilon with its own epsilon; homogeneity applies to G only',
  'target':'21/200000','accepted_H64_leaves':100804,'pending_from_queue':190,'pending_from_counter_inflight':4,'new_pending_epsilon':'1/128',
  'cache_key':'(canonical epsilon,canonical point); separate outer dictionaries, no flat-cache fallback',
  'new_H128_points_requested':194,'maximum_new_anchor_starts':194,'workers_proposed':12,'wall_deadline':None,'production_authorized':False,
  'new_H_physical_gradient_kernel_log_exp_Gamma_calls':0}
 write(out/'preparation_registration.json',reg);base=json.load(open(BASE));old=json.load(open(STOP))
 assert old['status']=='counterexample' and len(old['queue'])==190 and len(old['inflight'])==4 and len(base['passed'])+len(old['new_passed'])==100804
 assert base['epsilon']==base['cache_model_epsilon']=='1/64' and base['target']=='21/200000'
 old64=dict(base['cache']);assert len(old64)==180040
 for a,v in old['new_cache'].items():assert a not in old64;old64[a]=v
 assert len(old64)==184792
 inflight_records=[]
 for ordinal,entry in sorted(old['inflight'].items(),key=lambda x:int(x[0])):
  got=entry['returned'];assert got['kind']=='counterexample';record=got['record'];assert F(record['values'][0][1])<F(21,200000)
  a=kernel.pointkey(record['point']);assert a not in old64;old64[a]=record
  inflight_records.append({'node':entry['node'],'source089_inflight_ordinal':ordinal,'H64_record_sha256':digest(record),'H64_point':record['point'],'H64_interval':record['values'][0]})
 assert len(old64)==184796
 accepted=[];acceptedmap=[]
 for origin,leaves in [('source086_passed',base['passed']),('source089_new_passed',old['new_passed'])]:
  for i,leaf in enumerate(leaves):
   c=leaf['certificate'];assert c['passed'] is True and F(c['minimum_log_margin'])>0
   if 'epsilon' in c:assert c['epsilon']=='1/64'
   if 'target' in c:assert c['target']=='21/200000'
   accepted.append({**leaf,'leaf_epsilon':'1/64','lower_bound_object':'G','component_lower_function':'H_epsilon','target':'21/200000'})
   acceptedmap.append({'new_passed_index':len(accepted)-1,'node':leaf['node'],'source':origin,'source_index':i,'source_leaf_sha256':digest(leaf),'certificate_sha256':digest(c),'leaf_epsilon':'1/64'})
 pending=[]
 for i,node in enumerate(old['queue']):pending.append({'node':node,'source089_queue_index':i,'source_kind':'queue','leaf_epsilon_for_next_evaluation':'1/128'})
 for r in inflight_records:pending.append({**r,'source_kind':'retained_counter_inflight','leaf_epsilon_for_next_evaluation':'1/128'})
 pending.sort(key=lambda x:(0 if x['source_kind']=='retained_counter_inflight' else 1,
  F(x['H64_interval'][1]) if x['source_kind']=='retained_counter_inflight' else F(old['hardness'][str(x['node'][0])+':'+x['node'][1]][0]),*x['node']))
 assert len(pending)==194 and len({tuple(x['node']) for x in pending})==194
 for i,r in enumerate(pending):r['new_queue_index']=i
 queue=[r['node'] for r in pending];splits=base['splits']+old['new_splits'];assert len(splits)==48638
 mixed={'task':'upper_mixed098','status':'prepared_partial_G105_no_new_H128','model':'alpha1,16 unnormalised body Cauchy directions,(2pi)^-5',
  'lower_bound_object':'G','target':'21/200000','accepted_leaf_epsilon':'per leaf,currently all1/64','pending_evaluation_epsilon':'1/128',
  'roots':base['roots'],'splits':splits,'passed':accepted,'queue':queue,'inflight':{},'resume_jobs':{},
  'cache_key_format':'(canonical epsilon,canonical point)','cache_by_epsilon':{'1/64':old64,'1/128':{}},
  'source089_fixed_H64_obstruction':{'terminal_manifest':str(TERM),'terminal_manifest_sha256':TERM_SHA,'stop_overlay_sha256':STOP_SHA,'scope':'H64 target only;not a G refutation'},
  'source_refs':{'source086':str(BASE),'source086_sha256':BASE_SHA,'source089_stop':str(STOP),'source089_stop_sha256':STOP_SHA},
  'certificate_rule':'each accepted leaf uses its stored epsilon for H bound and Q=I4/(16epsilon^2);same full G>=H_epsilon',
  'global_extension':'G homogeneity only;do not assert a mixed-H same-epsilon scaling law','g105_closed_cover':False,'new_H128_records':0,'production_authorized':False}
 for r in pending:
  r['box_lower'],r['box_upper']=(list(map(str,a)) for a in box(r['node'],mixed['roots']))
 paths=defaultdict(list);internals=set()
 for leaf in accepted:paths[leaf['node'][0]].append(leaf['node'][1])
 for r,p in queue:paths[r].append(p)
 for split in splits:
  r,p=split['node'];assert (r,p) not in internals;internals.add((r,p));assert split['children']==[[r,p+'0'],[r,p+'1']] and split['axis']==len(p)%4
 for r in range(52360):
  ps=sorted(paths[r]);assert ps and all(not b.startswith(a) for a,b in zip(ps,ps[1:]));d=max(map(len,ps));assert sum(1<<(d-len(p)) for p in ps)==1<<d
  for p in ps:
   for j in range(len(p)):assert (r,p[:j]) in internals
 assert len(accepted)+len(queue)==52360+len(splits)==100998
 write(out/'accepted_H64_mapper.json',acceptedmap);write(out/'pending_transfer_mapper.json',pending)
 mixed['fingerprint098']=digest({'preparation_registration_sha256':sha(out/'preparation_registration.json'),'accepted_mapper_sha256':sha(out/'accepted_H64_mapper.json'),'pending_mapper_sha256':sha(out/'pending_transfer_mapper.json'),'source086_sha256':BASE_SHA,'source089_stop_sha256':STOP_SHA})
 write(out/'initial_mixed_state.json',mixed)
 # Integer-only generic and epsilon-key preflight, with zero evaluation.
 jet=kernel.FirstJetKernel(json.loads(gzip.decompress(TABLE.read_bytes())),F(1,128));requests=[];old64hits=0
 for i,node in enumerate(queue):
  lo,hi=box(node,mixed['roots']);mid=[(a+b)/2 for a,b in zip(lo,hi)];eta=min(b-a for a,b in zip(lo,hi))/4096
  candidates=[mid]+[[x+eta*beta**j for j,x in enumerate(mid)] for beta in (F(1,8),F(1,16),F(1,32))];tests=[];point=None
  for j,tail in enumerate(candidates):
   p=[F(1),*tail];test={'candidate':j,**jet.eligible(p)};tests.append(test)
   if test['eligible']:point=p;candidate=j;break
  assert point is not None and lookup(mixed,F(1,128),point) is None
  hit64=lookup(mixed,F(1,64),point) is not None;old64hits+=hit64
  requests.append({'ordinal':i,'node':node,'epsilon':'1/128','point':list(map(str,point)),'pair_cache_key':list(pair_key(F(1,128),point)),
   'candidate':candidate,'geometry_tests':tests,'H128_cache_hit':False,'H64_same_point_present':hit64,'target':'21/200000'})
 assert len({tuple(r['pair_cache_key']) for r in requests})==194 and not mixed['cache_by_epsilon']['1/128']
 write(out/'planned_H128_requests.json',requests)
 write(out/'first_batch_prepared_registration.json',{'task':'upper_mixed098_first_H128_batch','source_mixed_state_sha256':sha(out/'initial_mixed_state.json'),
  'selected_boxes':194,'maximum_new_anchor_starts':194,'epsilon':'1/128','target':'21/200000','workers':12,'wall_deadline':None,
  'planned_requests_sha256':sha(out/'planned_H128_requests.json'),'H128_cache_hits':0,'same_point_H64_present_but_not_reused':old64hits,
  'new_kernel_sha256':sha(ROOT/'epsilon_kernel.py'),'Q_penalty_coefficient':'1024','production_authorized':False,
  'stop':'finite194;strict H128 target obstacle;actual numerical/process/state fault;explicit STOP;all poolsclose/join'})
 assert sha(BASE)==BASE_SHA and sha(STOP)==STOP_SHA and sha(TERM)==TERM_SHA
 result={'status':'prepared_only_no_H128','pid':os.getpid(),'initial_mixed_state_sha256':sha(out/'initial_mixed_state.json'),'state_bytes':(out/'initial_mixed_state.json').stat().st_size,
  'accepted_H64_leaves':100804,'pending_H128_boxes':194,'requeued_counter_nodes':4,'old_queue_nodes_preserved':190,'source_closed_root_Kraft':52360,
  'H64_cache_entries':184796,'H128_cache_entries':0,'initial_H128_cache_hits':0,'same_point_H64_present_not_reused':old64hits,
  'new_H_physical_gradient_kernel_log_exp_Gamma_calls':0,'production_authorized':False,'wall_seconds':time.monotonic()-began,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
  'initial_source_refs_exact':True,'g105_closed_cover':False,'not_a_G_or_C5_refutation':True}
 write(out/'preparation_summary.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':main()
