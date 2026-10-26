"""Prepare153 curvature + g115 inputs only. No new q/log/corner evaluation."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter,defaultdict
import ast,gzip,hashlib,json,os,time,resource
D=Path('/Users/zhengzhihao/Desktop/T');S=D/'cube-finite-dim-scratch';H=Path(__file__).resolve().parent;W=S/'upper_retarget155';A=S/'upper_retarget149';BASE=H.parent
P153=D/'Math-Research-cube5-curvature128-sol/cube-dim-5-bounds/upper_curvature153';P154=D/'Math-Research-cube5-auditcurvature129-sol/cube-dim-5-bounds/audit_curvature154';ACCEPT=D/'Math-Research-cube-dim-5/cube-dim-5-bounds/round017/primary_global_acceptance.json'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def write(p,x):
 assert not Path(p).exists();Path(p).write_text(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n')
def line(f,x):f.write(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n')
ROOTS=[[a,b,c,d] for a in range(32) for b in range(a+1) for c in range(b+1) for d in range(c+1)]
def box(n):
 lo=[F(x,32) for x in ROOTS[n[0]]];hi=[x+F(1,32) for x in lo]
 for j,b in enumerate(n[1]):
  k=j%4;m=(lo[k]+hi[k])/2
  if b=='0':hi[k]=m
  else:assert b=='1';lo[k]=m
 return lo,hi
def main():
 began=time.monotonic();W.mkdir(exist_ok=True);assert not (W/'registration.json').exists();m149=BASE/'upper_retarget149/prepared_manifest.json';m151=BASE/'upper_global151/manifest.json';assert sha(m149)=='bcbd1f26218bc63018170463fac1b69ffbdad0bda003632249cf678ad0434f32' and sha(m151)=='2ad3625ea95d29447850e995465afc2480d4b5d61574180af23a080a4b0c98a8';assert sha(ACCEPT)=='41adc085a86bc462bcfd85404c1e2112c2d76781e2665a5b0633616b5539a192';assert sha(P153/'penalty.py')=='0f271bf3ffc6a52bb276d16a1144f4b0e748653054ce4d0690884d3219a13a4e'
 locks={str(m149):sha(m149),str(m151):sha(m151),str(ACCEPT):sha(ACCEPT)}
 for p in [P153/'penalty.py',P153/'manifest.json',P153/'report.md',P154/'manifest.json',P154/'report.md',BASE/'upper_retarget149/shared_log.py']:locks[str(p)]=sha(p)
 a=json.load(open(m149));p151=json.load(open(m151))
 for name in ['canonical_objects.jsonl','leaf_mapper.jsonl','source_provenance.jsonl']:
  path=A/name;assert sha(path)==a['prepared_external_files'][str(path)];locks[str(path)]=sha(path)
 for path,h in locks.items():assert sha(path)==h
 write(W/'registration.json',{'pid':os.getpid(),'code_sha256':sha(__file__),'source_locks':locks,'source_objects':107530,'source_leaves':111201,'old_target':'11/100000','new_target':'23/200000','theory':'153+154 fixedanchor samephysicalEuclidean q bounds','mass_policy':'only original qualified H64g105 for old129 EPS128;no new rawL or globalH128 gstar as H64','new_q_log_corner_H_calls':0,'worker_count_proposed':15,'source_ownership':'canonical_index modulo15;all aliases stay together;newq once per uniqueobject when authorized','wall_deadline':None})
 supports={}
 for part in p151['datasets']['supports']:
  path=BASE/'upper_global151'/part['path'];assert sha(path)==part['sha256'];locks[str(path)]=part['sha256']
  with gzip.open(path,'rt') as f:
   for r in map(json.loads,f):supports[tuple(r['node'])]=r
 assert len(supports)==3812;aliases={}
 for r in map(json.loads,open(A/'source_provenance.jsonl')):
  oid=r['canonical_id']
  if oid not in aliases or aliases[oid]['curvature_proof'] is None and r['curvature_proof'] is not None:aliases[oid]=r
 per_sources=[[] for _ in range(15)];per_leaves=[[] for _ in range(15)];owners={};classes=Counter();support_used=set();oldq_checks=0
 for obj in map(json.loads,open(A/'canonical_objects.jsonl')):
  oid=obj['canonical_id'];c=obj['core'];assert digest(c)==oid and F(c['target110'])==F(11,100000);e=F(c['epsilon']);oldq=F(c['penalty_coefficient']);upper=F(obj['H_interval'][1]);assert upper>0 and oldq>0;alias=aliases[oid];mass=None
  if e==F(1,64):assert c['curvature_rule_class']=='075_full_frame' and oldq==256;kind='153_EPS64_uniform_no_H32_mass'
  elif e==F(1,128) and c['curvature_rule_class']=='075_full_frame':assert oldq==1024;kind='153_EPS128_uniform_no_qualified_H64_mass_used'
  else:
   assert e==F(1,128) and c['curvature_rule_class']=='129_anchor_fixed_Jensen';assert F(21,200000)<upper and oldq==1024-768*F(21,200000)/upper;oldq_checks+=1;proof=alias['curvature_proof'];assert proof is not None;node=tuple(proof['node']);support=supports[node];cert=support['certificate'];assert support['typed_cache_key'][0]=='1/64' and F(cert['target'])==F(21,200000) and cert['passed'] and cert['checked_points']==16;lo,hi=box(node);assert all(a<=F(x)<=b for a,x,b in zip(lo,c['point'][1:],hi));support_used.add(node);kind='153_EPS128_min_old_uniform_split_using_only_old_H64g105';mass={'lower':'21/200000','lower_epsilon':'1/64','upper_epsilon':'1/128','upper':obj['H_interval'][1],'point':c['point'],'support_node':list(node),'support_source_leaf_sha256':support['source_leaf_sha256'],'support_source_state_sha256':support['source_state_sha256'],'support_certificate_sha256':digest(cert),'original_support_certificate_sha256':support['original_certificate_sha256'],'source_alias_index':alias['source_alias_index'],'source_alias_sha256':digest(alias),'same_point_same_model':True,'point_in_support_closedbox':True,'no_new_raw_lower_lookup':True,'global_H128_floor_not_used_as_H64':True}
  owner=obj['canonical_index']%15;assert oid not in owners;owners[oid]=owner;classes[kind]+=1;per_sources[owner].append({'canonical_id':oid,'canonical_index':obj['canonical_index'],'old_core':c,'old_H_interval':obj['H_interval'],'typed_cache_key':obj['typed_cache_key'],'raw_record_sha256':obj['raw_record_sha256'],'new_rule':kind,'new_q_status':'NOT_COMPUTED','mass_qualification':mass,'old_proof_alias_index':alias['source_alias_index'],'old_proof_alias_sha256':digest(alias)})
 assert len(owners)==107530;nodes=set();byroot=defaultdict(list);oldmin=None
 for r in map(json.loads,open(A/'leaf_mapper.jsonl')):
  node=tuple(r['node']);assert node not in nodes;nodes.add(node);byroot[node[0]].append(node[1]);margin=F(r['saved_g110_margin']);assert margin>0;oldmin=margin if oldmin is None else min(oldmin,margin);lo,hi=box(node);assert list(map(str,lo))==r['box_lower'] and list(map(str,hi))==r['box_upper'];per_leaves[owners[r['canonical_id']]].append({k:r[k] for k in ('leaf_index','node','canonical_id','source_alias_index','epsilon','box_lower','box_upper','closed_box','saved_g110_margin')})
 assert len(nodes)==111201 and len(byroot)==52360
 for root in range(52360):
  ps=sorted(byroot[root]);assert all(not b.startswith(a) for a,b in zip(ps,ps[1:]));d=max(map(len,ps));assert sum(1<<(d-len(p)) for p in ps)==1<<d
 acceptance=json.load(open(ACCEPT));assert F(acceptance['nominal_g'])==F(11,100000) and F(acceptance['minimum_saved_log_margin'])==oldmin
 shards=[]
 for i in range(15):
  p=W/f'shard{i:02d}_input.json';obj={'shard_id':i,'canonical_owner':'canonical_index modulo15','source_target':'11/100000','next_target':'23/200000','sources':per_sources[i],'leaves':per_leaves[i],'new_math_status':'NOT_STARTED;no real newq/log/corners','same_canonical_all_aliases_one_worker':True};write(p,obj);shards.append({'shard_id':i,'input_path':str(p),'input_sha256':sha(p),'objects':len(per_sources[i]),'leaves':len(per_leaves[i]),'bytes':p.stat().st_size})
 write(W/'used_mass_supports.jsonl',{}) if False else None
 with open(W/'used_mass_supports.jsonl','x') as f:
  for node in sorted(support_used):line(f,supports[node])
 write(W/'source_locks.json',locks);plan={'source_objects':107530,'source_leaves':111201,'worker_count':15,'shards':shards,'q_object_ownership_unique':True,'leaf_cover_exactly_once':True,'same_global_old_forest_preserved':True,'public_log_ratio':'23/22','public_log_status':'NOT_COMPUTED;one new shared constant only afterprimaryreview','new_q_status':'NOT_COMPUTED','corner_status':'NOT_STARTED','old_margin_positive_and_qualified':True,'future_leaf_required_fields':['margin110_new','margin115=margin110_new-Ulog','margin110_new>=saved_old_margin>0','pass115','q_new_rule153','single_anchor_closed16corners'],'future_full_domain_candidate':'g0*(1+min_all111201_margin110_new);remains complete evenif target115 haspending','automatic_H_split_or_cutoff_change':False,'wall_deadline':None,'source_mass_policy':'only oldH64g105 for EPS128old129;EPS64/EPS128full use uniform;no gstar substitution or extra rawL'};write(W/'finite_plan.json',plan)
 summary={'status':'155_joint153_and_g115_PREPARED_ONLY','pid':os.getpid(),'objects':107530,'leaves':111201,'roots':52360,'oldq129_checks':oldq_checks,'source_rule_counts':dict(classes),'old_mass_support_cells':len(support_used),'source_counts':[len(r) for r in per_sources],'leaf_counts':[len(r) for r in per_leaves],'serialized_shard_bytes':sum(r['bytes'] for r in shards),'actual_new_q_log_H_gradient_corner_calls':0,'oldq_corner_baseline_runs':0,'future_workers':15,'wall_seconds':time.monotonic()-began,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};write(W/'summary.json',summary)
 for path,h in locks.items():assert sha(path)==h
 print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
