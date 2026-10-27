"""Conditional terminal coverage/unknown/budget/source structure. Exact saved arithmetic only, no Arb/native/old validator."""
import sys,json,hashlib,collections,itertools
from pathlib import Path
from fractions import Fraction as Q
sys.set_int_max_str_digits(0);sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;A=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-numerics/cube-dim-6-bounds/upper/round129/recovered_continuation');REC=A.parents[2]/'audit/round127/crash_recovery'
rd=lambda p:json.loads(Path(p).read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
checks=0
def ck(ok,label):
 global checks
 if not ok:raise AssertionError(label)
 checks+=1
seed=rd(A/'continuation_seed129.json');forest=rd(A/'forest_v129.json');budget=rd(REC/'root_recovery_and_budget.json');geo=rd(A/'unknown_geometry_receipt.json');unknown=rd(REC/'unknown_no_reissue.json');oldcache=rd(REC/'combined_cache.json');newidx=rd(A/'semantic_runtime_index_v129.json')['entries'];cidx=rd(A/'certificate_index_v129.json')['entries'];prepared=rd(P/'prepared_new_jobs.json')
ck(sha(A/'terminal_source_pins.json')if False else prepared['terminal_source_pins_sha256']==sha(P/'terminal_source_pins.json'),'prepared source pin binding')
for row in rd(P/'terminal_source_pins.json')['files']:ck(sha(row['path'])==row['sha256'],'terminal source unchanged')
continued=set(forest['continued_parent_roots']);ck(len(seed['roots'])==len(forest['roots'])==325 and set(seed['roots'])==set(forest['roots'])and continued==set(budget)and len(continued)==55,'same325 roots and55 continued roots')
ck(all(forest[k]==seed[k]for k in ['inherited_passed_count','inherited_passed_leaf_ids','new117_passed_leaves']),'inherited coverage metadata unchanged')
ck(all(forest['roots'][k]==seed['roots'][k]for k in set(seed['roots'])-continued),'all270 unrelated closed roots exactly preserved')
newstates={};artifact={}
for k,r in newidx.items():
 s=rd(r['state_path'])
 if s['status']=='AVAILABLE':newstates[k]=s;artifact[s['artifact_key']]=s
newcert={e['key']:e for e in cidx};certvalues={};data={};old_read=set();new_read=set();windows=[{'kind':'B','a':a,'T':t}for a,t in [('1/32','1/2'),('1/64','1/32'),('1/128','1/64'),('1/250','1/128'),('1/1000','1/250')]];T=Q(1,65000)
def dy(v):return tuple(Q(int(v[k]))*Q(2)**int(v['common_binary_exponent'])for k in ['lower_mantissa','upper_mantissa'])
def source(ak):
 if ak in data:return data[ak]
 if ak in artifact:s=artifact[ak];e=s['expected'];p=s['native_path'];h=s['native_sha256'];new_read.add(ak)
 else:
  k=oldcache['artifact_to_semantic'][ak];s=oldcache['entries'][k];e=s['expected'];p=s['native_path'];h=s['native_sha256'];old_read.add(ak)
 ck(sha(p)==h,'saved numeric constant bytes');r=rd(p);value=(e,dy(r['band']),[dy(v)for v in r['gradient']],dy(r['p_a']),dy(r['p_T']));data[ak]=value;return value

def prefix(n):
 p=n.get('partial_prefix');
 if not p:return None
 keys=p['artifact_keys'];ck(1<=len(keys)<=5 and len(set(keys))==len(keys),'one consecutive native band per prefix')
 ds=[source(k)for k in keys];center=['1']+[str((Q(l)+Q(u))/2)for l,u in n['box']]
 ck(all(d[0]['center']==center and d[0]['window']==windows[j]for j,d in enumerate(ds)),'prefix true current midpoint and exact windows')
 F=tuple(sum(d[1][k]for d in ds)for k in(0,1));g=[tuple(sum(d[2][i][k]for d in ds)for k in(0,1))for i in range(6)];inter=[]
 for j in range(1,len(ds)):
  l=max(ds[j-1][3][0],ds[j][4][0]);u=min(ds[j-1][3][1],ds[j][4][1]);ck(l<=u,'prefix adjacent density intersection');inter.append(list(map(str,(l,u))))
 a=Q(windows[len(ds)-1]['a']);end=Q(1,2);terms=[(5*F[0],5*F[1]),(end*ds[0][4][0],end*ds[0][4][1]),(-a*ds[-1][3][1],-a*ds[-1][3][0])]+[tuple(sorted(Q(x)*z for z in v))for x,v in zip(center,g)];E=tuple(sum(t[k]for t in terms)for k in(0,1));ck(E[0]<=0<=E[1],'merged exact prefix Euler')
 ck(p=={'artifact_keys':keys,'F_interval':list(map(str,F)),'six_gradients':[list(map(str,v))for v in g],'outer_window':{'kind':'B','a':str(a),'T':'1/2'},'Euler_interval':list(map(str,E)),'internal_density_intersections':inter},'all saved exact prefix arithmetic')
 return keys,F

def bounds(n):
 r=n['certificate_refs'];ck(r['cache']=='PER_BAND_SUM','new C129 box uses explicit band sum');keys=n['partial_prefix']['artifact_keys'];ck(r['keys']==keys and len(r['certificates'])==len(keys),'same physical bands in box sum')
 l=u=Q(0)
 for j,(ak,ref)in enumerate(zip(keys,r['certificates'])):
  ck(ref['cache']=='NEW'and ref['key']in newcert,'only registered new geometry expressions here');e=newcert[ref['key']];ck(e['native_artifact_key']==ak and e['path']==ref['path']and e['expected']['face_box']==n['box']and e['expected']['window']==windows[j],'new certificate source/box/window reference')
  if e['key']not in certvalues:
   ck(sha(e['path'])==e['sha256'],'new certificate pinned bytes');c=rd(e['path']);certvalues[e['key']]=(Q(c['tangent_bound_lower']),Q(c['tangent_bound_upper']))
  a,b=certvalues[e['key']];l+=a;u+=b
 ck(list(map(Q,n['whole_box_interval']))==[l,u],'exact disjoint band box sum')
 return l,u

def tighten(raw):
 lo=[max(Q(raw[j][0])for j in range(i,5))for i in range(5)];hi=[min(Q(raw[j][1])for j in range(i+1))for i in range(5)];empty=any(a>b for a,b in zip(lo,hi));return None if empty else[[str(a),str(b)]for a,b in zip(lo,hi)]
def split(n,nodes):
 box=n['box'];axis=max(range(5),key=lambda j:(Q(box[j][1])-Q(box[j][0]),-j));m=(Q(box[axis][0])+Q(box[axis][1]))/2;lo=[b.copy()for b in box];hi=[b.copy()for b in box];lo[axis][1]=str(m);hi[axis][0]=str(m);s=n['split']
 ck(s['axis']==axis and Q(s['seam'])==m and s['lower_raw_box']==lo and s['upper_raw_box']==hi and s['closed_seam']is True,'closed exact raw split geometry')
 ck(n['children']==[n['id']+'0',n['id']+'1'],'exact two child IDs')
 for cid,raw in zip(n['children'],[lo,hi]):
  c=nodes[cid];ck(c['parent']==n['id']and c['depth']==n['depth']+1 and c['raw_box']==raw and c['box']==tighten(raw),'ordered tightening preserves each closed child domain')
  ck((c['status']=='empty')==(c['box']is None),'empty iff ordered domain excluded')

statuses=collections.Counter();changed=processed=created=0;newpassed=0;oldclosed=0;perroot=[];prefix_nodes=0
for rid in sorted(continued):
 old=seed['roots'][rid];r=forest['roots'][rid];nodes=r['nodes'];initial=set(old['queue']+old['blocked_leaves']);ck(len(initial)==len(old['queue'])+len(old['blocked_leaves']),'no duplicate seed frontier')
 ck(set(old['nodes'])<=set(nodes),'all seed nodes retained')
 for nid,n in old['nodes'].items():
  if nid not in initial:ck(nodes[nid]==n,'old passed/empty/split nodes immutable')
 if old['root_complete']:
  oldclosed+=1;ck(r['root_complete']and r['passed_leaves']==old['passed_leaves']and r['empty_leaves']==old['empty_leaves'],'six already closed roots stay closed')
 touched=set();leaf_pass=set();leaf_empty=set();leaf_unpaid=set();todo=list(initial);root_processed=0
 while todo:
  nid=todo.pop();ck(nid not in touched,'no cycle or duplicate child');touched.add(nid);n=nodes[nid];status=n['status'];statuses[status]+=1
  if n!=old['nodes'].get(nid):changed+=1
  if status not in['queued','empty']:root_processed+=1
  pref=prefix(n)
  if pref:prefix_nodes+=1
  if status=='split':
   ck(pref is not None and pref[1][0]>T,'center prefix passed before C129 split');l,u=bounds(n);ck(l<T,'C129 split follows failed box lower');split(n,nodes);todo.extend(n['children'])
  elif status=='passed':
   ck(pref is not None and pref[1][0]>T,'center prefix for new passed leaf');l,u=bounds(n);ck(l>=T and not n.get('children'),'complete new passed leaf bound');leaf_pass.add(nid);newpassed+=1
   cover=[v for v in r['closed_coverage']if v['node']==nid];ck(len(cover)==1 and cover[0]['box']==n['box']and list(map(Q,cover[0]['interval']))==[l,u]and cover[0]['refs']==n['certificate_refs'],'new closed coverage proof binding')
  elif status=='empty':leaf_empty.add(nid)
  else:ck(status in['queued','blocked'],'only honest unpaid frontier status');leaf_unpaid.add(nid)
 ck(set(nodes)-set(old['nodes'])<=touched,'no orphan new nodes');created+=len(set(nodes)-set(old['nodes']));processed+=root_processed
 ck(set(r['passed_leaves'])==set(old['passed_leaves'])|leaf_pass and set(r['empty_leaves'])==set(old['empty_leaves'])|leaf_empty,'paid/empty leaf lists exact old plus new')
 ck(set(r['queue']+r['blocked_leaves'])==leaf_unpaid and len(r['queue'])+len(r['blocked_leaves'])==len(leaf_unpaid)and r['residual_frontier']==r['queue']+r['blocked_leaves'],'all unpaid child domains retained exactly')
 ck(r['root_complete']==(not leaf_unpaid)and r['inflight']is None,'root closure iff no unpaid domain')
 ck(r['native_reserved_129']==budget[rid]['native_used_127']+prepared['per_root_new_reservations'].get(rid,0)<=512,'cumulative127+129 per-root reservation budget')
 debit=old.get('recovery_geometry_nodes',0);ck(r['nodes_processed_129']==old['nodes_processed_127']+debit+root_processed<=1024,'node budget includes prior127 and geometry exactly once')
 if not r['root_complete']:ck(r['status']=='NATIVE_BUDGET_EXHAUSTED'and r['native_reserved_129']==512,'35 honest root budget stops')
 perroot.append({'root_id':rid,'root_complete':r['root_complete'],'native127plus129':r['native_reserved_129'],'new_processed':root_processed,'new_passed_leaves':len(leaf_pass),'unpaid_leaves':len(leaf_unpaid)})
# Every old unknown parent was replaced by its exact geometry-seed children and never altered/reissued.
for record in geo['records']:
 rid=record['root_id'];nid=record['unknown_parent_node'];ck(forest['roots'][rid]['nodes'][nid]==seed['roots'][rid]['nodes'][nid]and seed['roots'][rid]['nodes'][nid]['unknown_attempt_preserved']in unknown,'15 untouched quarantined parent/closed-child bridges')
ck(len(certvalues)==18257 and len(new_read)==18257 and len(old_read)==4,'all new sources/certificates and only four reused old numeric constants accounted')
ck(sum(r['root_complete']for r in perroot)==20 and oldclosed==6 and len(perroot)==55,'20 closed continued roots include6 previously closed')
ck(forest['global_g_certificate']is False and forest['new_C6_bound']is False,'partial certificate no new global/C6 claim')
out={'status':'PASS_TERMINAL_STRUCTURE_AND_SAVED_ARITHMETIC_PENDING_NEW_ARB768_ENDPOINTS','assertions':checks,'continued_roots':55,'closed_continued_roots':20,'previously_closed_continued_roots_preserved':6,'newly_closed_roots':14,'budget_stopped_roots':35,'total_closed_of325':290,'new_processed_nodes':processed,'new_created_nodes':created,'new_passed_leaves':newpassed,'frontier_status_counts':dict(statuses),'all_registered_new_certificates_used':len(certvalues),'new_raw_constants_referenced':len(new_read),'old_raw_constants_referenced_only_not_revalidated':len(old_read),'unknown15_bridges_preserved':True,'new_Arb_native_old_baseline_validator_calls':0,'per_root':perroot}
(P/'structure_review.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k!='per_root'},indent=2))
