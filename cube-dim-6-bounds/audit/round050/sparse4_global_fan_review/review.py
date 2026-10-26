"""F050 independent saved-ledger/fan audit. No source writer or linear solve.
Active records use reduced affine 5x5 integer determinants, not original6x6.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
from functools import lru_cache
from collections import Counter,defaultdict
from math import gcd,lcm,prod,comb
import gzip,json,hashlib,time,os,platform
from flint import fmpz_mat
P=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
 p=Path(p);return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def ids(mask):
 while mask:
  b=mask&-mask;yield b.bit_length()-1;mask-=b

def main():
 assert not(P/'execution_receipt.json').exists(),'one complete ledger pass only'
 clock=read(P/'input_manifest.json');src=Path(clock['source_directory']);tic=time.monotonic();ops=Counter();parents=[]
 receipt={'status':'RUNNING','pid':os.getpid(),'started_unix':time.time(),'active_records_completed':0,'new_solves':0,'new_vertices':0,'first_parent_reenumerations':0}
 def save():
  receipt['operation_counts']=dict(ops);receipt['operation_total']=sum(ops.values());receipt['elapsed_seconds']=time.monotonic()-tic
  (P/'execution_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 def tick(k,n=1):
  ops[k]+=n
  assert sum(ops.values())<=2000000,'operation cap'
  assert time.time()-clock['started_unix']<1200,'time cap'
 def matrix_det(a):tick('other_determinants');return int(fmpz_mat(a).det())
 try:
  save()
  for f,h in clock['expected_sha256'].items():assert sha(f)==h,f
  manifest=read(src/'manifest.json')
  for f,h in manifest['files_sha256'].items():assert sha(src/f)==h,f
  reg=read(src/'registration.json')
  for f,h in reg['source_sha256'].items():assert sha(f)==h,f
  frame=read(reg['frame_path']);cols=frame['column_numerators'];assert len(cols)==152 and frame['common_denominator']==80
  first=read(reg['inherited_parent1_fan']);assert sha(reg['inherited_parent1_fan'])=='bae26cd08d31dac644e5ae0e77483ae912de0f687a03d6bc3d7b98a4a35a78b6'
  inherited_audit=P.parents[1]/'round046/sparse4_local_fan_review/manifest.json'
  assert sha(inherited_audit)=='6b347c5b9aa900f11893afbd0a7e706fda066c7891496d119eed7226bdd59edd'
  # Old geometry is inherited, not solved. Read exact original parent rays.
  oldpath=next(k for k in reg['source_sha256'] if k.endswith('fan_integer_input.txt'))
  tokens=iter(map(int,Path(oldpath).read_text().split()));assert next(tokens)==21;old=[]
  for p in range(21):
   d=next(tokens);rr=[]
   for j in range(6):next(tokens);rr.append(tuple(next(tokens) for k in range(6)))
   old.append((d,rr))
  assert next(tokens,None)is None
  globalfan=read(src/'global_fan.json.gz');assert hashlib.sha256(gzip.decompress((src/'global_fan.json.gz').read_bytes())).hexdigest()==clock['decompressed_global_sha256']
  gr=globalfan['rays'];assert len(gr)==478
  priority={tuple(v['primitive_physical_ray']):v['id'] for v in first['vertices']}
  globalrays=[tuple(v['primitive']) for v in gr];assert len(set(globalrays))==478 and [v['id'] for v in gr]==list(range(478))
  order=lambda r:(0,priority[r]) if r in priority else (1,r)
  assert globalrays==sorted(globalrays,key=order)
  rates={};formvalues={}
  for r,v in zip(globalrays,gr):
   tick('primitive_ray_checks');assert gcd(*r)==1 and r[0]>0 and r[5]>=0 and all(r[j]>=r[j+1] for j in range(5))
   values=[]
   for a in cols:tick('global_ray_form_dots');values.append(dot(a,r))
   q=Q(sum(map(abs,values)),80);assert q>0 and q==Q(v['q_true_152']);rates[r]=q;formvalues[r]=values
  rmap={r:k for k,r in enumerate(globalrays)};fans={1:first};active_all=Counter()
  for pid in range(2,22):
   fan=read(src/f'parents/{pid:02d}.fan.json.gz');fans[pid]=fan
   rays=list(map(tuple,fan['parent_rays']));pd=matrix_det(list(zip(*rays)));assert rays==old[pid-1][1] and abs(pd)==old[pid-1][0]==Q(fan['parent_determinant'])
   # Rebind H to all true new forms, rather than accepting saved boundaries.
   crossings=[];cuts=[]
   for k,a in enumerate(cols[32:]):
    tick('parent_crossing_form_tests');v=tuple(dot(a,r) for r in rays)
    if min(v)<0<max(v):
     crossings.append(k+1);g=gcd(*v);cuts.append(tuple(z//g for z in v))
   assert crossings==fan['crossing_new_form_indices_1_based']
   assert len(set(cuts))==len(cuts) and not any(tuple(-z for z in h) in cuts for h in cuts)
   H=[tuple(int(i==j) for i in range(6)) for j in range(6)]+cuts
   assert [list(h) for h in H]==fan['boundary_forms_in_u']
   vv=fan['vertices'];U=[tuple(map(Q,v['u'])) for v in vv];rr=[tuple(v['primitive_physical_ray']) for v in vv];N=len(U)
   assert [v['id'] for v in vv]==list(range(N)) and len(set(U))==len(set(rr))==N
   assert rr==sorted(rr,key=order)
   nums=[];dens=[];scales=[]
   for u,r,v in zip(U,rr,vv):
    tick('local_vertex_relations');assert min(u)>=0 and sum(u)==1 and r in rmap
    d=lcm(*(z.denominator for z in u));n=tuple(int(z*d) for z in u);nums.append(n);dens.append(d)
    scale=Q(v['slice_scale']);scales.append(scale);assert scale>0
    assert all(dot([ray[k] for ray in rays],u)==scale*r[k] for k in range(6))
    assert Q(v['q_true_152'])==rates[r]
   enum={v['enumeration_id']:k for k,v in enumerate(vv)};assert set(enum)==set(range(N))
   # Single full ordered saved ledger pass; no solve/new point is performed.
   ledger=src/f'parents/{pid:02d}.active.jsonl.gz';sequence=iter(combinations(range(len(H)),5));count=Counter();seen=set();digest=hashlib.sha256()
   with gzip.open(ledger,'rb') as handle:
    for line in handle:
     digest.update(line);act,kind,*pay=json.loads(line);assert tuple(act)==next(sequence)
     tick('active_det_and_witness_records');receipt['active_records_completed']+=1
     assert receipt['active_records_completed']<=783994
     reduced=[[H[k][j]-H[k][5] for j in range(5)] for k in act]
     determinant=int(fmpz_mat(reduced).det());count[kind]+=1
     if kind=='S':assert determinant==0 and not pay
     elif kind=='O':
      assert determinant and len(pay)==1;u=tuple(map(Q,pay[0]));assert len(u)==6 and sum(u)==1 and min(u)<0
      d=lcm(*(z.denominator for z in u));n=tuple(int(z*d) for z in u);assert all(dot(H[k],n)==0 for k in act)
     else:
      assert kind=='I' and determinant and len(pay)==1;eid=pay[0];k=enum[eid];assert all(dot(H[j],nums[k])==0 for j in act)
      if eid not in seen:assert eid==len(seen) and act==vv[k]['active_witness'];seen.add(eid)
   assert next(sequence,None)is None and sum(count.values())==comb(len(H),5) and seen==set(enum)
   st=fan['active_enumeration'];assert count==Counter(S=st['singular'],O=st['outside'],I=st['inside_systems'])
   assert sum(count.values())==st['expected_active_systems']==st['active_completed']
   assert digest.hexdigest()==st['active_ledger_uncompressed_sha256'] and sha(ledger)==st['active_ledger_sha256'];active_all.update(count)
   allmask=(1<<N)-1;zero=[];positive=[];negative=[]
   for h in H:
    tick('halfspace_masks');vals=[dot(h,n) for n in nums]
    zero.append(sum(1<<k for k,z in enumerate(vals) if z==0));positive.append(sum(1<<k for k,z in enumerate(vals) if z>=0));negative.append(sum(1<<k for k,z in enumerate(vals) if z<=0))
   @lru_cache(None)
   def dim(mask):
    tick('exact_ranks');return int(fmpz_mat([nums[k] for k in ids(mask)]).rank())-1 if mask else -1
   def maskof(v):return sum(1<<k for k in v)
   cells=fan['full_cells'];codes={v['sign_code']:v for v in cells};assert len(codes)==len(cells);m=len(cuts)
   graph={};cellmask={}
   for c in cells:
    tick('closed_cell_filters');code=c['sign_code'];assert 0<=code<(1<<m)
    mask=allmask
    for j in range(m):mask&=positive[6+j] if code>>j&1 else negative[6+j]
    assert list(ids(mask))==c['vertices'] and dim(mask)==5
    assert c['signs']==[1 if code>>j&1 else -1 for j in range(m)];cellmask[code]=mask
   for code,c in codes.items():
    adj=[];mask=cellmask[code]
    for j in range(m):
     t=mask&zero[6+j]
     if dim(t)==4:
      tick('forced_neighbor_checks');other=code^(1<<j);assert other in codes and mask&cellmask[other]==t
      adj.append({'cut_index':j,'neighbor_sign_code':other,'facet_vertices':list(ids(t))})
    assert adj==c['adjacency'];graph[code]=[a['neighbor_sign_code'] for a in adj]
   generic=st['generic_seed_positive_integers'];assert len(generic)==6 and min(generic)>0
   z=generic[1];assert 1<=z<=5*m+1 and generic==[z**j for j in range(6)]
   signvalues=[dot(h,generic) for h in cuts];assert all(signvalues)
   seed=sum((s>0)<<j for j,s in enumerate(signvalues));assert seed==fan['bfs_seed_sign_code'] and seed in codes and fan['bfs_queue_exhausted']
   reached={seed};todo=[seed]
   while todo:
    for other in graph[todo.pop()]:
     if other not in reached:reached.add(other);todo.append(other)
   assert reached==set(codes) and allmask==__import__('functools').reduce(int.__or__,cellmask.values(),0)
   # Topological proof DAG, exact facets and bottom-up pulling.
   nodes={maskof(v['vertices']):v for v in fan['pulling_face_certificates']};assert len(nodes)==len(fan['pulling_face_certificates'])
   triangles={}
   for mask,n in sorted(nodes.items(),key=lambda item:(item[1]['dimension'],item[0])):
    tick('pulling_nodes');d=dim(mask);vs=tuple(ids(mask));assert d==n['dimension'] and n['apex']==vs[0]
    proper={mask&z for z in zero}-{0,mask}
    facets={t for t in proper if not any(t!=u and t&u==t for u in proper)} if d else set()
    assert all(dim(t)==d-1 for t in facets) and facets=={maskof(v) for v in n['facets']}
    if not d:expected={(vs[0],)}
    else:
     opposite=[t for t in facets if not(t>>vs[0]&1)];assert opposite and all(t in triangles for t in opposite)
     expected={tuple(sorted((vs[0],)+simplex)) for t in opposite for simplex in triangles[t]}
    actual={tuple(sorted(t)) for t in n['simplices']};assert actual==expected and len(actual)==len(n['simplices'])
    assert all(len(t)==d+1 and dim(maskof(t))==d for t in actual);triangles[mask]=actual
   simp=fan['simplices'];seen_simp=set();total=Q(0)
   for code,c in codes.items():
    mask=cellmask[code];vertices=c['vertices'];sg=[]
    for j in range(152):
     tick('cell_true_form_sign_tests');values=[formvalues[rr[k]][j] for k in vertices];assert not(min(values)<0<max(values));sg.append(1 if max(values)>0 else -1 if min(values)<0 else 0)
    slope=tuple(Q(sum(sg[j]*cols[j][k] for j in range(152)),80) for k in range(6));assert slope==tuple(map(Q,c['linear_symbol_coefficient']))
    assert all(dot(slope,rr[k])==rates[rr[k]] for k in vertices)
    localset={tuple(sorted(simp[j]['vertices'])) for j in c['simplex_indices']};assert localset==triangles[mask]
    cv=Q(0)
    for sid in c['simplex_indices']:
     tick('child_inclusion_and_sign_inheritance');child=simp[sid];ks=child['vertices'];key=tuple(sorted(ks));assert len(ks)==6 and key not in seen_simp;seen_simp.add(key)
     assert child['id']==sid and child['cell']==c['id'] and all(k in vertices for k in ks)
     phys=matrix_det(list(zip(*(rr[k] for k in ks))));un=matrix_det(list(zip(*(nums[k] for k in ks))))
     du=Q(un,prod(dens[k] for k in ks));vol=abs(du)
     assert phys and vol>0 and Q(child['primitive_physical_determinant'])==phys and Q(child['relative_slice_volume'])==vol
     assert phys*prod(scales[k] for k in ks)==pd*du
     cv+=vol
    assert cv==Q(c['relative_slice_volume']);total+=cv
   assert total==1 and len(seen_simp)==len(simp)
   parents.append({'parent':pid,'active_records':sum(count.values()),'classifications':dict(count),'ledger_uncompressed_sha256':digest.hexdigest(),'vertices':N,'full_cells':len(cells),'pulling_nodes':len(nodes),'simplices':len(simp),'forced_adjacencies':sum(map(len,graph.values())),'relative_volume_exact':'1','status':'PASS'})
   receipt['completed_parents']=[x['parent'] for x in parents];save();print(json.dumps({'parent':pid,'active_records':receipt['active_records_completed'],'operations':sum(ops.values()),'seconds':time.monotonic()-tic}),flush=True)
  assert receipt['active_records_completed']==783994 and active_all==Counter(S=247176,O=319333,I=217485)
  # Inherited parent1 is only hash/certificate data; no ledger/rank/ray re-enumeration.
  facets=defaultdict(list);flat=[];slice_norm=Q(0);seen_global=set()
  for pid,fan in sorted(fans.items()):
   vv=fan['vertices'];pr=list(map(tuple,fan['parent_rays']));assert pr==old[pid-1][1]
   for child in fan['simplices']:
    tick('global_cone_assembly');rs=[tuple(vv[k]['primitive_physical_ray']) for k in child['vertices']];d=abs(int(child['primitive_physical_determinant']));assert d>0
    key=tuple(sorted(rs));assert key not in seen_global;seen_global.add(key)
    flat.append({'parent':pid,'local_id':child['id'],'abs_primitive_determinant':d,'ray_ids':[rmap[r] for r in rs]})
    assert all(Q(vv[k]['q_true_152'])==rates[r] for k,r in zip(child['vertices'],rs));slice_norm+=Q(d,prod(r[0] for r in rs))
   for j in range(6):
    tick('old_parent_facet_relations');on={k for k,v in enumerate(vv) if Q(v['u'][j])==0};induced=set()
    for child in fan['simplices']:
     v=tuple(sorted(tuple(vv[k]['primitive_physical_ray']) for k in set(child['vertices'])&on))
     if len(v)==5:induced.add(v)
    assert induced
    facets[tuple(sorted(pr[:j]+pr[j+1:]))].append((pid,j,induced))
  assert flat==globalfan['cones'] and len(flat)==13229 and slice_norm==1
  assert set(globalrays)==set(r for f in fans.values() for r in map(lambda v:tuple(v['primitive_physical_ray']),f['vertices']))
  walls=[tuple(int(k==j)-int(k==j+1) for k in range(6)) for j in range(5)]+[(0,0,0,0,0,1)]
  shared=[];boundary=[]
  for key,owners in sorted(facets.items()):
   tick('cross_parent_facet_checks');assert len(owners) in (1,2)
   if len(owners)==2:
    a,b=owners;assert a[2]==b[2]
    shared.append({'parents':[a[0],b[0]],'parent_facets':[a[1],b[1]],'old_facet_rays':[list(x) for x in key],'triangulation':[[list(r) for r in q] for q in sorted(a[2])]})
   else:
    a=owners[0];which=[j for j,w in enumerate(walls) if all(dot(w,r)==0 for r in key)];assert which
    assert all(all(dot(walls[j],r)==0 for r in simplex) for simplex in a[2] for j in which)
    boundary.append({'parent':a[0],'parent_facet':a[1],'chamber_walls':which,'old_facet_rays':[list(x) for x in key],'triangulation':[[list(r) for r in q] for q in sorted(a[2])]})
  assert len(shared)==28 and len(boundary)==70 and 2*len(shared)+len(boundary)==126
  assert read(src/'cross_parent_faces.json.gz')=={'shared_facets':shared,'boundary_facets':boundary,'all_matched':True}
  assert globalfan['full_signed_permutation_multiplicity']==46080 and globalfan['logical_terms_per_future_full_kernel']==609592320==46080*len(flat)
  assert globalfan['ordered_xi1_slice_volume']=='1/120' and globalfan['global_negative_pair_factor_not_applied'] and not globalfan['numeric_engine_int64_or_zero_node_compatibility_claimed']
  run=read(src/'run_receipt.json');checkpoint=read(src/'checkpoint.json');assert run['completed_parents']==list(range(1,22))
  for name in ['uncovered_parents','uncovered']: 
   if name in checkpoint:assert checkpoint[name]==[]
  for f,h in manifest['files_sha256'].items():assert sha(src/f)==h,('drift',f)
  out={'status':'PASS_COMPLETE_MATHEMATICAL_FAN_ONLY','parents':parents,'active_ledger_records':783994,'active_classifications':dict(active_all),'new_solves':0,'new_vertices':0,'first_parent_reenumerations':0,'operation_counts':dict(ops),'matrix_and_structure_operations':sum(ops.values()),'operation_accounting':'one exact determinant plus witness-validation bundle per ledger record; each cache-miss rank/determinant, ray/form dot, cell/form sign test and listed structural object separately counted; elementary arithmetic is not a matrix/structure operation','global_rays':478,'full_cells':12+sum(p['full_cells'] for p in parents),'full_simplices':13229,'shared_old_facets':28,'outer_wall_facets':70,'xi1_slice_volume_exact':'1/120','future_full_group_logical_terms':609592320,'numeric_engine_released':False,'new_S2_or_tube_inherited':False,'new_C6_bound':False,'pid':os.getpid(),'runtime_seconds':time.monotonic()-tic,'audit_elapsed_seconds':time.time()-clock['started_unix'],'source_sha256':sha(__file__),'author_manifest_sha256':sha(src/'manifest.json'),'kernel_calls':0,'MATLAB_calls':0}
  assert out['full_cells']==4266
  (P/'results.json').write_text(json.dumps(out,indent=2)+'\n');receipt.update(status='PASS',finished_unix=time.time());save();print(json.dumps({k:out[k] for k in ['status','matrix_and_structure_operations','runtime_seconds','audit_elapsed_seconds']}))
 except BaseException as exc:
  receipt.update(status='FAILED_STOP_NO_RETRY',error=repr(exc),finished_unix=time.time());(P/'partial_parent_results.json').write_text(json.dumps(parents,indent=2)+'\n');save();raise
if __name__=='__main__':main()
