"""F046 one 462-system audit using affine elimination and Fraction arithmetic.
No author code is imported or executed. No kernel or old fan enumeration.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations,product
from math import gcd,lcm,prod
from functools import lru_cache
import json,hashlib,time,os,platform
P=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def rref(rows,n):
 a=[[Q(x) for x in row] for row in rows];r=0;piv=[]
 for c in range(n):
  k=next((k for k in range(r,len(a)) if a[k][c]),None)
  if k is None:continue
  a[r],a[k]=a[k],a[r];z=a[r][c];a[r]=[x/z for x in a[r]]
  for k in range(len(a)):
   if k!=r and a[k][c]:
    z=a[k][c];a[k]=[x-z*y for x,y in zip(a[k],a[r])]
  piv.append(c);r+=1
  if r==len(a):break
 return a,piv

def det(rows):
 a=[[Q(x) for x in row] for row in rows];n=len(a);d=Q(1)
 for c in range(n):
  k=next((k for k in range(c,n) if a[k][c]),None)
  if k is None:return Q(0)
  if k!=c:a[c],a[k]=a[k],a[c];d=-d
  p=a[c][c];d*=p
  for k in range(c+1,n):
   z=a[k][c]/p
   for j in range(c+1,n):a[k][j]-=z*a[c][j]
 return d

def detcols(cols):return det(list(zip(*cols)))
def maximum_sets(items):
 items=set(tuple(sorted(x)) for x in items if x)
 return {x for x in items if not any(set(x)<set(y) for y in items)}

def main():
 if (P/'execution_receipt.json').exists():raise SystemExit('one full enumeration only')
 clock=json.loads((P/'input_manifest.json').read_text());src=Path(clock['source_directory']);tic=time.monotonic()
 receipt={'status':'RUNNING','pid':os.getpid(),'started_unix':time.time(),'active_systems_attempted':0,'structure_relation_groups':0,'kernel_calls':0}
 (P/'execution_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 try:
  for f,h in clock['source_sha256'].items():assert sha(f)==h,f
  m=json.loads((src/'manifest.json').read_text())
  for f,h in m['files_sha256'].items():assert sha(src/f)==h,f
  inp=json.loads((src/'input_manifest.json').read_text())
  for f,h in inp['source_sha256'].items():assert sha(f)==h,f
  bp=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round044/sparse4_frame_preflight')
  assert sha(bp/'manifest.json')=='b1a705ea2f4a06052a837642cfde6469907729baccba1cb58fc175aad5ad6306'
  bm=json.loads((bp/'manifest.json').read_text())
  for f,h in bm['whitelist_sha256'].items():assert sha(bp/f)==h,f
  frame=json.loads((bp/'frame.json').read_text());fan=json.loads((src/'local_fan.json').read_text());A=frame['column_numerators'];den=frame['common_denominator'];rels=[]
  assert den==80 and len(A)==152 and frame['beta']=='1/80' and fan['frame_sha256']==sha(bp/'frame.json')
  # 1. Legal columns, exact Gram and marginals; no old 2160 check.
  assert len(set(map(tuple,A)))==152
  assert all(a[0]==80 and all(abs(z)==80 for z in a) for a in A[:32])
  supports={}
  for a in A[32:]:
   ids=tuple(i for i,z in enumerate(a) if z)
   assert len(ids)==4 and a[ids[0]]==1 and all(abs(a[i])==1 for i in ids)
   supports[ids]=supports.get(ids,0)+1
  assert len(supports)==15 and set(supports.values())=={8}
  gram=[[Q(sum(a[i]*a[j] for a in A),den**2) for j in range(6)] for i in range(6)]
  assert gram==[[Q(2561,80) if i==j else Q(0) for j in range(6)] for i in range(6)]
  assert [Q(sum(abs(a[i]) for a in A),den) for i in range(6)]==[Q(33)]*6
  rels.append({'id':1,'name':'true_152_frame','Gram_diagonal_exact':'2561/80','marginal_rate_exact':'33','supports':15,'sign_patterns_per_support':8,'same_time_positive_Cauchy':'proved in proof_review'})
  # 2. Original first parent only; physical slice via independent edge Gram.
  old=Path(next(f for f in inp['source_sha256'] if f.endswith('fan_integer_input.txt')))
  lines=old.read_text().splitlines();rays=[tuple(map(int,l.split()[1:])) for l in lines[2:8]]
  assert fan['parent_rays']==[list(r) for r in rays] and detcols(rays)==1
  ell=(0,-1,1,1,0,0);assert all(dot(ell,r)==1 for r in rays)
  edges=[tuple(r[j]-rays[5][j] for j in range(6)) for r in rays[:5]]
  edgegram=[[dot(a,b) for b in edges] for a in edges]
  assert det(edgegram)==3 and dot(ell,ell)==3
  rels.append({'id':2,'name':'parent_and_physical_slice','parent_det':'1','physical_normal':list(ell),'physical_edge_Gram_det':'3','u_volume':'sqrt(6)/120','physical_volume':'sqrt(3)/120','area_Jacobian':'1/sqrt(2)'})
  # 3. All local new forms only; exact interior crossing witness.
  pd=[tuple(dot(a,r) for r in rays) for a in A[32:]]
  cross=[j+1 for j,v in enumerate(pd) if min(v)<0<max(v)]
  assert cross==[2,10,26,50,82]==fan['crossing_new_form_indices_1_based']
  H=[tuple(int(i==j) for i in range(6)) for j in range(6)]+[pd[j-1] for j in cross]
  assert [list(h) for h in H]==fan['boundary_forms_in_u']
  assert all(h[4:]==(0,0) for h in H[6:])
  plus=[2,1,1,1,1,1];minus=[1,1,2,1,1,1]
  assert dot(H[6],plus)==1 and dot(H[6],minus)==-1
  xp=tuple(dot([r[j] for r in rays],plus) for j in range(6));xm=tuple(dot([r[j] for r in rays],minus) for j in range(6))
  assert all(xp[j]>xp[j+1] for j in range(5)) and xp[5]>0
  assert all(xm[j]>xm[j+1] for j in range(5)) and xm[5]>0
  rels.append({'id':3,'name':'first_parent_true_split','crossings':cross,'local_new_form_dots':720,'strict_interior_plus':list(xp),'strict_interior_minus':list(xm),'symbol_midpoint_defect_lower':'1/80','q_only_update_valid':False})
  # 4. ONE full 462 enumeration, eliminate u6 first, solve 5x5 by pure Fraction RREF.
  vertices=set();audit=[];sing=outside=0
  for active in combinations(range(11),5):
   receipt['active_systems_attempted']+=1
   assert time.time()-clock['started_unix']<900
   equations=[[H[k][j]-H[k][5] for j in range(5)]+[-H[k][5]] for k in active]
   rr,piv=rref(equations,5)
   if len(piv)<5:sing+=1;audit.append({'active':active,'status':'singular'});continue
   y=tuple(row[5] for row in rr);u=y+(1-sum(y),)
   assert sum(u)==1 and all(dot(H[k],u)==0 for k in active)
   if min(u)<0:outside+=1;status='outside'
   else:vertices.add(u);status='feasible'
   audit.append({'active':active,'status':status,'u':[str(x) for x in u]})
  U=sorted(vertices);assert receipt['active_systems_attempted']==462 and sing==237 and outside==6 and len(U)==12
  savedU=[tuple(Q(x) for x in v['u']) for v in fan['vertices']];assert U==savedU
  rels.append({'id':4,'name':'complete_independent_active_enumeration','systems':462,'algorithm':'eliminate u6; Fraction 5x5 RREF, not author FLINT6x6 solve','singular':sing,'outside':outside,'distinct_vertices':len(U)})
  @lru_cache(None)
  def dimension(ids):return len(rref([U[i] for i in ids],6)[1])-1 if ids else -1
  # 5. Complete patterns, closure witnesses, common-edge ties retained.
  pats=fan['closed_patterns'];assert len(pats)==32
  full=[]
  for j,signs in enumerate(product((-1,1),repeat=5)):
   p=pats[j];assert p['id']==j and p['signs']==list(signs)
   ids=tuple(k for k,u in enumerate(U) if all(s*dot(h,u)>=0 for s,h in zip(signs,H[6:])))
   assert list(ids)==p['vertices'] and dimension(ids)==p['affine_dimension']
   assert ids and {0,1}<=set(ids)
   if dimension(ids)==5:full.append(p)
  assert len(full)==12
  for p in pats:
   if p['affine_dimension']<5:
    q=pats[p['containing_full_pattern']];assert q in full and set(p['vertices'])<=set(q['vertices'])
  rels.append({'id':5,'name':'all_closed_patterns','full':12,'lower_nonempty':20,'empty':0,'all_boundary_patterns_retained':True})
  # 6. Facets found by maximal zero sets, not author's dimension filter.
  records={tuple(d['vertices']):d for d in fan['pulling_face_certificates']};assert len(records)==66
  facets={}
  for ids,d in records.items():
   dim=dimension(ids);assert dim==d['dimension'] and d['apex']==min(ids)
   zeros=[tuple(i for i in ids if dot(h,U[i])==0) for h in H]
   fs=maximum_sets(T for T in zeros if T and len(T)<len(ids)) if dim else set()
   assert all(dimension(T)==dim-1 for T in fs)
   assert fs==set(map(tuple,d['facets']));facets[ids]=fs
  rels.append({'id':6,'name':'complete_facets_from_maximal_zero_faces','proof_nodes':66,'all_saved_facets_match':True})
  # 7. Bottom-up pulling recurrence, all roots connected to recursively complete bases.
  tri={}
  for ids in sorted(records,key=lambda z:(dimension(z),z)):
   apex=min(ids);d=dimension(ids)
   if d==0:computed={(apex,)}
   else:
    bases=[T for T in facets[ids] if apex not in T];assert bases
    assert all(T in tri for T in bases)
    computed={tuple(sorted((apex,)+t)) for T in bases for t in tri[T]}
   assert computed==set(tuple(sorted(t)) for t in records[ids]['simplices'])
   assert len(computed)==len(records[ids]['simplices'])
   assert all(len(t)==d+1 and dimension(t)==d for t in computed)
   tri[ids]=computed
  for p in full:
   expected={tuple(sorted(fan['simplices'][j]['vertices'])) for j in p['simplex_indices']}
   assert tri[tuple(p['vertices'])]==expected
  rels.append({'id':7,'name':'bottom_up_pulling_complete','root_full_cells':12,'simplex_count':14,'coverage_basis':'closed partition plus facet/pyramid dimension induction; not volume'})
  # 8. Stronger shared-face check: all maximal induced faces, including lower ties.
  shared=[]
  for p,q in combinations(full,2):
   common=tuple(sorted(set(p['vertices'])&set(q['vertices'])));assert common and dimension(common)<5
   def induced(cell):return maximum_sets(set(t)&set(common) for t in tri[tuple(cell['vertices'])])
   left,right=induced(p),induced(q);assert left==right and left
   assert all(dimension(t)==dimension(common) for t in left)
   shared.append({'patterns':[p['id'],q['id']],'vertices':common,'maximal_induced_simplices':sorted(left)})
  assert len(shared)==66
  rels.append({'id':8,'name':'shared_face_maximal_induced_triangulations','all_pairs':66,'exact_match':True})
  # 9. Primitive rays, scales and independent true 152-column rates.
  PR=[];scales=[];rates=[];dots=[]
  for k,u in enumerate(U):
   x=tuple(dot([r[j] for r in rays],u) for j in range(6));D=lcm(*(z.denominator for z in x))
   ints=[int(z*D) for z in x];g=gcd(*ints);r=tuple(z//g for z in ints);scale=Q(g,D)
   save=fan['vertices'][k]
   assert r==tuple(save['primitive_physical_ray']) and scale==Q(save['slice_scale'])
   assert gcd(*r)==1 and r[5]>=0 and all(r[j]>=r[j+1] for j in range(5))
   ds=tuple(dot(a,r) for a in A);rate=Q(sum(abs(x) for x in ds),80)
   assert rate>0 and rate==Q(save['q_true_152'])
   active=save['active_witness'];assert len(active)==5 and len(rref([[1]*6]+[H[j] for j in active],6)[1])==6
   assert all(dot(H[j],u)==0 for j in active)
   PR.append(r);scales.append(scale);rates.append(rate);dots.append(ds)
  rels.append({'id':9,'name':'primitive_rays_and_rates','rays':12,'true_form_dots':1824,'rates':[str(q) for q in rates]})
  # 10. Full-cell and simplex sign constancy; independently recomputed slopes.
  for p in full:
   ids=p['vertices'];sg=[]
   for j in range(152):
    values=[dots[k][j] for k in ids];assert not(min(values)<0<max(values))
    sg.append(1 if max(values)>0 else -1 if min(values)<0 else 0)
   slope=tuple(Q(sum(sg[j]*A[j][a] for j in range(152)),80) for a in range(6))
   assert slope==tuple(Q(x) for x in p['linear_symbol_coefficient'])
   assert all(dot(slope,PR[k])==rates[k] for k in ids)
  rels.append({'id':10,'name':'all_true_forms_and_slopes','full_cells':12,'forms_per_cell':152,'no_old_q_reuse':True})
  # 11. Signed determinants, physical/coordinate relation, redundant exact total.
  volume=Q(0);keys=set();childcert=[]
  for j,t in enumerate(fan['simplices']):
   ids=t['vertices'];assert t['id']==j and len(ids)==6
   key=tuple(sorted(ids));assert key not in keys;keys.add(key)
   du=detcols([U[k] for k in ids]);dr=detcols([PR[k] for k in ids]);ratio=abs(du)
   assert ratio>0 and ratio==Q(t['relative_slice_volume']) and dr==Q(t['primitive_physical_determinant'])
   assert dr*prod(scales[k] for k in ids)==du
   assert all(not(min(dots[k][a] for k in ids)<0<max(dots[k][a] for k in ids)) for a in range(152))
   assert all(k in pats[t['pattern']]['vertices'] for k in ids)
   volume+=ratio;childcert.append({'id':j,'signed_u_determinant':str(du),'signed_primitive_determinant':str(dr),'relative_volume':str(ratio)})
  assert len(keys)==14 and volume==1
  for p in full:assert sum(Q(fan['simplices'][j]['relative_slice_volume']) for j in p['simplex_indices'])==Q(p['relative_slice_volume'])
  rels.append({'id':11,'name':'primitive_determinants_and_corrected_measure','simplices':14,'relative_total':'1','physical_volume':'sqrt(3)/120','not_psi_section':True})
  # 12. Negative global sign is an exact even-symbol symmetry, not an inserted factor.
  for r,q in zip(PR,rates):assert Q(sum(abs(dot(a,tuple(-x for x in r))) for a in A),80)==q
  rels.append({'id':12,'name':'global_negative_sign_scope','symbol_even_on_all_child_rays':True,'local_weight':'abs(physical determinant), no factor2','future_pairing':'opposite cones produce conjugate Fourier integrals; 2Re only after global coverage','full_fan_certified':False})
  assert len(rels)==12
  for f,h in clock['source_sha256'].items():assert sha(f)==h,('drift',f)
  elapsed=time.time()-clock['started_unix'];assert elapsed<900
  out={'status':'PASS_FIRST_PARENT_LOCAL_FAN_ONLY','relation_groups':rels,'active_system_ledger':audit,'independent_vertices':[[str(x) for x in u] for u in U],'shared_face_certificates':shared,'child_determinants':childcert,'pid':os.getpid(),'python':platform.python_version(),'runtime_seconds':time.monotonic()-tic,'audit_elapsed_seconds':elapsed,'source_sha256':sha(__file__),'author_manifest_sha256':sha(src/'manifest.json'),'kernel_calls':0,'author_writer_calls':0,'old_2160_checks':0,'full_fan_certified':False,'new_C6_bound':False}
  (P/'results.json').write_text(json.dumps(out,indent=2)+'\n')
  receipt.update(status='PASS',structure_relation_groups=12,elapsed_seconds=time.monotonic()-tic,finished_unix=time.time(),audit_elapsed_seconds=elapsed)
  print(json.dumps({k:out[k] for k in ['status','pid','runtime_seconds','audit_elapsed_seconds']}|{'active_systems':462,'relation_groups':12,'vertices':12,'full_cells':12,'simplices':14,'shared_faces':66,'relative_volume':'1'}))
 except BaseException as exc:
  receipt.update(status='FAILED_STOP_NO_RETRY',error=repr(exc),elapsed_seconds=time.monotonic()-tic,finished_unix=time.time());raise
 finally:(P/'execution_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
if __name__=='__main__':main()
