#!/usr/bin/env python3
"""Single frozen psi7-004 tangent template and 3 exact coefficient gates."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[k]='1'
from fractions import Fraction as F
from itertools import product,accumulate
from pathlib import Path
from flint import fmpz_mat
import json,hashlib,time,resource,sys,math,datetime,traceback
D=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def require(x,m):
 if not x:raise ValueError(m)
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def verify(lock):
 for rootkey,filekey in [('theory_dir','theory'),('r010_dir','r010')]:
  for n,h in lock[filekey].items():require(sha(Path(lock[rootkey])/n)==h,'source hash '+n)
def main():
 wall=time.perf_counter();cpu=time.process_time();when=datetime.datetime.now(datetime.timezone.utc).isoformat()
 lock=json.loads((D/'sourcehash.json').read_text());verify(lock)
 theory=Path(lock['theory_dir']);base=Path(lock['r010_dir']);spec=json.loads((theory/'prototype.json').read_text())
 require(spec['orbit_id']=='psi7-004','wrong prototype')
 body=[(1,)+z for z in product((-1,1),repeat=6)]; p=spec['primitive_p']; v=tuple(F(t,spec['psi_rate']) for t in p)
 values=[dot(z,p) for z in body];require(sum(map(abs,values))==344 and [i for i,t in enumerate(values) if t==0]==spec['zero_indices'],'psi or zero indices')
 basis=[body[i] for i in spec['basis_indices']];require([list(r) for r in basis]==spec['basis_rows'],'basis mismatch')
 require([sum(c*basis[j][i] for j,c in enumerate(spec['extra_in_basis'])) for i in range(7)]==spec['extra_row']==list(body[34]),'extra linear relation')
 g=[sum(((a>0)-(a<0))*z[i] for a,z in zip(values,body)) for i in range(7)]
 require(dot(g,v)==1 and all(dot(z,v)==0 for z in basis),'normalization')
 L=[[-x for x in g]]+[list(z) for z in basis]; det=int(fmpz_mat(L).det());require(det!=0,'singular L')
 q=[(1,0,0),(0,1,0),(0,0,1),(0,0,0)]
 triples=[]
 for es in spec['paired_triples_edges']:
  vec=[tuple(q[i-1][k]-q[j-1][k] for k in range(3)) for i,j in es]
  require(abs(int(fmpz_mat(vec).det()))==1,'nonunimodular')
  triples.append(vec)
 require(len(triples)==10 and len({tuple(x) for x in triples})==10,'wrong triple count')
 rays=json.loads((base/'r009-rays.json').read_text()); rr=[]
 for ray in rays:
  x=ray['xi_primitive']; norm=sum(abs(dot(z,x)) for z in body);require(norm>0,'bad ray')
  rr.append(tuple(F(t,norm) for t in x))
 ids=[i for i,x in enumerate(rr) if x==v];require(len(ids)==1,'vertex ray missing or duplicated'); target=ids[0]
 cells=[json.loads(s) for s in (base/'cells.jsonl').read_text().splitlines()];simplices=[s for c in cells for s in c['simplices']]
 require(len(cells)==135 and len(simplices)==212,'R010 shape changed')
 star=[s for s in simplices if target in s['vertices']]
 require(star and 8*len(star)<=1696,'star term cap')
 perms=[]
 for bits in product((0,1),repeat=3):
  a=list(range(7))
  for bit,(i,j) in zip(bits,[(1,2),(3,4),(5,6)]):
   if bit:a[i],a[j]=a[j],a[i]
  require(tuple(v[i] for i in a)==v,'not stabilizer');perms.append(a)
 require(len({tuple(x) for x in perms})==8,'stabilizer duplicated')
 gates=[]
 for gate,raw in enumerate(spec['gate_beta']):
  b=[F(t) for t in raw];require(len(b)==6 and max(map(abs,b))<1,'exposing condition')
  x=[F(g[j])+sum(b[i]*basis[i][j] for i in range(6)) for j in range(7)]
  require(dot(x,v)==1,'physical gate normalization')
  a=F(-1);c=(-b[3],b[4],b[5],F(0));q4=((F(1),F(0),F(0),F(0)),(F(0),F(1),F(0),F(0)),(F(0),F(0),F(1),F(0)),(F(0),F(0),F(0),F(0)))
  # c on zero-sum space; root coordinates are 3-dimensional.
  paired=[]
  for vec in triples:
   ds=[dot(c[:3],r) for r in vec];den=[4*a*a-d*d for d in ds];require(all(t>0 for t in den),'template pole')
   paired.append((4*a*a+ds[0]*ds[1]+ds[0]*ds[2]+ds[1]*ds[2])/math.prod(den))
  first=[a*a-z*z for z in b[:3]];require(all(z>0 for z in first),'template first-factor pole')
  template=F(32)*a**4/F(abs(det))/math.prod(first)*sum(paired,F(0))
  accum=F(0);termcount=0;min_den=None
  for perm in perms:
   for s in star:
    verts=s['vertices'];require(len(verts)==7 and len(set(verts))==7,'bad simplex')
    den=[]
    for i in verts:
     if i==target:continue
     u=[rr[i][j] for j in perm];d=dot(x,[v[k]-u[k] for k in range(7)]);require(d>0,'star denominator nonpositive');den.append(d);min_den=d if min_den is None else min(min_den,d)
    require(len(den)==6,'star denominator count');accum+=F(s['alpha'])/math.prod(den);termcount+=1
  require(termcount==8*len(star) and termcount<=1696,'term count')
  require(template==accum and template>0,'gate mismatch')
  gates.append({'id':gate,'beta':list(map(str,b)),'physical_x':list(map(str,x)),'x_dot_v':str(dot(x,v)),'template_value':str(template),'star_value':str(accum),'star_terms':termcount,'template_paired_terms':10,'min_star_denominator':str(min_den),'all_denominators_positive':True,'exact_equal':True})
 require(len(gates)==3,'gate count');verify(lock)
 return {'status':'PASS','orbit_id':spec['orbit_id'],'gradient_g':g,'L':L,'det_L':det,'psi_unit_vertex':list(map(str,v)),'basis_indices':spec['basis_indices'],'zero_indices':spec['zero_indices'],'ordered_ray_id':target,'incident_ordered_simplices':len(star),'stabilizer_size':len(perms),'three_dimensional_sign_cones':14,'simplices_3d':20,'paired_terms':10,'triples_vectors':triples,'gates':gates,'source_hashes_unchanged':True,'started_utc':when,'wall_seconds':time.perf_counter()-wall,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'python':sys.version,'executable':sys.executable,'worker_processes':1,'no_full_G_or_global_g':True}
if __name__=='__main__':
 try:
  result=main();print(json.dumps(result,indent=2));sys.exit(0)
 except Exception as exc:
  traceback.print_exc(file=sys.stderr);print(json.dumps({'status':'FAIL','error':repr(exc)}));sys.exit(1)
