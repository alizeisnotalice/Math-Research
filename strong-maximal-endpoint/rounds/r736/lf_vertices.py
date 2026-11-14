#!/usr/bin/env python3
"""Exact bounded-polytope clipping of the full responsible LF face."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import json,time,math,itertools,sys,hashlib
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor,as_completed
import sympy as sp
HERE=Path(__file__).resolve().parent
SOURCE=Path(sys.argv[1]) if len(sys.argv)>1 else HERE.parent/'r732'/'capacity.json'

def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def norm(a,b):
 vals=list(a)+[b];lcm=1
 for x in vals:lcm=math.lcm(lcm,x.denominator)
 ints=[int(x*lcm) for x in vals];g=0
 for x in ints:g=math.gcd(g,abs(x))
 if not any(ints[:-1]):
  assert ints[-1]>=0,'infeasible constant constraint';return None
 return tuple(x//g for x in ints)
def rank(rows,d):
 basis={}
 for row in rows:
  v=list(map(F,row))
  for p,b in basis.items():
   if v[p]:
    c=v[p];v=[x-c*y for x,y in zip(v,b)]
  piv=next((i for i,x in enumerate(v) if x),None)
  if piv is not None:
   c=v[piv];basis[piv]=[x/c for x in v]
 return len(basis)

def run(r):
 start=time.monotonic();mass=list(map(F,r['atom_masses']));vol=list(map(F,r['row_volumes']));masks=r['atom_masks'];f=list(map(F,r['source_f']));pi=list(map(F,r['pi']));q=F(r['target']);om=sum(mass);k=len(vol);a=len(mass);nv=k+1
 inc=[[int(s>>i&1) for s in masks] for i in range(k)];av=[sum(mass[j]*f[j]*inc[i][j] for j in range(a))/vol[i] for i in range(k)];frac=[j for j in range(a) if 0<f[j]<1]
 dens=[sum(pi[i]*inc[i][j]/vol[i] for i in range(k)) for j in range(a)]
 eta=dens[frac[0]] if frac else (max([dens[j] for j in range(a) if f[j]==0] or [F(0)])+min(dens[j] for j in range(a) if f[j]==1))/2
 point=pi+[eta];dmax=max(1/v for v in vol)
 eq=[[F(1)]*k+[F(0),F(1)],av+[F(0),q]]
 for j in frac:eq.append([inc[i][j]/vol[i] for i in range(k)]+[F(-1),F(0)])
 assert all(dot(row[:-1],point)==row[-1] for row in eq)
 rr,piv=sp.Matrix([[sp.Rational(x.numerator,x.denominator) for x in row] for row in eq]).rref();assert nv not in piv
 free=[j for j in range(nv) if j not in piv];d=len(free);x0=[F(0)]*nv;B=[[F(0)]*d for _ in range(nv)]
 for col,j in enumerate(free):B[j][col]=F(1)
 for row,j in enumerate(piv):
  x0[j]=F(int(rr[row,nv].p),int(rr[row,nv].q))
  for col,t in enumerate(free):B[j][col]=-F(int(rr[row,t].p),int(rr[row,t].q))
 known=tuple(point[j] for j in free)
 def transform(c,b):return norm([sum(c[i]*B[i][j] for i in range(nv)) for j in range(d)],b-dot(c,x0))
 original=[]
 for i in range(k):
  c=[F(0)]*nv;c[i]=F(-1);original.append((c,F(0)))
 c=[F(0)]*nv;c[k]=F(-1);original.append((c,F(0)));c=[F(0)]*nv;c[k]=F(1);original.append((c,dmax))
 for j in range(a):
  if f[j] in [0,1]:
   c=[inc[i][j]/vol[i] for i in range(k)]+[F(-1)]
   if f[j]==1:c=[-x for x in c]
   original.append((c,F(0)))
 for s in range(1,1<<k):original.append(([F(int(s>>i&1)) for i in range(k)]+[F(0)],sum(mass[j] for j,t in enumerate(masks) if s&t)/om))
 processed=[];seen=set()
 for j,t in enumerate(free):
  bound=dmax if t==k else F(1)
  for sign,b in [(-1,F(0)),(1,bound)]:
   c=[F(0)]*d;c[j]=F(sign);row=norm(c,b);processed.append(row);seen.add(row)
 limits=[dmax if j==k else F(1) for j in free]
 vs={tuple(z):sum(1<<i for i,row in enumerate(processed) if dot(row[:-1],z)==row[-1]) for z in itertools.product(*[(F(0),v) for v in limits])}
 steps=[];cache={};maxvertices=len(vs)
 def edge(av,bv):
  mask=av&bv
  if mask.bit_count()<d-1:return False
  if mask not in cache:
   cache[mask]=rank([processed[i][:-1] for i in range(len(processed)) if mask>>i&1],d)
   if len(cache)>200000:cache.clear();cache[mask]=rank([processed[i][:-1] for i in range(len(processed)) if mask>>i&1],d)
  return cache[mask]==d-1
 for c,b in original:
  row=transform(c,b)
  if row is None or row in seen:continue
  assert dot(row[:-1],known)<=row[-1],'known point excluded'
  seen.add(row);bit=1<<len(processed);neg=[];pos=[];new={}
  for v,active in vs.items():
   dv=dot(row[:-1],v)-row[-1]
   if dv>0:pos.append((v,active,dv))
   else:
    new[v]=active|(bit if dv==0 else 0)
    if dv<0:neg.append((v,active,dv))
  added=0
  for v,av0,dv in neg:
   for u,au,du in pos:
    if edge(av0,au):
     t=-dv/(du-dv);z=tuple(x+t*(y-x) for x,y in zip(v,u));assert dot(row[:-1],z)==row[-1]
     new[z]=new.get(z,0)|(av0&au)|bit;added+=1
  before=len(vs);processed.append(row);vs=new;assert vs,'empty LFface';maxvertices=max(maxvertices,len(vs));assert maxvertices<=20000,'vertex cap exceeded'
  steps.append(dict(constraint=len(processed)-1,before=before,outside=len(pos),edge_intersections=added,after=len(vs)))
 full=[];pi_unique={}
 for z in sorted(vs):
  x=[x0[i]+dot(B[i],z) for i in range(nv)];pv=x[:k];et=x[k]
  assert all(dot(c,x)<=b for c,b in original)
  assert all(dot(row[:-1],x)==row[-1] for row in eq)
  assert sum(pv)==1 and min(pv)>=0 and 0<=et<=dmax
  dd=[sum(pv[i]*inc[i][j]/vol[i] for i in range(k)) for j in range(a)];remain=q;cost=F(0)
  for j in sorted(range(a),key=lambda j:dd[j],reverse=True):
   take=min(remain,mass[j]*dd[j]);cost+=take/dd[j] if take else 0;remain-=take
   if remain==0:break
  assert remain==0 and cost==F(r['exact_NP_price'])
  full.append(dict(free_coordinates=[str(x) for x in z],pi=[str(x) for x in pv],threshold=str(et)))
  pi_unique[tuple(pv)]=None
 return dict(case=r['case'],target=r['target'],free_dimension=d,free_variables=free,affine_origin=[str(x) for x in x0],affine_basis=[[str(x) for x in row] for row in B],constraint_normals=[list(row) for row in processed],extended_vertex_count=len(full),unique_pi_count=len(pi_unique),maximum_intermediate_vertices=maxvertices,exact_LF_membership_and_NP_all=True,vertices=full,unique_pi=[[str(x) for x in pv] for pv in sorted(pi_unique)],clipping_counts=steps,wall_seconds=time.monotonic()-start)
if __name__=='__main__':
 start=time.monotonic();raw=SOURCE.read_bytes();rows=json.loads(raw)['results'];done=[]
 with ProcessPoolExecutor(max_workers=12) as pool:
  pending={pool.submit(run,r):r['case'] for r in rows}
  for job in as_completed(pending):
   result=job.result();done.append(result)
   print(json.dumps(dict(case=result['case'],free_dimension=result['free_dimension'],extended_vertices=result['extended_vertex_count'],unique_pi=result['unique_pi_count'],wall_seconds=result['wall_seconds'])),flush=True)
   data=dict(status='complete' if len(done)==len(rows) else 'partial',source_sha256=hashlib.sha256(raw).hexdigest(),workers=12,method='Exact rational bounded-polytope clipping, common-active rank edge enumeration',scope='LFface vertices; not capacities or central asymptotics',wall_seconds=time.monotonic()-start,results=sorted(done,key=lambda x:x['case']))
   (HERE/'lf_vertices.json').write_text(json.dumps(data,indent=2)+'\n')
 print(json.dumps(dict(status='complete',cases=len(done),wall_seconds=time.monotonic()-start,total_unique_pi=sum(r['unique_pi_count'] for r in done),max_intermediate_vertices=max(r['maximum_intermediate_vertices'] for r in done))))
