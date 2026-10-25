"""Exact rational sign-region diagnostic; no sampled estimate is an integral bound."""
from fractions import Fraction as Q
from itertools import combinations
import importlib.util, hashlib, json, argparse, time
from pathlib import Path
SOURCE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-finite-dim/cube-finite-dim-bounds/upper/u11_green.py')
SOURCE_SHA='fdfa422f2105710c84e7e43fdff1a4dcbef20dcc4e8a69d5a42116db3010b99d'
EPS=Q(1,10000000000)
def cross(u,v):return (u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0])
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def value(h,v):return h[0]+dot(h[1:],v)
def det(a,b,c):return dot(a,cross(b,c))
def load_source():
 assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==SOURCE_SHA
 spec=importlib.util.spec_from_file_location('ci4_source',SOURCE);g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g);return g
# (point, set of incident supporting planes) representation
PLANES=[(Q(0),Q(0),Q(0),Q(1)),(Q(0),Q(0),Q(1),Q(-1)),(Q(0),Q(1),Q(-1),Q(0)),(Q(1),Q(-1),Q(0),Q(0))]
BASE_PLANES=tuple(PLANES)
def initial():
 vs=[(Q(0),Q(0),Q(0)),(Q(1),Q(0),Q(0)),(Q(1),Q(1),Q(0)),(Q(1),Q(1),Q(1))]
 return [(v,frozenset(i for i,h in enumerate(PLANES) if value(h,v)==0)) for v in vs]
def edges(poly):
 for a,b in combinations(range(len(poly)),2):
  ids=poly[a][1]&poly[b][1]
  if len(ids)>=2 and any(cross(PLANES[i][1:],PLANES[j][1:])!=(0,0,0) for i,j in combinations(ids,2)):
   yield a,b

def split(poly,h,idx):
 vals=[value(h,v) for v,active in poly]
 if min(vals)>=0 or max(vals)<=0:return [poly]
 cut=[]
 for i,j in edges(poly):
  a,b=vals[i],vals[j]
  if a*b<0:
   va,aa=poly[i];vb,ab=poly[j];t=a/(a-b)
   v=tuple(x+t*(y-x) for x,y in zip(va,vb))
   cut.append((v,(aa&ab)|{idx}))
 result=[]
 for sg in (-1,1):
  pp=[(v,aa|({idx} if a==0 else set())) for (v,aa),a in zip(poly,vals) if sg*a>=0]+cut
  result.append(list({v:(v,frozenset(a)) for v,a in pp}.values()))
 return result

def planar_hull(points):
 # Exact Andrew hull. A coordinate projection is injective on a plane after
 # dropping any coordinate with nonzero plane normal.
 pts=sorted(set(points))
 def turn(a,b,c):return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
 lower=[]
 for p in pts:
  while len(lower)>=2 and turn(lower[-2],lower[-1],p)<=0:lower.pop()
  lower.append(p)
 upper=[]
 for p in reversed(pts):
  while len(upper)>=2 and turn(upper[-2],upper[-1],p)<=0:upper.pop()
  upper.append(p)
 return lower[:-1]+upper[:-1]

def volume(poly):
 vs=[v for v,a in poly];c=tuple(sum(v[j] for v in vs)/len(vs) for j in range(3))
 total=Q(0);seen=set()
 for idx in set().union(*(a for v,a in poly)):
  face=[v for v,a in poly if idx in a]
  if len(face)<3:continue
  h=PLANES[idx];drop=next(j for j in range(3) if h[j+1]);keep=[j for j in range(3) if j!=drop]
  lookup={tuple(v[j] for j in keep):v for v in face};hull=planar_hull(lookup)
  if len(hull)<3:continue
  verts=[lookup[p] for p in hull];key=frozenset(verts)
  assert key not in seen,'duplicate facet';seen.add(key)
  vals=[value(h,v) for v in vs];assert min(vals)>=0 or max(vals)<=0,'not supporting'
  for k in range(1,len(verts)-1):
   a,b,d=verts[0],verts[k],verts[k+1]
   total+=abs(det(tuple(a[j]-c[j] for j in range(3)),tuple(b[j]-c[j] for j in range(3)),tuple(d[j]-c[j] for j in range(3))))/6
 assert total>0
 return total

def build():
 PLANES[:]=BASE_PLANES
 g=load_source();root=initial();tv=[(1,)+v for v,a in root]
 hp=[]
 for r in g.RAYS:
  if next(x for x in r if x)<0:continue
  vals=[dot(r,v) for v in tv]
  if min(vals)<0<max(vals):hp.append(tuple(map(Q,r)))
 assert len(hp)==30
 # Remove only the simplex sum deficits<epsilon; it lies in radius epsilon cap.
 cap=(Q(3)-EPS,Q(-1),Q(-1),Q(-1));idx=len(PLANES);PLANES.append(cap)
 ps=[p for p in split(root,cap,idx) if value(cap,tuple(sum(v[j] for v,a in p)/len(p) for j in range(3)))>=0]
 counts=[]
 for h in hp:
  idx=len(PLANES);PLANES.append(h);ps=[q for p in ps for q in split(p,h,idx)]
  counts.append(len(ps));assert len(ps)<=4096
 return g,ps,hp,counts

def main(save):
 start=time.time();g,ps,hp,counts=build()
 vols=[volume(p) for p in ps];expected=Q(1,6)-EPS**3/36
 # ordered cap is one sixth of octant simplex of volume eps^3/6.
 assert sum(vols)==expected,(sum(vols),expected)
 worst=0;zero_cross=0;bad=[]
 q=g.rates()
 for ix,p in enumerate(ps):
  bounds=[(min(dot((Q(1),)+v,r)/qq for v,a in p),max(dot((Q(1),)+v,r)/qq for v,a in p)) for r,qq in zip(g.RAYS,q)]
  assert all(a>=0 or b<=0 for a,b in bounds)
  # Potential unresolved cross-sign denominator at a node-zero intersection.
  for ids,d in g.CONES:
   z=[bounds[i] for i in ids]
   for a,b in combinations(z,2):
    if a[0]==0<a[1] and b[0]<0==b[1]:
     zero_cross+=1
     if len(bad)<12:bad.append({'region':ix,'cone':list(ids),'pair':[list(map(str,a)),list(map(str,b))]})
    if b[0]==0<b[1] and a[0]<0==a[1]:
     zero_cross+=1
     if len(bad)<12:bad.append({'region':ix,'cone':list(ids),'pair':[list(map(str,a)),list(map(str,b))]})
  worst=max(worst,len(p))
 out={'status':'geometry certified; integral not enclosed','source_sha256':SOURCE_SHA,'epsilon':str(EPS),'region_count':len(ps),'counts_after_hyperplanes':counts,'node_hyperplanes':[[int(z) for z in h] for h in hp],'total_volume':str(sum(vols)),'expected_volume':str(expected),'max_vertices':worst,'cross_sign_zero_endpoint_pairs':zero_cross,'examples':bad,'seconds':time.time()-start}
 if save:Path(save).write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--save');args=a.parse_args();main(args.save)
