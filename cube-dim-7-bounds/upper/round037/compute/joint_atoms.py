"""Exact inactive-coordinate radial slab atomization for A/B rational cubes."""
from fractions import Fraction as Q
from collections import defaultdict
class GeometryCensored(Exception):pass

def atoms(cubes,n,a,limit=8_000_000):
 k=n-2;N=len(cubes);E=a+4
 cuts=[]
 for coord in range(2):
  values=set()
  for z,r in cubes:
   for kappa in (1,E):values.add(z[coord]-kappa*r);values.add(z[coord]+kappa*r)
  cuts.append(sorted(values))
 radii=sorted({Q(0)}|{r for _,r in cubes}|{E*r for _,r in cubes})
 potential=(len(cuts[0])-1)*(len(cuts[1])-1)*(len(radii)-1)
 if potential>=limit:raise GeometryCensored('elementary_cells_limit preflight '+str(potential))
 def masks(coord,lo,hi):
  mid=(lo+hi)/2; qm=0;em=0
  for i,(z,r) in enumerate(cubes):
   if abs(mid-z[coord])<r:qm|=1<<i
   if abs(mid-z[coord])<E*r:em|=1<<i
  return qm,em
 xx=[(cuts[0][i+1]-cuts[0][i],*masks(0,cuts[0][i],cuts[0][i+1])) for i in range(len(cuts[0])-1)]
 yy=[(cuts[1][i+1]-cuts[1][i],*masks(1,cuts[1][i],cuts[1][i+1])) for i in range(len(cuts[1])-1)]
 slabs=[]
 for lo,hi in zip(radii,radii[1:]):
  rho=(lo+hi)/2;qm=0;em=0
  for i,(_,r) in enumerate(cubes):
   if rho<r:qm|=1<<i
   if rho<E*r:em|=1<<i
  slabs.append((qm,em,(2*hi)**k-(2*lo)**k))
 out=defaultdict(Q)
 for dx,xq,xe in xx:
  if not xq:continue
  for dy,yq,ye in yy:
   active=xq&yq
   if not active:continue
   em=xe&ye; area=dx*dy
   for rq,re,vol in slabs:
    qmask=active&rq
    if qmask:out[(qmask,em&re)]+=area*vol
 return dict(out),potential

def increment(cells,index,selected):return sum((vol for (qm,_),vol in cells.items() if qm>>index&1 and not qm&selected),Q())

def initial_flow(cells,cubes,n,a,order,tau=None):
 selected=0;ps={};vs={};h=[(2*a*r)**n for _,r in cubes]
 for i in order:
  centre=cubes[i][0]
  current=sum((p for j,p in ps.items() if all(abs(x-y)<(a+3)*cubes[j][1] for x,y in zip(centre,cubes[j][0]))),Q())
  if tau is None or current<=tau:
   v=increment(cells,i,selected);selected|=1<<i;vs[i]=v;ps[i]=v/h[i]
 total=sum(cells.values(),Q());d={}
 for (qm,em),vol in cells.items():
  if qm&selected:continue
  mask=em&selected
  if not mask:raise AssertionError('uncovered_rejected_atom')
  d[mask]=d.get(mask,Q())+vol
 charges={i:vs[i] for i in vs}
 for mask,vol in d.items():
  ptotal=sum((p for j,p in ps.items() if mask>>j&1),Q())
  if tau is not None and not ptotal>tau:raise AssertionError('rejected_atom_payment_not_greater_than_tau')
  if ptotal<=0:raise AssertionError('zero_rejected_atom_payment')
  for j,p in ps.items():
   if mask>>j&1:charges[j]+=vol*p/ptotal
 if sum(charges.values(),Q())!=total:raise AssertionError('mass_conservation')
 return selected,ps,[charges.get(i,Q())/h[i] for i in range(len(cubes))],total,len(d)
