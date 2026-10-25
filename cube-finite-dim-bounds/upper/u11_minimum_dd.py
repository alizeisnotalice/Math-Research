"""Stable DD3(sign) and affine-cell bounds; shared g/I interface.
Input vertex nodes must be normalized b=(x.r)/rate, exact Fractions.
Bounds are exact rational coefficients C with G=C/pi^3.
None lower/upper explicitly means unbounded at an individual cone; never zero.
"""
from fractions import Fraction as Q
from math import prod
import importlib.util
from pathlib import Path
SOURCE=Path(__file__).with_name('u11_green.py')
_sp=importlib.util.spec_from_file_location('u11_minimum_source',SOURCE)
g=importlib.util.module_from_spec(_sp);_sp.loader.exec_module(g)

QBODY=g.rates();ALPHA=tuple(Q(det)/prod(QBODY[i] for i in ids) for ids,det in g.CONES)

def dd_sign_stable(nodes):
 z=list(map(Q,nodes))
 if any(b==0 for b in z):raise ValueError('zero node: full-sum/confluent handling required')
 p=[b for b in z if b>0];n=[-b for b in z if b<0]
 if not p or not n:return Q(0)
 if len(p)==1:return Q(2)/prod(p[0]+a for a in n)
 if len(n)==1:return Q(2)/prod(n[0]+a for a in p)
 return -Q(2)*sum(p+n)/prod(a+b for a in p for b in n)

def normalized_vertex_nodes(vertices):
 return [[g.dot(tuple(map(Q,x)),r)/q for r,q in zip(g.RAYS,QBODY)] for x in vertices]

def green_coefficient_interval(node_values,signs=None,details=False):
 """Exact one-sided intervals for G*pi^3 on a sign-fixed affine polytope.
 Shared vertex rows preserve p+n denominator correlation. Signs describe INTERIOR.
 q=0 vertex values permitted as boundary limits; true cross0 flags poles.
 """
 vv=[tuple(map(Q,row)) for row in node_values]
 if not vv or any(len(row)!=len(g.RAYS) for row in vv):raise ValueError('need Nx240 nodes')
 if signs is None:
  sums=[sum(row[i] for row in vv) for i in range(len(g.RAYS))]
  signs=[1 if x>0 else -1 if x<0 else 0 for x in sums]
 if any(sg not in (-1,1) for sg in signs):raise ValueError('full-dimensional sign region required')
 if any(signs[i]*row[i]<0 for row in vv for i in range(len(g.RAYS))):raise ValueError('cell crosses a node hyperplane')
 lo=Q(0);hi=Q(0);lbad=[];ubad=[]
 cache={}
 def cross_bounds(i,j):
  key=(i,j)
  if key not in cache:
   z=[row[i]-row[j] for row in vv];a,b=min(z),max(z)
   if a<0:raise ValueError('cross sum inconsistent with signs')
   cache[key]=(a,b)
  return cache[key]
 for ci,((ids,det),alpha) in enumerate(zip(g.CONES,ALPHA)):
  pos=[i for i in ids if signs[i]>0];neg=[i for i in ids if signs[i]<0]
  if not pos or not neg:continue
  pairs=[(i,j) for i in pos for j in neg]
  if len(pos)==1 or len(neg)==1:
   ds=[cross_bounds(i,j) for i,j in pairs]
   if any(a==0 for a,b in ds):
    lbad.append({'cone_index':ci,'ids':ids,'zero_pairs':[ij for ij,(a,b) in zip(pairs,ds) if a==0]})
   else:lo-=alpha/(16*prod(a for a,b in ds))
   if any(b==0 for a,b in ds):raise ValueError('full dimensional positive cross identically zero')
   hi-=alpha/(16*prod(b for a,b in ds))
  else:
   ds=[cross_bounds(i,j) for i,j in pairs]
   ns=[sum(row[i] for i in pos)-sum(row[j] for j in neg) for row in vv]
   nlo,nhi=min(ns),max(ns)
   if any(b==0 for a,b in ds):raise ValueError('cross identically zero')
   lo+=alpha*nlo/(16*prod(b for a,b in ds))
   if any(a==0 for a,b in ds):
    ubad.append({'cone_index':ci,'ids':ids,'zero_pairs':[ij for ij,(a,b) in zip(pairs,ds) if a==0]})
   else:hi+=alpha*nhi/(16*prod(a for a,b in ds))
 return {'lower_pi3':None if lbad else lo,'upper_pi3':None if ubad else hi,
         'lower_pole_count':len(lbad),'upper_pole_count':len(ubad),
         'first_lower_pole':lbad[0] if lbad else None,'first_upper_pole':ubad[0] if ubad else None}
