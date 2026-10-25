#!/usr/bin/env python3
"""Reconstruct point geometry and union directly from saved rational witnesses."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from pathlib import Path
from math import prod
import json,hashlib
HERE=Path(__file__).resolve().parent
def verify(path):
 data=json.loads(path.read_text());a=list(map(Q,data['spacing_rationals']));eps=Q(data['epsilon']);shift=data['half_shift'];rows=data['witnesses'];boxes=[];ineq=0
 for row in rows:
  atoms={};r=Q(row['radius']);box=tuple(map(Q,row['box']))
  for orbit,key,weight in [(0,'lengths_primary',Q(1)),(1,'lengths_satellite',eps)]:
   ls=row[key]
   if ls is None:continue
   coords=[]
   for d,n in enumerate(ls):
    s=shift[d] if orbit else 0
    # Sort a generous coordinate cloud by distance to an interior receiver;
    # gives nearest consecutive set without importing producer index formula.
    vals=sorted(range(-5,6),key=lambda k:abs(Q(k)+Q(s,2)-Q(1,4)))[:n]
    coords.append(sorted(vals))
   for inds in product(*coords):
    pt=tuple((Q(k)+Q(shift[d],2)*orbit)*a[d] for d,k in enumerate(inds))
    assert pt not in atoms;atoms[pt]=weight
  assert sum(atoms.values())==Q(row['mass'])
  assert (2*r)**5<Q(row['mass'])
  for p in atoms:
   for d in range(5):
    assert abs(p[d]-box[2*d])<=r and abs(p[d]-box[2*d+1])<=r
    ineq+=2
  boxes.append(box)
 axes=[]
 for d in range(5):
  cuts=sorted(set(x for b in boxes for x in b[2*d:2*d+2]));axis=[]
  for lo,hi in zip(cuts,cuts[1:]):
   ids=tuple(i for i,b in enumerate(boxes) if b[2*d]<=lo and hi<=b[2*d+1])
   if ids:axis.append((hi-lo,frozenset(ids)))
  axes.append(axis)
 @lru_cache(None)
 def vol(d,ids):
  if d==5:return Q(1)
  return sum((length*vol(d+1,ids&other) for length,other in axes[d] if ids&other),Q(0))
 volume=vol(0,frozenset(range(len(boxes))))
 ratio=32*volume/(1+eps);assert ratio==Q(data['ratio_exact'])
 C=json.loads((HERE.parent/'lower_round003/a2_equals_a3_results.json').read_text());assert ratio<Q(C['ratio_interval'][0])
 return dict(file=path.name,ratio_matches=True,below_complete_C=True,checked_point_endpoint_inequalities=ineq,atom_geometry_reconstructed_without_producer_import=True,source_result_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
out={'status':'saved_rational_witnesses_verified','results':[verify(HERE/f'exact_{kind}_results.json') for kind in ('fixed','codilation')]}
(HERE/'saved_verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
