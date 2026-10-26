"""Prepared H128 closed-box certificate for a full-G lower bound.

Functions are definitions only; preparation does not evaluate H or logs.
"""
from fractions import Fraction as F
from itertools import product
import epsilon_kernel as api
from flint import ctx
EPSILON=F(1,128)
TARGET=F(21,200000)
PENALTY=F(1)/(16*EPSILON*EPSILON)
assert PENALTY==1024
def make_certificate(box,record):
 assert F(record['epsilon'])==EPSILON,'H64 values cannot be interpreted as H128'
 ctx.prec=record['bits'];hl,hu=map(F,record['values'][0]);assert TARGET<hl<=hu
 gradients=[]
 for pair in record['values'][1:]:
  ratio=api.unpack(pair)/api.unpack(record['values'][0]);gl,gu=map(F,api.pack(ratio));a,b=map(F,pair)
  assert gl<=min(a/hl,a/hu,b/hl,b/hu)<=max(a/hl,a/hu,b/hl,b/hu)<=gu
  gradients.append((gl,gu))
 delta=api.endpoint((api.ab(hl)/api.ab(TARGET)).log());point=list(map(F,record['point']));lo,hi=box
 assert point[0]==1 and all(a<=p<=b for a,p,b in zip(lo,point[1:],hi))
 rows=[]
 for bits in product((0,1),repeat=4):
  v=[(hi[j] if bits[j] else lo[j])-point[j+1] for j in range(4)]
  rows.append((delta+sum(min(a*x,b*x) for x,(a,b) in zip(v,gradients))-PENALTY*sum(x*x for x in v),bits))
 margin,bits=min(rows)
 return {'lower_bound_object':'G','component_function':'H_epsilon','leaf_epsilon':'1/128','target':'21/200000',
  'mechanism':'075_full_frame_16','quadratic':'I4/(16epsilon^2)','penalty_coefficient':'1024','checked_points':16,
  'point_key':api.pointkey(record['point']),'Delta_lower':str(delta),'loggradient_directed':[[str(a),str(b)] for a,b in gradients],
  'minimum_log_margin':str(margin),'worst_corner_bits':list(bits),'passed':margin>0,
  'inference':'H_1/128>target on this closed box implies G>target on this box;no H homogeneity claim'}
