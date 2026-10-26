"""Target-parameterized H128 fullQ certificate. No H on import."""
from fractions import Fraction as F
from itertools import product
from flint import ctx
from binding import api
from cache_model import validate_record
def make_certificate(box,record,target):
 target=F(target);validate_record(record,target);epsilon=F(record['epsilon']);penalty=1/(16*epsilon**2);assert penalty==1024
 ctx.prec=record['bits'];hl,hu=map(F,record['values'][0]);assert target<hl<=hu
 gradients=[]
 for pair in record['values'][1:]:
  ratio=api.unpack(pair)/api.unpack(record['values'][0]);gl,gu=map(F,api.pack(ratio));a,b=map(F,pair)
  assert gl<=min(a/hl,a/hu,b/hl,b/hu)<=max(a/hl,a/hu,b/hl,b/hu)<=gu;gradients.append((gl,gu))
 delta=api.endpoint((api.ab(hl)/api.ab(target)).log());point=list(map(F,record['point']));lo,hi=box
 assert point[0]==1 and all(a<=p<=b for a,p,b in zip(lo,point[1:],hi))
 rows=[]
 for bits in product((0,1),repeat=4):
  v=[(hi[j] if bits[j] else lo[j])-point[j+1] for j in range(4)]
  rows.append((delta+sum(min(a*x,b*x) for x,(a,b) in zip(v,gradients))-penalty*sum(x*x for x in v),bits))
 margin,bits=min(rows)
 return {'lower_bound_object':'G','component_function':'H_epsilon','leaf_epsilon':str(epsilon),'target':str(target),'mechanism':'075_full_frame_16','quadratic':'I4/(16epsilon^2)',
  'penalty_coefficient':str(penalty),'checked_points':16,'point_key':record['pair_cache_key'][1],'Delta_lower':str(delta),'loggradient_directed':[[str(a),str(b)] for a,b in gradients],
  'minimum_log_margin':str(margin),'worst_corner_bits':list(bits),'passed':margin>0,'inference':'same-epsilon H bound then G>=H_epsilon;G homogeneity only'}
