"""Saved fulljet g110 proof reuse, global075 anchor, no H or logarithm."""
from fractions import Fraction as F
from itertools import product
from cache_model import validate_record
def make_certificate(box,record,target,saved):
 target=F(target);validate_record(record,target);c=saved['certificate'];assert c['target']==str(target) and c['leaf_epsilon']=='1/128' and c['point_key']==record['pair_cache_key'][1]
 assert F(c['penalty_coefficient'])==1024 and F(record['values'][0][0])>0
 point=list(map(F,record['point']));lo,hi=box;gradient=[tuple(map(F,p)) for p in c['loggradient_directed']];delta=F(c['Delta_lower']);rows=[]
 for bits in product((0,1),repeat=4):
  v=[(hi[j] if bits[j] else lo[j])-point[j+1] for j in range(4)]
  rows.append((delta+sum(min(a*x,b*x) for x,(a,b) in zip(v,gradient))-1024*sum(x*x for x in v),bits))
 margin,bits=min(rows)
 return {'lower_bound_object':'G','component_function':'H_epsilon','leaf_epsilon':'1/128','target':str(target),'mechanism':'075_saved_global_fulljet_16','point_key':c['point_key'],'Delta_lower':c['Delta_lower'],'loggradient_directed':c['loggradient_directed'],'penalty_coefficient':'1024','minimum_log_margin':str(margin),'worst_corner_bits':list(bits),'checked_points':16,'passed':margin>0,'anchor_inside_box':all(a<=x<=b for a,x,b in zip(lo,point[1:],hi)),'saved_source_certificate_sha256':saved['certificate_sha256'],'scope':'single saved EPS128 anchor all16 corners;G only,not integral mean code'}
