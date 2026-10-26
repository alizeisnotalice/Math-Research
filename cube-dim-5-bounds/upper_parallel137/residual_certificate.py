"""New-jet Delta and physical gradient, then129 anchored Jensen curvature."""
from fractions import Fraction as F
from itertools import product
def make(box,record,target,task,legacy):
 legacy.validate_record(record,target);assert F(target)==F(11,100000) and record['epsilon']=='1/128'
 api=legacy.worker.api;api.ctx.prec=record['bits'];hl,hu=map(F,record['values'][0]);assert F(target)<hl<=hu
 proof=task['H64_g105_ancestor'];ancestor=proof['node'];node=task['node'];assert ancestor[0]==node[0] and node[1].startswith(ancestor[1]) and F(proof['lower'])==F(21,200000)
 point=list(map(F,record['point']));lo,hi=box;assert point[0]==1 and all(a<x<b for a,x,b in zip(lo,point[1:],hi))
 q=1024-768*F(proof['lower'])/hu;assert 256<=q<1024
 gradients=[]
 for pair in record['values'][1:]:
  gl,gu=map(F,api.pack(api.unpack(pair)/api.unpack(record['values'][0])));a,b=map(F,pair);assert gl<=min(a/hl,a/hu,b/hl,b/hu)<=max(a/hl,a/hu,b/hl,b/hu)<=gu;gradients.append((gl,gu))
 delta=api.endpoint((api.ab(hl)/api.ab(target)).log());rows=[]
 for bits in product((0,1),repeat=4):
  v=[(hi[j] if bits[j] else lo[j])-point[j+1] for j in range(4)];rows.append((delta+sum(min(a*x,b*x) for x,(a,b) in zip(v,gradients))-q*sum(x*x for x in v),bits))
 margin,bits=min(rows)
 return {'lower_bound_object':'G','component_function':'H_1/128','leaf_epsilon':'1/128','target':str(F(target)),'mechanism':'129_anchor_Jensen_g105_U128','penalty_coefficient':str(q),'H64_lower':proof['lower'],'H64_lower_source_node':ancestor,'H64_g105_ancestor_proof':proof,'H128_upper':str(hu),'checked_points':16,'point_key':record['pair_cache_key'][1],'Delta_lower':str(delta),'loggradient_directed':[[str(a),str(b)] for a,b in gradients],'minimum_log_margin':str(margin),'worst_corner_bits':list(bits),'passed':margin>0,'inference':'single fixed anchor all16 closed corners;own EPS128 physical jet;129 H64 mass from containing certified ancestor;G>=H128'}
