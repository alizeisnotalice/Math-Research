import ci4_partition as c
from itertools import combinations
from fractions import Fraction as Q
import json,time
from pathlib import Path
start=time.time();g,ps,hp,counts=c.build();rates=g.rates();regbad=[];denoms=0;zero=0;false=0;examples=[]
for ix,p in enumerate(ps):
 vs=[(Q(1),)+v for v,a in p];nodes=[[c.dot(x,r)/q for x in vs] for r,q in zip(g.RAYS,rates)]
 signs=[1 if min(a)>=0 else -1 if max(a)<=0 else 0 for a in nodes]
 assert all(signs)
 bad=False
 for ids,d in g.CONES:
  for a,b in combinations(ids,2):
   if signs[a]*signs[b]!=-1:continue
   denoms+=1
   exact=[signs[a]*nodes[a][j]+signs[b]*nodes[b][j] for j in range(len(vs))]
   assert min(exact)>=0
   marginalzero=min(map(abs,nodes[a]))==min(map(abs,nodes[b]))==0
   if min(exact)==0:
    zero+=1;bad=True
    if len(examples)<8:
     j=exact.index(0);examples.append({'region':ix,'cone':list(ids),'rays':[a,b],'joint_zero_vertex':list(map(str,vs[j])),'coincident_zeros':[str(nodes[a][j]),str(nodes[b][j])]})
   elif marginalzero:false+=1
 if bad:regbad.append(ix)
vols=[c.volume(p) for p in ps]
out={'status':'one correlation revision diagnostic; not an integral enclosure','regions':len(ps),'cross_denominator_pairs':denoms,'true_joint_zero_pairs':zero,'false_marginal_zero_pairs_resolved':false,'unresolved_regions':len(regbad),'unresolved_volume':str(sum(vols[i] for i in regbad)),'resolved_regions':len(ps)-len(regbad),'resolved_volume':str(sum(vols[i] for i in range(len(ps)) if i not in regbad)),'examples':examples,'seconds':time.time()-start}
x=(Q(1),Q(1),Q(1),Q(0))
rr=g.confluent_residues(x);C=g.green_allnodes_pi3_coefficient(x)
_,_,point_nodes=g.point_nodes(x)
assert C==Q(17,640) and not any(b==0 for b,k in rr) and point_nodes[14]==point_nodes[61]==0
out['finite_boundary_point']={'x':list(map(str,x)),'distance_to_body_lines_squared':str(g.body_line_distance_squared(x)),'G_times_pi3':str(C),'ray14':g.RAYS[14],'ray61':g.RAYS[61],'node14':str(point_nodes[14]),'node61':str(point_nodes[61]),'full_sum_zero_poles':0}
Path(__file__).with_name('ci4_correlation_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
