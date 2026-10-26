"""Future separate registration/release required. Not executed by C053 prep."""
from fractions import Fraction as Q
from itertools import combinations
from prepare_adapter import frame,det
def exact_eligibility(numerators=(50,50,50,45,27,0),denominator=50):
 x=tuple(Q(v,denominator) for v in numerators);D=[v for v,beta in frame()]
 def dot(a,b):return sum(u*v for u,v in zip(a,b))
 r2=dot(x,x);line_min=None;pair_min=None;zero_pairs=[]
 for v in D:
  distance=r2-dot(x,v)**2/dot(v,v)
  assert distance>=0
  line_min=distance if line_min is None else min(line_min,distance)
 for i,j in combinations(range(152),2):
  v,w=D[i],D[j];a,b,c=dot(v,v),dot(v,w),dot(w,w);gram=a*c-b*b;assert gram>0
  xv,xw=dot(x,v),dot(x,w)
  distance=r2-(c*xv*xv-2*b*xv*xw+a*xw*xw)/gram;assert distance>=0
  pair_min=distance if pair_min is None else min(pair_min,distance)
  if distance==0:zero_pairs.append((i,j))
 return dict(direction_count=152,pair_count=11476,line_min_squared=str(line_min),pair_min_squared=str(pair_min),eligible_off_true_S2=not zero_pairs,zero_pairs=zero_pairs)
if __name__=='__main__':raise SystemExit('NOT RUN: separate structural-preflight registration/release required;never inherit body32 eligibility')
