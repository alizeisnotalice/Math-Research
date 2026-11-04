from fractions import Fraction as F
# A unit-cube tiling of R^3, phase 1/2 only in column (x,y)=(1,0) mod 2.
# q1-q2 corner i=(0,0,0), j=(1,1,0) has transverse contact [-1/2,1/2].
i=(-F(1,2),F(1,2)); j=i
pminus=(-F(1),F(0)); pplus=(F(0),F(1))
def meet(*intervals):
    lo=max(a for a,b in intervals); hi=min(b for a,b in intervals)
    return (lo,hi,max(F(0),hi-lo))
print('source_target=',meet(i,j))
print('pminus_piece=',meet(i,j,pminus))
print('pplus_piece=',meet(i,j,pplus))
print('neither_single_neighbor_covers=',all(meet(i,j,p)[2]<meet(i,j)[2] for p in [pminus,pplus]))
print('pieces_sum_exactly=',sum(meet(i,j,p)[2] for p in [pminus,pplus])==meet(i,j)[2])
