#!/usr/bin/env python3
"""Exact rational finite-instance checks of Lipschitz small-box Gram inequalities."""
from fractions import Fraction as F
import json, random, time
from functools import reduce
from operator import mul

SEED=20261029

def check(n,L,r,s,y):
 assert len(y)==n and r>0 and s>0 and L>=0
 lam=1/(1+L/2); d=max(map(abs,y)); assert abs(r-s)<=L*d
 widths=[max(F(0), min(lam*r/2, yi+lam*s/2)-max(-lam*r/2,yi-lam*s/2)) for yi in y]
 overlap=reduce(mul,widths,F(1)); G=overlap/((lam*r)**n*(lam*s)**n)
 if G>0:assert d<=min(r,s)
 for z in (r,s):
  rhs=lam**(-n)*z**(-n)*max(F(0),1-d/z)
  assert G<=rhs, (n,L,r,s,y,G,rhs)
  left=z**(-n)*(d<=z)-lam**n*G
  right=d/z**(n+1)*(d<=z)
  assert left>=right,(n,L,r,s,y,left,right)
 return G,d

def main():
 st=time.time();rng=random.Random(SEED);counts={}; total=0; positive=0; zero=0; d0=0
 for n in range(2,7):
  for L in map(F,('0','1/100','1/10','1','10')):
   key=f'n={n},L={L}';yes=0;tries=0
   for i in range(400):
    r=F(rng.randint(1,500),100); d=r*F(rng.randint(0,200),100)
    while True:
     u=F(rng.randint(-100,100),100);s=r+u*L*d;tries+=1
     if s>0:break
    y=[d*F(rng.randint(-100,100),100) for _ in range(n)];y[rng.randrange(n)]=d
    G,dist=check(n,L,r,s,y);yes+=bool(G);positive+=bool(G);zero+=not bool(G);d0+=dist==0;total+=1
   counts[key]={'samples':400,'positive_overlap':yes,'radius_draws':tries}
 # Exact boundary tests: coincident centers, common-radius touching, and L maximal slope.
 boundaries=0
 for n in range(2,7):
  for L in map(F,('0','1/100','1/10','1','10')):
   r=F(3,2)
   for d in (F(0),r,F(1,2),F(3)):
    for s in {r, r+L*d}:
     check(n,L,r,s,[d]+[F(0)]*(n-1));boundaries+=1
 out={'seed':SEED,'random_checks':total,'boundary_checks':boundaries,'positive_overlap':positive,'zero_overlap':zero,'coincident_centers':d0,'dimensions':list(range(2,7)),'L_values':['0','1/100','1/10','1','10'],'all_pass':True,'elapsed_seconds':time.time()-st,'per_group':counts,'scope':'finite exact rational instances; not a proof of universal inequalities'}
 print({k:v for k,v in out.items() if k!='per_group'})
if __name__=='__main__':main()
