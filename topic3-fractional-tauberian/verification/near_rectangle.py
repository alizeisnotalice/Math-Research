#!/usr/bin/env python3
"""Exact common-interval capacity certificate and anisotropic obstruction."""
from fractions import Fraction as F
from itertools import permutations
from collections import Counter
from pathlib import Path
import random,json,platform
from near_square import flat_R

def T(v):return 4*(flat_R(F(1,2),[z/2 for z in v])+sum(v)-2)

def certificate():
 rows=[];hist=Counter()
 for order in permutations(range(4)):
  V=[]
  for k in range(5):
   v=[F(0)]*4
   for j in range(4):v[order[j]]=F(-1 if j<k else 1)
   V.append(v)
  def q(v):return 8-T(v)
  coeff=[]
  for i in range(5):
   for j in range(i,5):
    b=q(V[i]) if i==j else 2*q([(a+b)/2 for a,b in zip(V[i],V[j])])-(q(V[i])+q(V[j]))/2
    assert b>=0;coeff.append(str(b));hist[str(b)]+=1
  rows.append({'order':order,'bernstein_upper_triangle':coeff})
 return rows,dict(sorted(hist.items(),key=lambda z:F(z[0])))

def envelope_integral(p):
 # Symmetry reduces to [0,1/2]; both pieces are affine after cancellation.
 a=(1-p)/2
 def affine(A,B,l,r):return A*(r*r-l*l)/2+B*(r-l)
 return 2*(affine(8*(1+p),2*p*p,F(0),a)+affine(-8*(1-p),8*(1-p/2)**2,a,F(1,2)))

def phi(x,u):
 p=F(1,2);kx=F(3,2);ky=F(1,4);s=x-u+p/2;q=1-x-u-p/2
 return u*u+u+x*(1-x)+p*(x-F(1,2))+p*p/4+kx*ky-max(kx-s,0)*ky-kx*max(ky-q,0)-max(kx-s,0)*max(ky-q,0)

def exact_example():
 # Exhaust all u contacts and all affine x branches, then integrate their envelope.
 p=F(1,2);kx=F(3,2);ky=F(1,4);pieces=[]
 for L,R in [(F(0),F(1,4)),(F(1,4),F(3,4)),(F(3,4),F(1))]:
  mid=(L+R)/2;lo=(-F(1),p/2) if mid<F(3,4) else (F(1),-1-p/2);hi=(F(1),p/2) if mid<F(1,4) else (-F(1),1-p/2)
  cand=[lo,hi,(F(1),p/2-kx),(-F(1),1-p/2-ky)];cuts={L,R}
  for a,b in cand:
   for A,B in [(a-lo[0],b-lo[1]),(a-hi[0],b-hi[1]),(a-1,b+kx-p/2),(a+1,b+ky-1+p/2)]:
    if A and L<-B/A<R:cuts.add(-B/A)
  cuts=sorted(cuts)
  for l,r in zip(cuts,cuts[1:]):
   mid=(l+r)/2;lines=[]
   for a,b in cand:
    if not lo[0]*mid+lo[1]<=a*mid+b<=hi[0]*mid+hi[1]:continue
    z0,z1=phi(l,a*l+b),phi(r,a*r+b);A=(z1-z0)/(r-l);B=z0-A*l
    assert phi(mid,a*mid+b)==A*mid+B;lines.append((A,B))
   splits={l,r}
   for A,B in lines:
    for C,D in lines:
     if A!=C and l<(D-B)/(A-C)<r:splits.add((D-B)/(A-C))
   splits=sorted(splits)
   for a,b in zip(splits,splits[1:]):
    A,B=max(lines,key=lambda z:z[0]*(a+b)/2+z[1])
    if pieces and pieces[-1][1]==a and pieces[-1][2:]==[A,B]:pieces[-1][1]=b
    else:pieces.append([a,b,A,B])
 cubic=4*sum(A*(r*r-l*l)/2+B*(r-l) for l,r,A,B in pieces);assert cubic==F(29,16)
 pair=[];empty_pair=[]
 for x in [F(1,10),F(9,10)]:
  lo=max(p/2-x,-p/2-1+x);hi=min(p/2+x,-p/2+1-x);us=[lo,hi,x+p/2-kx,1-x-p/2-ky]
  got=max(phi(x,u) for u in us if lo<=u<=hi);assert got==F(17,40);pair.append(got);empty_pair.append(hi*hi+hi+x*(1-x)+p*(x-F(1,2))+p*p/4)
 assert empty_pair==[F(17,40),F(9,40)] and 2*sum(pair)==F(17,10) and 2*sum(empty_pair)==F(13,10)
 return {'rho':'1/2','hole_rectangle':['3/2','1/4'],'G_pieces':[[str(z) for z in r] for r in pieces],'actual_cubic':str(cubic),'pair_G':[str(v) for v in pair],'empty_pair_G':[str(v) for v in empty_pair],'pair_average':str(2*sum(pair)),'empty_pair_average':str(2*sum(empty_pair))}

def main():
 rows,hist=certificate();rng=random.Random(20261105);samples=1024
 for _ in range(samples):
  x=F(rng.randrange(1,65),128);h=x*F(rng.randrange(65),64);mu=F(rng.randrange(-32,33),64);v=[F(rng.randrange(-32,33),32) for _ in range(4)];u=[mu+h*z for z in v]
  actual=flat_R(x,u)+2*sum(u);identity=8*x*(1-x)+h*h*T(v)+4*mu*h*sum(v)+8*mu*mu
  assert actual==identity and actual<=8*x*(1-x)+8*(abs(mu)+h)**2
 for i in range(65):
  p=F(i,64);assert envelope_integral(p)==2+2*p*(1-p)<=F(5,2)
 out={'python':platform.python_version(),'arithmetic':'Fraction exact','seed':20261105,'bernstein_coefficients':360,'coefficient_histogram':hist,'certificate_rows':rows,'common_interval_regressions':samples,'exact_envelope_integral_cases':65,'example':exact_example(),'scope':'360 coefficients prove T<=8 on all 24 closed order simplices; finite regressions check identities; actual-set transfer is analytic.'}
 target=Path(__file__).resolve().parents[1]/'results/near_rectangle_python.json';target.write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k not in ['certificate_rows','example']});print(out['example'])
if __name__=='__main__':main()
