#!/usr/bin/env python3
"""Independent exact verification of the rational tau=1/10 two-atom witness.

No DE/grid search is performed.  Strict halo blocks are enumerated from the
two atoms, root bounds use Fraction powers, and the single-atom bound is the
recorded concavity/tangent certificate.
"""
from fractions import Fraction as Q
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parent
TAU=Q(1,10);PTS=(Q(-1),Q(1));MASSES=(Q(1,4),Q(15,8))
def blocks(points,masses):
 out=[]
 for i in range(len(points)):
  total=Q(0)
  for j in range(i,len(points)):
   total+=masses[j];radius=total/2;left,right=points[i],points[j]
   if radius>(right-left)/2:out.append((right-radius,left+radius,i,j,total))
 return out
def main():
 B=blocks(PTS,MASSES)
 intervals=[(x[0],x[1]) for x in B]
 expected=[(Q(-9,8),Q(-7,8)),(Q(-1,16),Q(1,16)),(Q(1,16),Q(31,16))]
 assert intervals==expected
 # Preserve strict adjacency: central and positive components share an excluded endpoint.
 assert expected[1][1]==expected[2][0] and all(a<b for a,b in expected)
 # Both atoms are at |x|=1, so the power-density factor is exactly one.
 cost=sum(MASSES,Q(0));assert cost==Q(17,8)
 # Formula for tau*halo/cost after signed primitive integration.
 a,b,c,d=Q(5059,5000),Q(2467,2500),Q(3789,5000),Q(10683,10000)
 root_checks={'positive_negative':a**10<Q(9,8),'negative_upper':b**10>Q(7,8),'central_lower':c**10<Q(1,16),'positive_lower':d**10<Q(31,16)}
 assert all(root_checks.values());lower=a-b+c+d;two=Q(8,17)*lower
 assert lower==Q(18511,10000) and two==Q(18511,21250)>Q(87,100)
 # Concavity/tangent all-single certificate (all arithmetic exact).
 A=Q(535887,500000);Bcoef=Q(53589,1000000);D=Q(173,100)-A;C=Q(173,100)-Bcoef;s0=Q(4361,100000)
 tangent_checks={'A_upper_2root':A**10>2,'B_above_derivative':(Q(1)/(10*Bcoef))**10<2**9,
                 'C_above_tangent_slope':s0**9>(Q(1)/(10*C))**10,'D_above_tangent_intercept':s0<(10*D/Q(9))**10}
 assert all(tangent_checks.values())
 single_upper=Q(173,200);assert A/2<single_upper and single_upper<two
 out={'input':{'tau':'1/10','points':[str(x) for x in PTS],'masses':[str(x) for x in MASSES],'cost':str(cost)},
      'strict_blocks':[{'left':str(a),'right':str(b),'atoms':[i,j],'mass':str(m)} for a,b,i,j,m in B],
      'root_checks':root_checks,'tau_two_lower':str(two),'single_tangent_checks':tangent_checks,
      'single_upper':str(single_upper),'gap':str(two-single_upper),
      'scope':'Endpoint/Fatou smoothing is analytic, not certified by this atomic computation.'}
 (ROOT/'independent_certificate.json').write_text(json.dumps(out,indent=2));print('independent_fraction_certificate PASS: tau R_two > 18511/21250 > 173/200 >= single certificate')
if __name__=='__main__':main()
