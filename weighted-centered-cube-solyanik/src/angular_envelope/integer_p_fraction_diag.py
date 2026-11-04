#!/usr/bin/env python3
"""Small exact Fraction diagnostics for integer-p general-envelope identities."""
from fractions import Fraction as Q
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parent
def values(p,d,q):
 delta=d/2;a=1-delta/q;b=1+delta/(1-q)
 G=q*a**p+(1-q)*b**p-1
 g2=p*(p-1)*delta**2*(a**(p-2)/q**3+b**(p-2)/(1-q)**3)
 dd=p*Q(1,2)*(b**(p-1)-a**(p-1))
 return a,b,G,g2,dd
def main():
 rows=[];bad=[]
 for p in (2,3,4,8,16):
  switch=Q(2,p)
  for d in (switch/Q(4),switch*Q(3,4),min(Q(3,2),switch+Q(1,8))):
   # choose q strictly above delta, keeping q<1.
   q=(d/Q(2)+Q(1))/2
   a,b,G,g2,dd=values(p,d,q);ok=(0<a<1<b and G>0 and g2>0 and dd>0)
   rows.append({'p':p,'d':str(d),'q':str(q),'a':str(a),'G':str(G),'G2':str(g2),'partial_d_G':str(dd),'pass':ok})
   if not ok:bad.append(rows[-1])
  d=switch;H=(1-d/2)**(1-p)-1
  rows.append({'p':p,'d_switch':str(d),'outer_H_exact':str(H),'pass':H>0})
 (ROOT/'integer_p_fraction.json').write_text(json.dumps({'rows':rows,'failed':bad},indent=2));assert not bad
 print('integer_p_fraction_rows=%d failed=0' % len(rows))
if __name__=='__main__':main()
