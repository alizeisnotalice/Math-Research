from pathlib import Path
from fractions import Fraction as F
from math import comb
import json,hashlib
H=Path(__file__).resolve().parent;rows=[];r=F(3,5);eps=F(1,4)
for n in [1,2,4,8,16,64]:
 atoms=[(F(comb(n,k),2**n),(1+r)**k*(1-r)**(n-k),(1-r)**k*(1+r)**(n-k))for k in range(n+1)]
 assert sum(p*a for p,a,b in atoms)==sum(p*b for p,a,b in atoms)==1
 E=F(n)*r*r/(2*(1+r));m=E.numerator//E.denominator;delta=2*F(1,2**m)
 overlap=sum(p*a for p,a,b in atoms if b/a>eps);assert overlap<=delta
 for lam in [F(1,10),F(1),F(2),F(10),F(100)]:
  W=lam*sum(p for p,a,b in atoms if (a+b)/2>lam)
  def wb(t):return t*(sum(p for p,a,b in atoms if a>t)+sum(p for p,a,b in atoms if b>t))/2
  lo=wb(2*lam)-delta;hi=(1+eps)*wb(2*lam/(1+eps))+delta
  assert lo<=W<=hi
  rows.append(dict(n=n,lambda_exact=str(lam),weak_exact=str(W),lower_exact=str(lo),upper_exact=str(hi),anchor_overlap_exact=str(overlap),anchor_overlap_bound_exact=str(delta)))
out=dict(cases=len(rows),rows=rows,script_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Exact abstract finite probability verification of profile lemma, not an HL example.')
(H/'FINITE.json').write_text(json.dumps(out,indent=2)+'\n');print({'cases':len(rows),'passed':True})
