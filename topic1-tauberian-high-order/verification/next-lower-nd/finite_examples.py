from fractions import Fraction as F
from math import comb
from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
from pathlib import Path
import json,time

def dec(x,r):
 with localcontext() as c:c.prec=30;c.rounding=r;return str(Decimal(x.numerator)/Decimal(x.denominator))
def calc(n):
 u=F(1,2**16);g=1/(1+2*u);kap=F(8*n,3)*g*g;C=10**6;pairs=comb(n,2);lam0=pairs*kap**3/(8*n*n*g);lam=lam0+C*u*u;w=u*g+lam*u**6;rho=1-2*w;a=(2*u)**n;e=kap*u**3;d=(1-e)**(n-1)*(1+(n-1)*e);M=a*rho**n;m0=a*g**n;D=(2*w)**n-M;K=a/(2*w)**(n-2);mult=3**n;width=w-e
 FF=lambda h:K*sum((F(comb(n,j))*rho**(n-j)*h**j for j in range(1,n+1)),F(0))
 guards={'w_lt_quarter':w<F(1,4),'e_lt_w_over_2':e<w/2,'strain':(n-1)*e<F(1,12),'det_gt_outer':d>(2*w)**n,'double_high_budget':FF(width)>2*mult*e**2}
 if not all(guards.values()):return {'n':n,'guards':guards,'status':'guard_failed'}
 def bracket(m):
  l=F(0);h=width
  for _ in range(400):
   mid=(l+h)/2
   if FF(mid)<m*e**2:l=mid
   else:h=mid
  i=l*2**400/width;assert i.denominator==1
  return l,h,hex(i.numerator)
 l,h,i=bracket(1);lh,hh,j=bracket(mult)
 prim=lambda z:sum((F(comb(n,i))*rho**(n-i)*z**(i+1)/F(i+1) for i in range(1,n+1)),F(0))
 Al=2*(e*l-K/e*prim(h));Au=2*(e*h-K/e*prim(l));bl=pairs*2*Al*(2*w)**(n-2);bu=pairs*2*Au*(2*w)**(n-2)
 if n>=3:
  bl+=8*n*comb(n-1,2)*lh*(e+lh)**2*(2*w)**(n-3);bu+=8*n*comb(n-1,2)*hh*(e+hh)**2*(2*w)**(n-3)
 G=(M-d*m0)/((d-M)*(1-m0)*a*(2*u)**6)
 guards.update(hstar_enclosed=FF(l)<e**2<=FF(h),hH_enclosed=FF(lh)<mult*e**2<=FF(hh),hH_inside=hh<width,area_positive=Al>0,budget_strict=D-bu>0,G_strict=G>0)
 return {'n':n,'u':str(u),'C':C,'kappa':str(kap),'lambda':str(lam),'roots':{'hstar':i,'hH':j,'scale':'w-e','pow2':400},'budget_bounds':[dec(D-bu,ROUND_FLOOR),dec(D-bl,ROUND_CEILING)],'G':[dec(G,ROUND_FLOOR),dec(G,ROUND_CEILING)],'guards':guards,'status':'pass' if all(guards.values()) else 'arithmetic_fail','scope':'uses supplied sharpened three-face bound; finite periodic formula only'}
if __name__=='__main__':
 t=time.perf_counter();rows=[calc(n) for n in (2,3,4)];obj={'cases':rows,'all_pass':all(x['status']=='pass' for x in rows),'seconds':time.perf_counter()-t};Path('finite_examples.json').write_text(json.dumps(obj,indent=2)+'\n');print(json.dumps({'all_pass':obj['all_pass'],'cases':[(r['n'],r['status']) for r in rows],'seconds':obj['seconds']}))
