"""Exact rational polynomial certificate, no floating arithmetic."""
from fractions import Fraction as F
from math import comb
import json
from pathlib import Path
import argparse

def add(a,b):
 c=[F(0)]*max(len(a),len(b))
 for i,v in enumerate(a):c[i]+=v
 for i,v in enumerate(b):c[i]+=v
 return trim(c)
def trim(a):
 while len(a)>1 and not a[-1]:a.pop()
 return a
def mul(a,b):
 c=[F(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return trim(c)
def diff(a):return [F(i)*a[i] for i in range(1,len(a))] or [F(0)]
def scale(a,s):return [x*s for x in a]
def shift(a,origin,width=F(1)):
 c=[F(0)]*len(a)
 for j,x in enumerate(a):
  for i in range(j+1):c[i]+=x*comb(j,i)*origin**(j-i)*width**i
 return trim(c)
def poly(b,j):
 ans=[F(0)]
 for k in range(j+1):ans=add(ans,[F(comb(b,k)*comb(b-1,i))*(-k)**(b-1-i) for i in range(b)])
 return ans
def defect(p):return add(mul(diff(p),diff(p)),scale(mul(p,diff(diff(p))),-1))
def bernstein(a):
 d=len(a)-1
 return [sum((a[i]*F(comb(k,i),comb(d,i)) for i in range(k+1)),F(0)) for k in range(d+1)]
def pieces(q,a,b,depth=0):
 bs=bernstein(shift(q,a,b-a))
 if min(bs)>=0:return [{'interval':[str(a),str(b)],'bernstein_coefficients':[str(x) for x in bs]}]
 assert depth<20,(a,b,depth)
 mid=(a+b)/2
 return pieces(q,a,mid,depth+1)+pieces(q,mid,b,depth+1)

rows=[]
for b in range(4,8):
 cells=[]
 for j in range(b):
  p=poly(b,j);q=defect(p)
  cells.append({'density_piece':j,'power_defect':[str(x) for x in q],'certificates':pieces(q,F(j),F(j+1))})
 q=defect(poly(b,b));tail=b
 while min(shift(q,F(tail)))<0:
  cells.append({'density_piece':b,'power_defect':[str(x) for x in q],'certificates':pieces(q,F(tail),F(tail+1))});tail+=1;assert tail<100
 rows.append({'b':b,'cells':cells,'tail_from':tail,'tail_nonnegative_power_coefficients':[str(x) for x in shift(q,F(tail))]})
 print('b',b,'finitepieces',sum(len(x['certificates']) for x in cells),'tail',tail)
result={'arithmetic':'Python Fraction exact rational; Q=Rprime^2-R Rsecond','base_orders':[4,5,6,7],'rows':rows,'b3_right_defect_at_1':str(sum(defect(poly(3,1))))}
assert result['b3_right_defect_at_1']=='-4'
counter_x=F(101,100);counter_q=sum((a*counter_x**i for i,a in enumerate(defect(poly(3,1)))),F(0))
assert counter_q==F(-2398,625)
result['b3_interior_counterexample']={'x':str(counter_x),'Q':str(counter_q)}
parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
target=Path(__file__).resolve().parents[1]/'rounds'/'r027'/'logconcavity_certificate.json'
if args.verify:
 assert json.loads(target.read_text())==result, 'Stored certificate differs from exact reconstruction'
 print('Stored exact certificate verified.')
else:
 target.write_text(json.dumps(result,indent=2)+'\n')
