"""Finite Arb box diagnostic; no implicit extension across denominator zeros."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict,Counter
import runpy,sympy as s,json,heapq,hashlib,time
from flint import arb,ctx
from sympy.printing.pycode import PythonCodePrinter
class ArbPrinter(PythonCodePrinter):
 def _print_Rational(self,e): return '(arb(%d)/%d)'%(e.p,e.q)
def fn(e):return s.lambdify((y,z),e,modules=[{'arb':arb}],printer=ArbPrinter(),cse=True)
ctx.prec=160
OUT=Path(__file__).resolve().parent
SOURCE=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/cb11_frame_green_locked.py')
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()=='ffbacd4f1a466f8f20d11dc59c458be2980b45be964751b8ea1a980cbbf1d849'
ns=runpy.run_path(str(SOURCE));y,z=s.symbols('y z'); groups=defaultdict(list)
for cone,ks,w in ns['_prepared'](F(1,3)):
 bs=[sum(r[i]*(1,y,z)[i] for i in range(3))/s.Rational(k.numerator,k.denominator) for r,k in zip(cone,ks)]
 for j,b in enumerate(bs):
  c=s.Rational(w.numerator,w.denominator)/s.prod(b-bs[k] for k in range(3) if k!=j)
  if str(b).startswith('-'):b=-b
  groups[b].append(c)
rows=[];factor_set=set()
for b,terms in groups.items():
 c=s.factor(sum(terms));r=s.factor(c/b**2);num,den=s.fraction(r)
 factors=s.factor_list(den)[1]
 for p,n in factors: factor_set.add(str(p))
 rows.append((b,r,fn(b),fn(r),[fn(p) for p,n in factors]))
(OUT/'grouped_formula.json').write_text(json.dumps({'formula':'8*pi^3*G=sum r(y,z)*b(y,z)^2*log(abs(b(y,z)))','rows':[{'b':str(b),'r':str(r)} for b,r,*_ in rows],'denominator_factors':sorted(factor_set)},indent=2)+'\n')
def aq(x):return arb(x.numerator)/x.denominator
# b^2 log|b| is continuous at zero and bounded via scalar extrema.
def flog(b):
 if not b.contains(0):return b*b*abs(b).log()
 a=max(abs(b.lower()),abs(b.upper()))
 if a==0:return arb(0)
 q=a*a*a.log();v=-arb(1)/(2*arb(1).exp())
 lo=v  # global minimum of s^2 log|s|; safe also for narrow zero-crossing balls
 hi=max(arb(0),q)
 return lo.union(hi)
def evaluate(box):
 a,b,c,d=box;Y=aq((a+b)/2)+arb(0,aq((b-a)/2));Z=aq((c+d)/2)+arb(0,aq((d-c)/2))
 for _,_,_,_,fs in rows:
  if any(f(Y,Z).contains(0) for f in fs):return None,'coalescence_denominator'
 total=arb(0)
 for _,_,bf,rf,_ in rows:total+=rf(Y,Z)*flog(bf(Y,Z))
 v=total/(8*arb.pi()**3)
 if not v.is_finite():return None,'nonfinite'
 return v,'certified' if v.lower()>=arb(13)/2500 else 'interval_width_or_low_value'
def main():
 start=time.time(); todo=[];serial=0;heapq.heappush(todo,(-F(1),serial,(F(0),F(1),F(0),F(1))));tested=0;accepted=[];counts=Counter();best=None
 while todo and tested<4095:
  _,_,box=heapq.heappop(todo);a,b,c,d=box
  if c>b:continue
  v,status=evaluate(box);tested+=1;counts[status]+=1
  if status=='certified':accepted.append((box,str(v)));continue
  if b-a>=d-c:m=(a+b)/2;children=[(a,m,c,d),(m,b,c,d)]
  else:m=(c+d)/2;children=[(a,b,c,m),(a,b,m,d)]
  for ch in children:
   if ch[2]>ch[1]:continue
   serial+=1;heapq.heappush(todo,(-(ch[1]-ch[0])*(ch[3]-ch[2]),serial,ch))
 unresolved=[]
 for _,_,box in todo:
  # Final leaves were generated but not evaluated; do not claim leaf statuses.
  unresolved.append([str(t) for t in box])
 data={'status':'INCOMPLETE: not a global certificate','evaluated_regions':tested,'counts':dict(counts),'certified_regions':[{'box':[str(t) for t in b],'enclosure':v} for b,v in accepted],'unresolved_leaf_boxes':unresolved,'unresolved_leaf_count':len(todo),'seconds':time.time()-start,'precision_bits':ctx.prec,'domain':'0<=z<=y<=1; rectangles intersected with domain','initial_representation':'25 aggregated exact rational coefficient logarithms; zero projections continued by b^2 log|b|'}
 (OUT/'diagnostic.json').write_text(json.dumps(data,indent=2)+'\n');print(json.dumps({k:v for k,v in data.items() if k not in ('certified_regions','unresolved_leaf_boxes')},indent=2),flush=True)
if __name__=='__main__':main()
