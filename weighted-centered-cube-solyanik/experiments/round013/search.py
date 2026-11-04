import numpy as np, math, json, sys, time
from fractions import Fraction as F
from scipy.optimize import differential_evolution
from pathlib import Path
P=Path(__file__).resolve().parent

def merge(q):
 q=sorted((a,b) for a,b in q if b>a); out=[]
 for a,b in q:
  if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(b,out[-1][1]))
  else: out.append((a,b))
 return out

def halo(E,alpha,exact=False):
 zero=F(0) if exact else 0.; two=F(2) if exact else 2.
 ends=[v for ab in E for v in ab]; length=sum((b-a for a,b in E),zero)
 lo=ends[0]-length/(two*alpha); hi=ends[-1]+length/(two*alpha)
 def mass(a,b): return sum((max(zero,min(b,d)-max(a,c)) for c,d in E),zero)
 result=list(E)
 for e in ends:
  cuts=sorted(set([lo,hi]+[(e+v)/two for v in ends if lo<(e+v)/two<hi]))
  def g(x): return mass(min(e,two*x-e),max(e,two*x-e))-two*alpha*abs(x-e)
  for a,b in zip(cuts,cuts[1:]):
   x1=(two*a+b)/3; x2=(a+two*b)/3
   y1,y2=g(x1),g(x2); slope=(y2-y1)/(x2-x1); ga=y1+slope*(a-x1); gb=y1+slope*(b-x1)
   if not exact:
    eps=1e-12
    if abs(ga)<eps: ga=0.
    if abs(gb)<eps: gb=0.
   if ga>0 and gb>0: result.append((a,b))
   elif ga>0 and gb<=0: result.append((a,a+(b-a)*ga/(ga-gb)))
   elif gb>0 and ga<=0: result.append((a+(b-a)*(-ga)/(gb-ga),b))
 return merge(result)

def wm(E,tau):
 def f(x): return math.copysign(abs(float(x))**tau,float(x)) if x else 0.
 return sum(f(b)-f(a) for a,b in E)
def score(E,t,tau,exact=False): return (wm(halo(E,1-F(str(t)) if exact else 1-t,exact),tau)/wm(E,tau)-1)/t**tau

def mk(v,t,n):
 L,a,b=v[:3]; R=2-L; holes=[(-t*a,t*b)]; idx=3
 if n==4:
  c,p=v[idx:idx+2]; idx+=2; w=t*c
  if w>=L-t*a: return None
  available=L-t*a-w; distance=t*math.expm1(p*math.log1p(available/t)); start=-t*a-w-distance; holes.append((start,start+w))
 c,p=v[idx:idx+2]; w=t*c
 if w>=R-t*b: return None
 available=R-t*b-w; distance=t*math.expm1(p*math.log1p(available/t)); start=t*b+distance; holes.append((start,start+w))
 holes.sort(); endpoints=[-L]+[x for h in holes for x in h]+[R]
 if any(endpoints[i]>=endpoints[i+1] for i in range(len(endpoints)-1)): return None
 return list(zip(endpoints[::2],endpoints[1::2]))

def control():
 for t in [.0625,.01,1/256]:
  for A,h1,h2,B in [(2.,t,t,2.),(1.4,t*.7,t*1.3,2.1),(1.2,t*2,t,1.6)]:
   E=[(-A,-h1),(h2,B)]; H=halo(E,1-t)
   l=max(-h1,-A+(h1+h2)/(2*t)); r=min(h2,B-(h1+h2)/(2*t))
   expected=merge(E+([(l,r)] if l<r else []))
   assert len(H)==len(expected) and all(abs(float(a-c))<1e-8 and abs(float(b-d))<1e-8 for (a,b),(c,d) in zip(H,expected)),(H,expected)
   EF=[tuple(F(str(x)) for x in ab) for ab in E]; HF=halo(EF,1-F(str(t)),True)
   assert len(HF)==len(expected) and all(abs(float(a)-c)<1e-8 and abs(float(b)-d)<1e-8 for (a,b),(c,d) in zip(HF,expected))
 print('control: 9 float and fraction comparisons passed',flush=True)

def main():
 control(); records=[]; seed=23001
 for tau in [.25,.5,.75]:
  for k in [4,8,12,16]:
   t=2**-k
   baseline=1/((1+t)**tau-t**tau)
   for n in [3,4]:
    bounds=[(.15,1.85),(.01,3.5),(.01,3.5)]+[(.001,3.5),(.00001,.99999)]*(n-2)
    count=[0]; start=time.time()
    def obj(v):
     count[0]+=1; E=mk(v,t,n)
     if E is None:return 1e3
     return -score(E,t,tau)
    result=differential_evolution(obj,bounds,seed=seed,popsize=12,maxiter=180,tol=0,polish=False)
    E=mk(result.x,t,n); rational=[(F(float(a)).limit_denominator(10**9),F(float(b)).limit_denominator(10**9)) for a,b in E]
    H=halo(rational,1-F(1,2**k),True); exactscore=(wm(H,tau)/wm(rational,tau)-1)/t**tau
    rec={'tau':tau,'k':k,'segments':n,'seed':seed,'evaluations':count[0],'float_score':-result.fun,'rational_score':exactscore,'baseline_score':baseline,'seconds':round(time.time()-start,2),'parameters':result.x.tolist(),'intervals':[[str(a),str(b)] for a,b in rational],'halo':[[str(a),str(b)] for a,b in H]}
    records.append(rec); seed+=1
    (P/'candidates.json').write_text(json.dumps(records,indent=2))
    print(tau,k,n,count[0],round(-result.fun,7),round(exactscore,7),'baseline',round(baseline,7),'s',rec['seconds'],flush=True)
if __name__=='__main__':main()
