#!/usr/bin/env python3
"""Exact finite audits for farthest-source translations and full-block budgets."""
from pathlib import Path
from fractions import Fraction as F
from math import factorial,comb
import json,random,sys,hashlib,importlib.util,tempfile
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[2]
def module(path,name):
 s=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def norm(a):return dot(a,a)
def lens(n,t,a):
 # |B(-ae1,t) intersect B(ae1,t)| / omega_n, for odd n.
 if t<=a:return F(0)
 m=(n-1)//2;c=F(factorial(2*m+2),4**(m+1)*factorial(m+1)*factorial(m))
 return 2*c*sum(((-1)**j*comb(m,j)*t**(2*(m-j))*(t**(2*j+1)-a**(2*j+1))/F(2*j+1) for j in range(m+1)),F(0))
def root_interval(value,n,bits=100):
 # Dyadic outward enclosure, integer powers only.
 scale=2**bits;lo=0;hi=2*scale
 while F(hi,scale)**n<value:hi*=2
 assert F(hi,scale)**n>=value
 while hi-lo>1:
  mid=(lo+hi)//2
  if F(mid,scale)**n<=value:lo=mid
  else:hi=mid
 return F(lo,scale),F(hi,scale)
def cover_obstruction():
 records=[]
 for batch,Ns in enumerate(((2,3,4),(8,16,32),(64,128,256))):
  for N in Ns:
   n=max(N+1,9);u=F(1,2*N)
   dl,dh=root_interval(u,n);cl,ch=root_interval(1+u,n)
   tl,th=(cl*cl-1)/4,(ch*ch-1)/4
   corelo,_=root_interval(1+tl,2);_,corehi=root_interval(1+th,2)
   privlo,_=root_interval(dl*dl+tl,2);_,privhi=root_interval(dh*dh+th,2)
   foreignlo,_=root_interval(1+(1+dl)**2+tl,2)
   # roots for farthest distances can exceed 2; foreign radius is >2 here.
   r1lo,r1hi=root_interval(F(16),n);r4lo,r4hi=root_interval(F(2),n)
   rmaxlo,rmaxhi=root_interval(F(4),n)
   eta=F(1,2**(20+3*(N-1).bit_length()));shift=2*eta
   margins=[corelo-shift-1,privlo-shift-dh,cl-corehi-shift,
            r4lo-max(corehi,privhi)-shift,foreignlo-shift-r1hi,
            foreignlo-shift-rmaxhi]
   assert all(v>0 for v in margins)
   q=F(1,4)*(1+u)
   assert q/4<F(1,8)<q/2
   assert q>F(1,4) and q/F(1,2)>F(1,2)
   assert (1+u)/(corelo-shift)**n<2
   # Far source at 100*e_(N+1): distance at observation neighborhoods > 90.
   assert th<1 and r1hi<2
   records.append(dict(batch=batch,N=N,n=n,core_mass='1/4',private_mass=str(F(1,8*N)),
                       source_and_observation_radius=str(eta),minimum_distance_margin=str(min(margins)),
                       necessary_overlap=str(F(N)-F(1,2))))
 return records
def main():
 maps=[];shells=[];blocks=[]
 for batch,dimensions in enumerate(((1,2,4),(8,16,32),(64,128,256))):
  for n in dimensions:
   seed=640000+1000*batch+n;rng=random.Random(seed);N=(4,8,16)[batch]
   # All coordinates integral: strong monotonicity and expansion residuals exact.
   centers=[tuple(rng.randint(-20,20) for _ in range(n)) for _ in range(N)]
   points=[tuple(rng.randint(-40,40) for _ in range(n)) for _ in range(24)]
   mapped=[]
   for x in points:
    k=max(range(N),key=lambda i:(norm(sub(x,centers[i])),-i));mapped.append((x,centers[k],sub(x,centers[k])))
   residuals=[]
   for j,(x,a,tx) in enumerate(mapped):
    for y,b,ty in mapped[:j]:
     dx=sub(x,y);dt=sub(tx,ty)
     assert dot(sub(a,b),dx)<=0
     residual=dot(dt,dx)-norm(dx);assert residual>=0
     assert norm(dt)>=norm(dx)
     residuals.append(residual)
   # Candidate same translated vector violates at least one farthest-cell condition.
   disjoint=0
   for i,a in enumerate(centers):
    for b in centers[:i]:
     if a==b:continue
     d=sub(a,b)
     for z in points[:4]:
      plus=norm(sub(z,tuple(-v for v in d)))-norm(z)
      minus=norm(sub(z,d))-norm(z)
      assert plus+minus==2*norm(d)>0;assert max(plus,minus)>0;disjoint+=1
   maps.append(dict(batch=batch,n=n,N=N,seed=seed,pairs=len(residuals),translation_checks=disjoint,min_monotonicity_residual=min(residuals)))
 for batch,dimensions in enumerate(((1,3,5),(9,17,33),(65,97,129))):
  for n in dimensions:
   count=0;best=F(0)
   for a in (F(0),F(1,1000),F(1,4),F(3,4),F(99,100)):
    for r in (F(0),F(1,4),F(1,2),F(99,100),F(999,1000)):
     v=lens(n,F(1),a)-lens(n,r,a);budget=1-r**n
     assert 0<=v<=budget
     # Full source sphere of radius a has farthest distance |x|+a.
     concentric=max(F(0),1-a)**n-max(F(0),r-a)**n
     assert 0<=concentric<=budget
     best=max(best,v/budget);count+=2
   assert lens(n,F(1),F(0))==1
   shells.append(dict(batch=batch,n=n,checks=count,max_shell_ratio=str(best)))
 prev=module('rounds/058/verify.py','prev64');flow=module('rounds/047/flow_probe.py','flow64');rec=module('rounds/046/threshold_probe.py','rec64')
 with tempfile.TemporaryDirectory(prefix='euclidean64_') as tmp:
  runtime=rec.load(ROOT,Path(tmp)/'scratch.json')
  for batch,N in enumerate((4,7,12)):
   for serial in range(5):
    seed=641000+100*batch+serial;rng=random.Random(seed)
    xs=sorted(rng.sample(range(-1024,1025),N));raw=[rng.randint(1,64) for _ in xs]
    atoms=[(F(x,2048),F(3,8)*F(w,sum(raw))) for x,w in zip(xs,raw)]+[(F(10),F(5,8))]
    alpha=rng.choice((F(1,2),F(1),F(2),F(4)))
    saved,radii,cells,_=flow.context(runtime,atoms,alpha)
    family=[(tuple(range(i,j+1)),F(1)) for i in range(N) for j in range(i,N)]
    family += [(tuple(sorted(rng.sample(range(N),rng.randint(1,N)))),F(rng.randint(1,3),4)) for _ in range(10)]
    maxload=F(0);active=0;checks=0
    for S,factor in family:
     q=sum((atoms[i][1] for i in S),F(0))*factor
     left=max(atoms[i][0] for i in S)-q/(2*alpha)
     right=min(atoms[i][0] for i in S)+q/(2*alpha)
     if left>=right:continue
     demand=F(0);volume=F(0)
     for lo,hi,J,trace in cells:
      l=max(lo,left);h=min(hi,right)
      if l>=h:continue
      x=(l+h)/2;volume+=h-l
      assert prev.maximal(atoms,x)<=2*alpha
      for k,(blo,bhi) in rec.intervals(trace,alpha/4,alpha/2).items():
       assert all(abs(x-atoms[i][0])<radii[k] for i in S)
       demand+=(h-l)*(bhi-blo)/(alpha/4)*trace[k-1];checks+=1
     assert volume<=q/(2*alpha);assert demand<=F(3,8)*q
     if demand:active+=1
     maxload=max(maxload,demand/q)
    blocks.append(dict(batch=batch,N=N,seed=seed,alpha=str(alpha),blocks=len(family),positive_blocks=active,task_checks=checks,max_load=str(maxload)))
 deps=json.loads((ROOT/'rounds/063/verification.json').read_text())['sha256']
 deps['rounds/064/verify.py']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
 for p,h in deps.items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
 print(json.dumps(dict(round=64,status='passed',python=sys.version.split()[0],translation=maps,shells=shells,actual_blocks=blocks,cover_obstruction=cover_obstruction(),sha256=deps,main_theorem_proved=False,scope='Exact finite checks; arbitrary dimension and general source statements require the report proofs. No global block cover is asserted.'),indent=2))
if __name__=='__main__':main()
