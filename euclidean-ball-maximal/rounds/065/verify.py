#!/usr/bin/env python3
"""Exact birth budgets, actual private-source routing, and clique margins."""
from pathlib import Path
from fractions import Fraction as F
from math import factorial,comb
import importlib.util,hashlib,json,random,sys,tempfile
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[2]
def module(path,name):
 s=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def length(a):return max(F(0),a[1]-a[0]) if a else F(0)
def overlap(a,b):return max(F(0),min(a[1],b[1])-max(a[0],b[0])) if a and b else F(0)
def born(a,b):return length(a)-overlap(a,b)
def halo(atoms,weights,alpha):
 active=[x for (x,w),q in zip(atoms,weights) if q>0];q=sum(weights)
 if not active:return None
 a=(max(active)-q/(2*alpha),min(active)+q/(2*alpha))
 return a if length(a)>0 else None
def lens_unequal(n,r,s,a):
 # Centers -a,+a, radii r,s; volume divided by omega_n.
 if a==0:return min(r,s)**n
 l=max(-a-r,a-s);h=min(-a+r,a+s)
 if l>=h:return F(0)
 m=(n-1)//2;c=F(factorial(2*m+2),4**(m+1)*factorial(m+1)*factorial(m))
 def primitive(z,t):return sum(((-1)**j*comb(m,j)*t**(2*(m-j))*z**(2*j+1)/F(2*j+1) for j in range(m+1)),F(0))
 pivot=(r*r-s*s)/(4*a);cut=min(h,max(l,pivot))
 return c*(primitive(cut-a,s)-primitive(l-a,s)+primitive(h+a,r)-primitive(cut+a,r))
def main():
 old=module('rounds/064/verify.py','old65');rec=module('rounds/046/threshold_probe.py','rec65');flow=module('rounds/047/flow_probe.py','flow65')
 geometry=[];chains=[];stars=[];cliques=[]
 for batch,ns in enumerate(((1,3,5),(9,17,33),(65,97,129))):
  for n in ns:
   count=0;maximum=F(0)
   for a in (F(0),F(1,1000),F(1,8),F(1,2),F(9,10)):
    full=old.lens(n,F(1),a);assert full==lens_unequal(n,F(1),F(1),a)
    for r in (F(1,4),F(1,2),F(3,4),F(99,100)):
     added=full-lens_unequal(n,r,F(1),a);budget=1-r**n
     assert 0<=added<=budget;count+=1;maximum=max(maximum,added/budget)
   geometry.append(dict(batch=batch,n=n,checks=count,max_birth_ratio=str(maximum)))
 for batch,N in enumerate((4,16,64)):
  for serial in range(6):
   seed=650000+100*batch+serial;rng=random.Random(seed)
   xs=sorted(rng.sample(range(-2048,2049),N));raw=[rng.randint(1,128) for _ in xs]
   atoms=[(F(x,1024),F(w,sum(raw))) for x,w in zip(xs,raw)];order=list(range(N));rng.shuffle(order)
   qs=[F(0)]*N;prev=None;variable_prev=None;qprev=F(0);vprev=F(0);births=F(0);variation=F(0);fixedalpha=F(2**batch,2)
   for i in order:
    qs[i]=atoms[i][1];q=sum(qs);cur=halo(atoms,qs,fixedalpha)
    birth=born(cur,prev);assert birth<=(q-qprev)/fixedalpha
    births+=birth;variation+=born(cur,prev)+born(prev,cur)
    assert births<=q/fixedalpha and variation<=2*q/fixedalpha
    varying=rng.choice((F(1,4),F(1,2),F(1),F(2),F(4),F(8)))
    vh=halo(atoms,qs,varying);v=q/varying
    assert born(vh,variable_prev)<=max(F(0),v-vprev)
    prev=cur;qprev=q;variable_prev=vh;vprev=v
   chains.append(dict(batch=batch,N=N,seed=seed,alpha=str(fixedalpha),birth_mass=str(fixedalpha*births),variation_mass=str(fixedalpha*variation)))
 with tempfile.TemporaryDirectory(prefix='euclidean65_') as tmp:
  runtime=rec.load(ROOT,Path(tmp)/'scratch.json')
  for batch,N in enumerate((3,7,15)):
   for serial in range(5):
    seed=651000+100*batch+serial;rng=random.Random(seed);xs=sorted(rng.sample([i for i in range(-512,513) if i],N));raw=[rng.randint(1,64) for _ in xs]
    atoms=[(F(0),F(1,4))]+[(F(x,2**18),F(w,8*sum(raw))) for x,w in zip(xs,raw)]+[(F(10),F(5,8))]
    alpha=rng.choice((F(1,2),F(1),F(2),F(4)));saved,radii,cells,_=flow.context(runtime,atoms,alpha)
    root=[F(0)]*len(atoms);root[0]=atoms[0][1];parts=[]
    for i in range(1,N+1):
     row=[F(0)]*len(atoms);row[i]=atoms[i][1];parts.append(row)
    source=[root]+parts;blocks=[root]+[[a+b for a,b in zip(root,row)] for row in parts]
    regions=[halo(atoms,q,alpha) for q in blocks];edges=sorted({t for a in regions if a for t in a})
    demand=[F(0)]*(N+1);taskchecks=0
    for lo,hi,J,trace in cells:
     nodes=[lo]+[t for t in edges if lo<t<hi]+[hi]
     for l,h in zip(nodes,nodes[1:]):
      x=(l+h)/2;eligible=[i for i,a in enumerate(regions) if a and a[0]<x<a[1]]
      if not eligible:continue
      chosen=min(eligible)
      for k,(blo,bhi) in rec.intervals(trace,alpha/4,alpha/2).items():
       assert all(abs(x-y)<radii[k] for (y,w),q in zip(atoms,source[chosen]) if q)
       demand[chosen]+=(h-l)*(bhi-blo)/(alpha/4)*trace[k-1];taskchecks+=1
    assert demand[0]<=F(3,8)*sum(root)
    for d,part in zip(demand[1:],parts):assert d<=F(3,4)*sum(part)
    assert sum(demand)<=F(3,8)*sum(root)+F(3,4)*sum(sum(row) for row in parts)
    stars.append(dict(batch=batch,N=N,seed=seed,alpha=str(alpha),task_checks=taskchecks,root_demand=str(demand[0]),private_demand=str(sum(demand[1:])),positive_private_branches=sum(d>0 for d in demand[1:]),maximum_private_load=str(max(d/sum(row) for d,row in zip(demand[1:],parts)))))
 for batch,ms in enumerate(((3,4),(5,6),(7,8))):
  for m in ms:
   N=2**m;n=N;eta=F(1,2**(20+3*m));shift=2*eta
   cl,ch=old.root_interval(F(3,2),n);al,ah=old.root_interval(F(4,3),n)
   rl,rh=old.root_interval(F(2),n);_,r1h=old.root_interval(F(4*N),n)
   flo,_=old.root_interval(3*cl*cl,2)
   margins=[cl-shift-ah,rl-ch-shift,flo-shift-r1h]
   assert min(margins)>0
   compressed=(rh+eta)**2-cl**2
   assert 0<compressed<F(1,n)
   w=F(1,N);vol_before=F(8,2**(m-1));vol_J=F(8,2**m)
   assert 2*w/vol_before==F(1,8) and 2*w/vol_J==F(1,4)
   assert 2*w/F(8,2**(m+1))==F(1,2) and 2*w/F(8,2**(m+2))==1
   cliques.append(dict(batch=batch,N=N,n=n,eta_over_R=str(eta),minimum_margin=str(min(margins)),compressed_radius_squared_upper=str(compressed),increment_overlap_lower=str(F(N-1,6)),pair_union_allocation_upper=str(F(3,8)*F(N-1,N**(N//2)))))
 # Purely one-sided estimate: tiny remote addition destroys almost all old halo.
 losses=[]
 for m in (8,16,32):
  eps=F(1,2**m);a=[(F(0),1-eps),(F(10),eps)]
  before=halo(a,[1-eps,F(0)],F(1));after=halo(a,[1-eps,eps],F(1));assert after is None
  assert born(before,after)==1-eps
  losses.append(dict(epsilon=str(eps),loss_over_added_mass=str((1-eps)/eps),total_variation=str(2*(1-eps))))
 deps=json.loads((ROOT/'rounds/064/verification.json').read_text())['sha256'];deps['rounds/065/verify.py']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
 for p,h in deps.items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
 print(json.dumps(dict(round=65,status='passed',python=sys.version.split()[0],geometry=geometry,chains=chains,actual_stars=stars,cliques=cliques,one_sided=losses,sha256=deps,main_theorem_proved=False,scope='Exact finite checks of birth charging and clique margins; no coverage of all actual E by a bounded-increment tree is asserted.'),indent=2))
if __name__=='__main__':main()
