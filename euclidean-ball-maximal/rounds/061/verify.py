#!/usr/bin/env python3
"""Exact tests of load monotonicity and logarithmic height-boundary loss."""
from pathlib import Path
from fractions import Fraction as F
import importlib.util,sys,tempfile,random,json,hashlib
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[2]
def module(path,name):
 s=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def maximal(pieces,x):
 radii=sorted({abs(t-x) for a,b,h in pieces for t in (a,b)}-{F(0)})
 def density(r):return sum((h*max(F(0),min(b,x+r)-max(a,x-r)) for a,b,h in pieces),F(0))/(2*r)
 return max([density(r) for r in radii]+[sum((h for a,b,h in pieces if a<x<b),F(0))])

def main():
 p=module('rounds/060/verify.py','prof61');a=module('rounds/056/allocation_probe.py','alloc61')
 prev=module('rounds/058/verify.py','prev61');flow=module('rounds/047/flow_probe.py','flow61');rec=module('rounds/046/threshold_probe.py','rec61')
 records=[]
 with tempfile.TemporaryDirectory(prefix='euclidean61_') as tmp:
  old=rec.load(ROOT,Path(tmp)/'scratch.json')
  witness=json.loads((ROOT/'rounds/060/verification.json').read_text())['actual_frozen_counterexample']
  inputs=[(0,'round60_witness',[(F(z['x']),F(z['w'])) for z in witness['atoms']],F(1))]
  for batch,N in enumerate((3,5,8)):
   for serial in range(10):
    seed=610100+100*batch+serial;rng=random.Random(seed)
    xs=sorted(rng.sample(range(-512,513),N));ws=[rng.randint(1,64) for _ in xs]
    atoms=[(F(x,1024),F(3,8)*F(w,sum(ws))) for x,w in zip(xs,ws)]+[(F(10),F(5,8))]
    inputs.append((batch,f'seed{seed}',atoms,rng.choice((F(1,2),F(1),F(2),F(4)))))
  for batch,label,atoms,alpha in inputs:
   saved,r,groups,tasks=prev.capture_data(flow,rec,old,atoms,alpha);weights=[w for x,w in atoms]
   keys=sorted(groups,key=lambda S:(len(S),tuple(sorted(S))))
   gates=[{S:d*F(1+(i+shift)%4,4) for i,(S,d) in enumerate((S,groups[S]) for S in keys)} for shift in (0,2)]
   pa,_,_=p.profile(gates[0],weights,a);pb,_,_=p.profile(gates[1],weights,a)
   discrepancy=sum((w*abs(x-y) for w,x,y in zip(weights,pa,pb)),F(0))
   rhs=sum((abs(gates[0][S]-gates[1][S]) for S in keys),F(0));assert discrepancy<=rhs
   last=[F(0)]*len(weights);prior=F(0);variation=F(0);steps=0
   for stage in range(1,4):
    selected={S:groups[S] for i,S in enumerate(keys) if i%3<stage}
    cur,_,_=p.profile(selected,weights,a);mass=sum(selected.values(),F(0))
    assert all(x>=y for x,y in zip(cur,last))
    inc=sum((w*abs(x-y) for w,x,y in zip(weights,cur,last)),F(0));assert inc==mass-prior
    variation+=inc;last=cur;prior=mass;steps+=1
   assert variation==sum(groups.values(),F(0))
   records.append(dict(batch=batch,label=label,N=len(weights),D=saved['D'],groups=len(keys),monotone_steps=steps,contraction_left=str(discrepancy),contraction_right=str(rhs),total_variation=str(variation)))
 chain=[]
 for batch,Ls in enumerate(((1,2,3),(4,5,6),(7,8,9))):
  for L in Ls:
   N=2**L
   eps=F(1);c=F(N-1,N);cp=F(N-1+eps,N)
   old=[F(N-i,N) for i in range(1,N)];new=[(1-eps)*v for v in old]
   for vals,target,addition in ((old,c,F(0)),(new,cp,eps)):
    loads=[F(0)]*N;loads[0]+=addition
    for i,b in enumerate(vals):assert 0<=b<=1;loads[i]+=b;loads[i+1]+=1-b
    assert all(v==target for v in loads)
   moved=sum((x-y for x,y in zip(old,new)),F(0))/N
   assert moved==eps*F(N-1,2*N)
   # Independent canonical solver at representative sizes in each batch.
   solved=N in (4,16,128)
   if solved:
    groups={frozenset({i,i+1}):F(1,N) for i in range(N-1)};weights=[F(1,N)]*N
    ph,_,_=p.profile(groups,weights,a);assert ph==[c]*N
    groups[frozenset({0})]=eps/N;ph,_,_=p.profile(groups,weights,a);assert ph==[cp]*N
   # Embed the chain into actual first-exit subtask boxes of a common P.
   spacing=F(4,7*N);w=2*spacing/5;k=L+3;D=k+2
   atoms=[(i*spacing,w) for i in range(N)]+[(F(10),F(27,35))]
   assert sum((v for x,v in atoms),F(0))==1
   demand=6*spacing/4375;boxchecks=0
   for edge in sorted({0,(N-2)//2,N-2}):
    for d in (F(49,320),F(13,80),F(55,320)):
     for endpoint in (False,True):
      if endpoint and edge!=0:continue
      x=-d*spacing if endpoint else (edge+d)*spacing
      dist={}
      for y,v in atoms:dist[abs(y-x)]=dist.get(abs(y-x),F(0))+v
      mass=F(0);mp=F(0)
      for r,v in sorted(dist.items()):mass+=v;mp=max(mp,mass/(2*r))
      assert mp==F(1,5)/d and 1<mp<=2
      trace=[sum((v for y,v in atoms if abs(y-x)<F(4,2**j)),F(0))/(2*F(4,2**j)) for j in range(1,D+1)]
      assert trace[0]<=F(1,8)
      assert next(j for j,g in enumerate(trace,1) if g>F(1,2))==D
      for beta in (F(41,100),F(17,40),F(11,25)):
       selected=next(j for j,g in enumerate(trace,1) if g>beta)
       assert selected==k+int(endpoint) and trace[selected-1]==F(16,35)
       S={i for i,(y,v) in enumerate(atoms) if abs(y-x)<F(4,2**selected)}
       assert S==({0} if endpoint else {edge,edge+1})
       boxchecks+=1
   eta=spacing/1000
   for edge in sorted({0,(N-2)//2,N-2}):
    for endpoint in (False,True):
     if endpoint and edge!=0:continue
     ends=[-d*spacing if endpoint else (edge+d)*spacing for d in (F(3,20),F(7,40))]
     xa,xb=sorted(ends)
     for j in range(1,D+1):
      r=F(4,2**j)
      for y,v in atoms:
       lo=max(F(0),xa-y,y-xb);hi=max(abs(y-xa),abs(y-xb))
       assert hi+eta<r or lo-eta>r
   assert F(1,5)/(F(7,40)+F(1,1000))>1
   assert F(1,5)/(F(3,20)-F(1,1000))<2
   assert spacing/F(40)*F(3,100)*F(16,35)/F(1,4)==demand
   chain.append(dict(batch=batch,N=N,added_mass=str(eps/N),old_rerouting=str(moved),ratio=str(moved/(eps/N)),canonical_solver=solved,
     actual_spacing=str(spacing),actual_atom_mass=str(w),actual_D=D,actual_box_demand=str(demand),actual_rerouting=str(demand*F(N-1,2)),actual_trajectory_checks=boxchecks))
 boundary=[]
 for batch,ks in enumerate(((12,16,20),(24,32,48),(64,96,128))):
  for k in ks:
   rho=F(1,2**k);Q=[(F(-1,8),F(1,8),F(1)),(F(10),F(11),F(3,4)-rho)];P=Q+[( -rho,rho,F(1,2))]
   assert sum((h*(b-a) for a,b,h in P),F(0))==1
   lower=F(0);checks=0
   for j in range(3,k-5):
    left=2**j*rho;right=2**(j+1)*rho;assert right<=F(1,32)
    lower+=(right-left)/(right+rho)
    for x in (left,(left+right)/2,right):
     mp=maximal(P,x);mq=maximal(Q,x);assert mq==1 and mp==1+rho/(2*(x+rho))<F(9,8)
     mr=maximal([(-rho,rho,F(1,2))],x)
     assert mp-mq==mr and mp-1<=min(F(1,8),2*mr)
     h=(1+mp)/2
     trace=[sum((height*max(F(0),min(b,x+F(4,2**z))-max(a,x-F(4,2**z))) for a,b,height in P),F(0))/(2*F(4,2**z)) for z in range(1,7)]
     assert trace[0]<=h/8 and trace[-1]>=1>h/2 and h<mp<=2*h
     intervals=rec.intervals(trace,h/4,h/2)
     density=sum(((hi-lo)*trace[z-1] for z,(lo,hi) in intervals.items()),F(0))/(h/4)
     assert F(3,8)*h<=density<=F(3,4)*h
     checks+=1
   assert lower>=F(8,17)*(k-8)
   # Smooth-envelope bounds use only total mass and exact plateau containment.
   assert (F(3,8)+rho)/4<F(1,8) and F(1)+F(1,14)<F(9,8)
   boundary.append(dict(batch=batch,k=k,rho=str(rho),annuli=k-8,center_checks=checks,
      lost_volume_over_rho_lower=str(lower),simpler_lower=str(F(8,17)*(k-8)),
      lost_task_over_rho_lower=str(F(3,8)*lower)))
 deps=json.loads((ROOT/'rounds/060/verification.json').read_text())['sha256'];deps['rounds/061/verify.py']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
 for path,h in deps.items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h
 print(json.dumps(dict(round=61,status='passed',python=sys.version.split()[0],actual_inputs=records,chain=chain,boundary=boundary,
   main_theorem_proved=False,scope='Exact one-dimensional diagnostics; smooth counterexample proved analytically, not numerically integrated.',sha256=deps),indent=2,ensure_ascii=False))
if __name__=='__main__':main()
