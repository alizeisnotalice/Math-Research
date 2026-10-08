#!/usr/bin/env python3
"""Canonical exact source profiles and frozen-child allocation audit."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict,deque
import importlib.util,json,sys,tempfile,hashlib,random
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def module(path,name):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def profile(groups,weights,alloc):
 """Repeated maximal tight cut, with exact primal/KKT certificate."""
 remaining=set(range(len(weights)));pending={S:S for S in groups};routes={};densities=[None]*len(weights);levels=[]
 while pending:
  order=sorted(remaining);ix={i:j for j,i in enumerate(order)};merged=defaultdict(F);members=defaultdict(list)
  for S,A in pending.items():
   key=frozenset(ix[i] for i in A);merged[key]+=groups[S];members[key].append(S)
  ws=[weights[i] for i in order];C,U,_,rows,loads=alloc.solve(merged,ws)
  # A positive allocated source i can reroute its task to any j in S.
  reverse=[set() for i in order]
  for S,row in rows.items():
   for i,v in row.items():
    if v>0:
     for j in S:reverse[j].add(i)
  reachable={i for i,w in enumerate(ws) if loads[i]<C*w};queue=deque(reachable)
  while queue:
   j=queue.popleft()
   for i in reverse[j]:
    if i not in reachable:reachable.add(i);queue.append(i)
  T=frozenset(set(range(len(ws)))-reachable);assert T
  demand=sum((d for S,d in merged.items() if S<=T),F(0));mass=sum((ws[i] for i in T),F(0));assert demand==C*mass
  assert all(rows[S].get(i,F(0))==0 for S in rows if not S<=T for i in T)
  original_T=frozenset(order[i] for i in T)
  if levels:assert C<F(levels[-1]['density'])
  levels.append(dict(sources=sorted(original_T),density=str(C),mass=str(mass),demand=str(demand)))
  finished=[]
  for S,row in rows.items():
   if not S<=T:continue
   for old in members[S]:
    routes[old]={order[i]:v*groups[old]/merged[S] for i,v in row.items() if v}
    finished.append(old)
  for old in finished:del pending[old]
  for i in original_T:densities[i]=C
  remaining-=original_T
  pending={S:A-original_T for S,A in pending.items()}
  assert all(pending.values())
 for i in remaining:densities[i]=F(0)
 for S,d in groups.items():
  assert sum(routes[S].values(),F(0))==d
  minimum=min(densities[i] for i in S)
  assert all(i in S and v>=0 and (v==0 or densities[i]==minimum) for i,v in routes[S].items())
 assert all(sum((row.get(i,F(0)) for row in routes.values()),F(0))==weights[i]*densities[i] for i in range(len(weights)))
 # Tightness of every superlevel, and exact convex-energy identity.
 for t in sorted(set([F(0)]+densities)):
  U=frozenset(i for i,x in enumerate(densities) if x>t)
  assert sum((weights[i]*densities[i] for i in U),F(0))==sum((d for S,d in groups.items() if S<=U),F(0))
 energy=sum((w*x*x for w,x in zip(weights,densities)),F(0))
 assert energy==sum((d*min(densities[i] for i in S) for S,d in groups.items()),F(0))
 return densities,routes,levels

def frozen_tree(groups,weights,alloc,mode):
 def recurse(U):
  local={S:d for S,d in groups.items() if S<=U}
  if len(U)==1:
   i=next(iter(U));d=local.get(U,F(0));return ({U:{i:d}} if d else {})
  order=sorted(U)
  if mode=='index':k=len(order)//2
  elif mode=='mass':
   total=sum((weights[i] for i in order),F(0));k=min(range(1,len(order)),key=lambda k:abs(sum((weights[i] for i in order[:k]),F(0))-total/2))
  else:order.sort(key=lambda i:(-weights[i],i));k=1
  children=[frozenset(order[:k]),frozenset(order[k:])];child=[recurse(A) for A in children]
  output={S:row.copy() for rows in child for S,row in rows.items()}
  crossing={S:d for S,d in local.items() if not any(S<=A for A in children)}
  order=sorted(U);ix={i:j for j,i in enumerate(order)};augmented={}
  for i in order:
   used=sum((row.get(i,F(0)) for row in output.values()),F(0))
   if used:augmented[frozenset({ix[i]})]=used
  for S,d in crossing.items():augmented[frozenset(ix[i] for i in S)]=d
  if augmented:
   phi,rows,_=profile(augmented,[weights[i] for i in order],alloc)
   for S in crossing:output[S]={order[i]:v for i,v in rows[frozenset(ix[i] for i in S)].items()}
   assert all(sum((row.get(i,F(0)) for row in output.values()),F(0))==weights[i]*phi[ix[i]] for i in order)
  assert set(output)==set(local)
  return output
 rows=recurse(frozenset(range(len(weights))))
 phi=[sum((row.get(i,F(0)) for row in rows.values()),F(0))/w for i,w in enumerate(weights)]
 for S,d in groups.items():assert set(rows[S])<=S and sum(rows[S].values(),F(0))==d
 return phi,rows

def main():
 alloc=module(ROOT/'rounds/056/allocation_probe.py','alloc60');previous=module(ROOT/'rounds/058/verify.py','prev60')
 flow=module(ROOT/'rounds/047/flow_probe.py','flow60');rec=module(ROOT/'rounds/046/threshold_probe.py','rec60')
 prefix=module(ROOT/'rounds/046/prefix_probe.py','prefix60');lag=module(ROOT/'rounds/054/lag_probe.py','lag60');obs=module(ROOT/'rounds/053/obstruction_probe.py','obs60')
 rows=[];witness=None;height_checks=[]
 with tempfile.TemporaryDirectory(prefix='euclidean60_') as tmp:
  old=rec.load(ROOT,Path(tmp)/'scratch.json')
  inputs=[c for c in flow.cases(old,prefix) if c[1] in lag.ORIGINAL|lag.EXTRA]
  inputs.append((0,'round53_four_atom_theta1',sorted(obs.ATOMS),F(1),{}))
  for L in (4,6,10,12):inputs.append((2,f'chain_L{L}_alpha3/2',old.old.old.chain(L,den=4,position=F(9,8),weight=F(27,16)),F(3,2),{}))
  for batch,nnear in enumerate((3,4,6)):
   for serial in range(30):
    seed=600100+batch*100+serial;rng=random.Random(seed)
    positions=sorted(rng.sample(range(-256,257),nnear));raw=[rng.randint(1,32) for _ in positions]
    atoms=[(F(y,512),F(3,8)*F(w,sum(raw))) for y,w in zip(positions,raw)]+[(F(10),F(5,8))]
    alpha=rng.choice((F(1,2),F(1),F(2),F(4)))
    inputs.append((batch,f'random_split_s{seed}',atoms,alpha,{}))
  for batch,label,atoms,alpha,_ in inputs:
   saved,radii,groups,tasks=previous.capture_data(flow,rec,old,atoms,alpha);weights=[w for y,w in atoms]
   if not groups:continue
   C,_,_,_,_=alloc.solve(groups,weights);phi,routes,levels=profile(groups,weights,alloc);assert max(phi)==C
   for index in sorted({i*(len(tasks)-1)//5 for i in range(6)}):
    x,J,tr,k,lo,hi,S=tasks[index];m=previous.maximal(atoms,x);u=tr[0];z=max(tr)
    a=alpha*F(3,4);b=alpha*F(3,2);LP=max(m/2,8*u);UP=min(m,2*z)
    original=max(F(0),min(b,UP)-max(a,LP))
    relaxed=max(F(0),min(2*b,m,4*z)-max(a/2,m/2,4*u))
    for eps in (F(1,4),F(1,64),F(1,4096)):
     UQ=(1-eps)*UP
     lost=max(F(0),min(b,UP)-max(a,LP,UQ))
     assert lost<=eps*m and lost<=2*eps*relaxed
     height_checks.append((batch,lost>0))

   frozen={}
   for mode in ('index','mass','peel'):
    fp,fr=frozen_tree(groups,weights,alloc,mode);cap=max(fp);assert cap>=C
    # Every stop-loss functional is minimized by the canonical profile.
    for t in sorted(set(phi+fp+[F(0)])):
     assert sum((w*max(F(0),x-t) for w,x in zip(weights,phi)),F(0))<=sum((w*max(F(0),x-t) for w,x in zip(weights,fp)),F(0))
    frozen[mode]=dict(C=str(cap),over_optimal=str(cap/C))
    if cap>C and (witness is None or len(atoms)<len(witness['atoms'])):
     witness=dict(label=label,alpha=str(alpha),D=saved['D'],mode=mode,atoms=[dict(x=str(x),w=str(w)) for x,w in atoms],
                  groups=[dict(S=sorted(S),d=str(d)) for S,d in groups.items()],optimal_C=str(C),optimal_profile=list(map(str,phi)),frozen_C=str(cap),frozen_profile=list(map(str,fp)),
                  optimal_routes=[dict(S=sorted(S),row={str(i):str(v) for i,v in row.items()}) for S,row in routes.items()],
                  frozen_routes=[dict(S=sorted(S),row={str(i):str(v) for i,v in row.items()}) for S,row in fr.items()])
   rows.append(dict(batch=batch,label=label,N=len(atoms),D=saved['D'],C=str(C),levels=levels,frozen=frozen))
 abstract=[]
 for batch,powers in enumerate(((2,4,6),(8,12,16),(24,32,48))):
  for power in powers:
   eps=F(1,2**power);W=2+eps;weights=[1/W,1/W,eps/W]
   gs={frozenset({0,1}):1/W,frozenset({0,2}):(1+eps)/W}
   phi,rt,_=profile(gs,weights,alloc);assert phi==[F(1)]*3
   child={frozenset({0,1}):1/W}
   cphi,_,_=profile(child,weights,alloc);assert cphi==[F(1,2),F(1,2),F(0)]
   augmented={frozenset({0}):F(1,2)/W,frozenset({1}):F(1,2)/W,frozenset({0,2}):(1+eps)/W}
   fphi,_,_=profile(augmented,weights,alloc);cap=(F(3,2)+eps)/(1+eps)
   assert fphi==[cap,F(1,2),cap]
   abstract.append(dict(batch=batch,epsilon=str(eps),global_C='1',frozen_C=str(cap)))
 recursive=[];recursive_concrete=[]
 for batch,Ls in enumerate(((1,2,4),(8,12,16),(24,32,48))):
  for L in Ls:
   eps=F(1,4*L*L);N=2**L;t=F(1,2)
   assert F(N,2)+sum((F(2**(L-j),2) for j in range(1,L)),F(0))+1==N
   for j in range(1,L+1):
    mass=2**(L-j);d=F(mass,2) if j<L else F(1)
    new=(t*mass+eps+d)/(mass+eps)
    assert new>=max(t,F(1)) and mass*(new-t)+eps*(new-1)==d
    t=new
   assert F(L+2,2)-F(1,8)<=t<=F(L+2,2)
   recursive.append(dict(batch=batch,L=L,epsilon=str(eps),global_C='1',frozen_C=str(t),lower_bound=str(F(L+2,2)-F(1,8)),scope='Abstract allocation only'))
   if L in (1,2,4):
    weights=[F(1)]*N;gs={frozenset(range(N)):F(N,2)};fphi=[F(1,2)]*N
    for j in range(1,L+1):
     weights.append(eps);idx=N+j-1;S=frozenset(range(2**(L-j)))|{idx};d=F(2**(L-j),2) if j<L else F(1)
     augmented={frozenset({i}):weights[i]*value for i,value in enumerate(fphi)}
     augmented[frozenset({idx})]=eps;augmented[S]=d
     fphi,_,_=profile(augmented,weights,alloc)
     gs[frozenset({idx})]=eps;gs[S]=d
    full,_,_=profile(gs,weights,alloc)
    assert max(fphi)==t and all(value==1 for value in full)
    recursive_concrete.append(dict(L=L,N=len(weights),global_C='1',frozen_C=str(t)))
 canonical=[]
 for batch,powers in enumerate(((6,8,10),(12,16,20),(24,32,48))):
  for power in powers:
   eps=F(1,2**power)
   for h in (F(3,4),F(4,5),F(9,10)):
    # On the full central interval all canonical balls j<=4 capture the near mass.
    pp=[h*(1+eps)*F(2**j,32) for j in range(1,5)]
    qq=[h*F(2**j,32) for j in range(1,5)]
    assert h<1<=1+eps<=2*h and pp[0]<=h/8
    assert next(j for j,g in enumerate(pp,1) if g>h/2)==4 and max(qq)==h/2
   canonical.append(dict(batch=batch,epsilon=str(eps),height_window=['3/4','9/10'],D=4,weighted_lost_volume='3/80'))
 # Absence of a witness is a valid finite-search outcome, not a theorem.
 deps=sorted((ROOT/'runtime').rglob('*.py'))+[ROOT/'rounds'/s for s in
 ('041/square_probe.py','042/row_probe.py','043/global_probe.py','044/tail_probe.py','046/prefix_probe.py','046/threshold_probe.py','047/flow_probe.py','053/obstruction_probe.py','054/lag_probe.py','056/allocation_probe.py','058/verify.py')]+[Path(__file__).resolve()]
 print(json.dumps(dict(status='passed',round=60,date='2026-10-08',python=sys.version.split()[0],cases=rows,actual_frozen_counterexample=witness,abstract_frozen_family=abstract,recursive_frozen_family=recursive,recursive_concrete=recursive_concrete,fixed_grid_height_checks=dict(count=len(height_checks),positive_loss=sum(v for b,v in height_checks),batch_counts=[sum(b==i for b,v in height_checks) for i in range(3)]),canonical_cutoff_counterexample=canonical,
 scope='Exact finite n=1 full-E source profiles, finite convex optimality and frozen-child routing. No uniform high-dimensional budget.',main_weak_type_theorem_proved=False,
 sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in deps}),ensure_ascii=False,indent=2))
if __name__=='__main__':main()
