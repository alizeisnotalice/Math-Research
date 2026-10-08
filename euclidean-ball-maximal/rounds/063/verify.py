#!/usr/bin/env python3
"""Full actual task completion and exact conditional minimum-rerouting audit."""
from pathlib import Path
from fractions import Fraction as F
import importlib.util,json,sys,tempfile,hashlib,random
from collections import defaultdict
from bisect import bisect_left,bisect_right
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[2]
def module(path,name):
 s=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
from collections import deque
from fractions import Fraction as F

def minmove(groups,B,weights,C,extra_source,extra):
 keys=list(groups);G=len(keys);N=len(weights);V=G+N+2;start=V-2;end=V-1;adj=[[] for _ in range(V)]
 def edge(u,v,capacity,cost,initial=F(0)):
  a=[v,capacity-initial,cost,len(adj[v])];b=[u,initial,-cost,len(adj[u])];adj[u].append(a);adj[v].append(b);return a,capacity
 handles={};load=[F(0)]*N
 for t,S in enumerate(keys):
  hs={}
  for i in S:
   old=B[S].get(i,F(0));load[i]+=old
   hs[i]=(edge(t,G+i,old,0,old),edge(t,G+i,groups[S],1))
  handles[S]=hs
 for i in range(N):edge(G+i,end,C*weights[i],0,load[i])
 edge(start,G+extra_source,extra,0)
 remaining=extra;cost=F(0);iterations=0;max_distance=0
 while remaining:
  ds=[None]*V;ds[start]=0;parent=[None]*V;q=deque([start]);inq={start}
  while q:
   u=q.popleft();inq.remove(u)
   for idx,e in enumerate(adj[u]):
    v,cap,c,rev=e
    if cap and (ds[v] is None or ds[v]>ds[u]+c):
     ds[v]=ds[u]+c;parent[v]=(u,idx)
     if v not in inq:q.append(v);inq.add(v)
  assert ds[end] is not None
  amt=remaining;v=end
  while v!=start:
   u,idx=parent[v];amt=min(amt,adj[u][idx][1]);v=u
  v=end
  while v!=start:
   u,idx=parent[v];e=adj[u][idx];e[1]-=amt;adj[v][e[3]][1]+=amt;v=u
  remaining-=amt;cost+=amt*ds[end];iterations+=1;max_distance=max(max_distance,ds[end])
 # A global potential certifies no negative residual cycle, hence exact optimum.
 potential=[0]*V
 for step in range(V):
  changed=False
  for u in range(V):
   for v,cap,c,rev in adj[u]:
    if cap and potential[v]>potential[u]+c:potential[v]=potential[u]+c;changed=True
  if not changed:break
 else:raise AssertionError('negative residual cycle')
 after=[F(0)]*N;movement=F(0)
 for S,hs in handles.items():
  row={i:sum((capacity-e[1] for e,capacity in pairs),F(0)) for i,pairs in hs.items()}
  assert sum(row.values())==groups[S] and all(v>=0 for v in row.values())
  for i,v in row.items():after[i]+=v;movement+=abs(v-B[S].get(i,F(0)))/2
 after[extra_source]+=extra
 assert all(v<=C*w for v,w in zip(after,weights)) and movement==cost
 return dict(cost=str(cost),ratio=str(cost/extra),iterations=iterations,max_path_cost=max_distance)

def independent_lattice(N,s,w,D):
 radii={j:F(4,2**j) for j in range(1,D+1)};positions=[i*s for i in range(N)]
 edges=sorted({y+sign*r for y in positions for r in radii.values() for sign in (-1,1)})
 groups=defaultdict(F);nearest=[F(0)]*N
 for i,y in enumerate(positions):
  for sign in (-1,1):
   a,b=sorted((y+sign*s/10,y+sign*s/5));nodes=[a]+edges[bisect_right(edges,a):bisect_left(edges,b)]+[b]
   for l,h in zip(nodes,nodes[1:]):
    x=(l+h)/2;trace=[];members=[]
    for j,r in radii.items():
     lo=bisect_right(positions,x-r);hi=bisect_left(positions,x+r);S=frozenset(range(lo,hi));members.append(S);trace.append(w*len(S)/(2*r))
    assert trace[0]<=F(1,8) and max(trace)>F(1,2)
    cuts=sorted({F(1,4),F(1,2)}|{v for v in trace if F(1,4)<v<F(1,2)})
    for lo,hi in zip(cuts,cuts[1:]):
     beta=(lo+hi)/2;k=next(j for j,g in enumerate(trace) if g>beta);assert i in members[k]
     d=(h-l)*(hi-lo)*4*trace[k];groups[members[k]]+=d;nearest[i]+=d
 assert all(v<=F(328,1225)*w for v in nearest)
 assert nearest[0]==nearest[-1]==F(328,1225)*w
 return groups,nearest

def winners(oldruntime,flow,rec,atoms,alpha):
 saved,radii,cells,_=flow.context(oldruntime,atoms,alpha);loads=[F(0)]*len(atoms);collective=F(0);total=F(0)
 boundaries=sorted({y+sign*w/(2*alpha) for y,w in atoms for sign in (-1,1)})
 for a,b,J,trace in cells:
  nodes=[a]+boundaries[bisect_right(boundaries,a):bisect_left(boundaries,b)]+[b]
  for l,h in zip(nodes,nodes[1:]):
   x=(l+h)/2;choices=[i for i,(y,w) in enumerate(atoms) if abs(y-x)<w/(2*alpha)]
   selected=min(choices) if choices else None
   for k,(lo,hi) in rec.intervals(trace,alpha/4,alpha/2).items():
    d=(h-l)*(hi-lo)/(alpha/4)*trace[k-1];total+=d
    if selected is None:collective+=d
    else:
     y,w=atoms[selected];assert abs(y-x)<radii[k] and abs(y-x)>=w/(4*alpha)
     loads[selected]+=d
 assert sum(loads)+collective==total
 assert all(v<=F(3,8)*w for v,(y,w) in zip(loads,atoms))
 return dict(total=str(total),single_witness=str(sum(loads)),collective=str(collective),maximum_atom_load=str(max((v/w for v,(y,w) in zip(loads,atoms)),default=F(0))))

def main():
 p=module('rounds/060/verify.py','profile63');alloc=module('rounds/056/allocation_probe.py','alloc63');prev=module('rounds/058/verify.py','prev63');flow=module('rounds/047/flow_probe.py','flow63');rec=module('rounds/046/threshold_probe.py','rec63')
 records=[];clusters=[];winner_checks=[]
 with tempfile.TemporaryDirectory(prefix='euclidean63_') as tmp:
  oldruntime=rec.load(ROOT,Path(tmp)/'scratch.json')
  for batch,Ls in enumerate(((1,2),(3,4),(5,6,7))):
   for L in Ls:
    N=2**L;s=F(4,7*N);w=2*s/5;atoms=[(i*s,w) for i in range(N)]+[(F(10),F(27,35))]
    saved,r,groups,tasks=prev.capture_data(flow,rec,oldruntime,atoms,F(1));weights=[v for y,v in atoms]
    direct,nearest=independent_lattice(N,s,w,L+5);assert direct==groups
    winner=winners(oldruntime,flow,rec,atoms,F(1));assert F(winner['collective'])==0
    winner_checks.append(dict(batch=batch,label=f'lattice{N}',**winner))
    phi,T,levels=p.profile(groups,weights,alloc);C=max(phi)
    demand=6*s/4375;before=dict(groups);before[frozenset({0})]-=demand;assert before[frozenset({0})]>=0
    bphi,B,_=p.profile(before,weights,alloc)
    assert all(v==C for v in phi[:N]) and phi[N]==0
    assert all(v==C-demand/(N*w) for v in bphi[:N]) and bphi[N]==0
    for S,d in groups.items():assert d<=sum((weights[i] for i in S),F(0))
    result=minmove(before,B,weights,C,0,demand)
    elementary=1-F(1,N)+(F(1,8) if N>=8 else F(0))
    assert F(result['ratio'])>=elementary
    frozen=bphi.copy();frozen[0]+=demand/w;assert max(frozen)>C
    records.append(dict(batch=batch,N=N,D=saved['D'],capture_classes=len(groups),full_capacity=str(C),new_demand=str(demand),before_capacity=str(max(bphi)),frozen_capacity=str(max(frozen)),geometric_lower_ratio=str(elementary),subtask_chain_ratio=str(F(N-1,2)),**result))
  for batch,powers in enumerate(((7,8,9),(12,16,20),(24,32,48))):
   for power in powers:
    eta=F(1,2**power);atoms=[(-eta,F(1,8)),(eta,F(1,8)),(F(10),F(3,4))]
    saved,r,groups,tasks=prev.capture_data(flow,rec,oldruntime,atoms,F(1))
    assert groups=={frozenset({0,1}):F(1,16)}
    assert F(saved['eligible_volume'])==F(1,8)
    phi,rows,levels=p.profile(groups,[F(1,8),F(1,8),F(3,4)],alloc)
    assert phi==[F(1,4),F(1,4),F(0)]
    for x,J,trace,k,lo,hi,S in tasks:
     assert J==3 and k==4 and S==frozenset({0,1})
     assert prev.maximal(atoms,x)==F(1,8)/(abs(x)+eta)
     assert next(j for j,g in enumerate(trace,1) if g>F(1,2))==5
    # Smooth support envelopes; kappa=eta/4 and smaller eta if needed.
    e=eta/4;kap=e/4;diam=e+kap
    assert diam<F(1,128)
    assert F(1,8)+2*diam<F(1,4)
    assert F(1,4)/(2*(F(7,64)+diam))>1
    assert F(1,4)/(2*(F(5,64)-diam))<2
    assert F(7,64)+diam<F(1,8)
    winner_checks.append(dict(batch=batch,label=f'cluster_eta2minus{power}',**winners(oldruntime,flow,rec,atoms,F(1))))
    clusters.append(dict(batch=batch,eta=str(eta),mass_near='1/4',demand='1/16',capacity='1/4',singleton_demand='0',completion_surplus='0',smooth_eta=str(e),smooth_kappa=str(kap)))
  for batch,N in enumerate((3,5,8)):
   for serial in range(5):
    seed=630100+100*batch+serial;rng=random.Random(seed);xs=sorted(rng.sample(range(-512,513),N));ws=[rng.randint(1,64) for _ in xs]
    atoms=[(F(x,1024),F(3,8)*F(w,sum(ws))) for x,w in zip(xs,ws)]+[(F(10),F(5,8))]
    alpha=rng.choice((F(1,2),F(1),F(2),F(4)))
    winner_checks.append(dict(batch=batch,label=f'seed{seed}',**winners(oldruntime,flow,rec,atoms,alpha)))
 deps=json.loads((ROOT/'rounds/062/verification.json').read_text())['sha256'];deps['rounds/063/verify.py']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
 for path,h in deps.items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h
 print(json.dumps(dict(round=63,status='passed',python=sys.version.split()[0],full_completion=records,zero_surplus=clusters,winner_checks=winner_checks,sha256=deps,main_theorem_proved=False,scope='Full E in n=1. Minimum cost is conditional on the recorded deterministic old optimal routing; no uniform recourse conclusion.'),ensure_ascii=False,indent=2))
if __name__=='__main__':main()
