from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict,deque
import importlib.util,json,argparse,sys,time
sys.dont_write_bytecode=True

def module(p,n):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def flow_at(groups,weights,C):
 sets=list(groups);G=len(sets);N=len(weights);sink=G+N+1;graph=[{} for _ in range(sink+1)];total=sum(groups.values(),F(0))
 def edge(a,b,v):graph[a][b]=graph[a].get(b,F(0))+v;graph[b].setdefault(a,F(0))
 for a,S in enumerate(sets,1):
  edge(0,a,groups[S])
  for i in S:edge(a,G+1+i,total+1)
 for i,w in enumerate(weights):edge(G+1+i,sink,C*w)
 value=F(0)
 while True:
  parent={0:None};queue=deque([0])
  while queue and sink not in parent:
   a=queue.popleft()
   for b,c in graph[a].items():
    if c>0 and b not in parent:parent[b]=a;queue.append(b)
  if sink not in parent:break
  path=[];b=sink
  while parent[b] is not None:a=parent[b];path.append((a,b));b=a
  delta=min(graph[a][b] for a,b in path)
  for a,b in path:graph[a][b]-=delta;graph[b][a]+=delta
  value+=delta
 reachable=set(parent);U=frozenset(i for i in range(N) if G+1+i in reachable)
 allocations={S:{i:graph[G+1+i][a] for i in S} for a,S in enumerate(sets,1)}
 return value,U,allocations

def solve(groups,weights):
 total=sum(groups.values(),F(0));C=total/sum(weights,F(0));U=frozenset(range(len(weights)));history=[]
 for iteration in range(10000):
  value,cut,allocation=flow_at(groups,weights,C)
  if value==total:
   assert sum((v for S,row in allocation.items() for v in row.values()),F(0))==total
   for S,row in allocation.items():assert sum(row.values(),F(0))==groups[S]
   loads=[sum((row.get(i,F(0)) for row in allocation.values()),F(0)) for i in range(len(weights))]
   assert all(loads[i]<=C*weights[i] for i in range(len(weights)))
   witness=sum((d for S,d in groups.items() if S<=U),F(0));wu=sum((weights[i] for i in U),F(0))
   assert witness==C*wu
   return C,U,history,allocation,loads
  demand=sum((d for S,d in groups.items() if S<=cut),F(0));mass=sum((weights[i] for i in cut),F(0));assert mass>0
  new=demand/mass;assert new>C
  history.append(dict(C=str(C),max_flow=str(value),cut=list(sorted(cut)),cut_demand=str(demand),cut_mass=str(mass),new_C=str(new)))
  C,U=new,cut
 raise AssertionError('ratio iteration bound exceeded')

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();root=args.repo_root
 flow=module(root/'rounds/047/flow_probe.py','flow56a');record=module(root/'rounds/046/threshold_probe.py','rec56a');prefix=module(root/'rounds/046/prefix_probe.py','pref56a');lag=module(root/'rounds/054/lag_probe.py','lag56a');obs=module(root/'rounds/053/obstruction_probe.py','obs56a')
 old=record.load(root,args.output);inputs=[c for c in flow.cases(old,prefix) if c[1] in lag.ORIGINAL|lag.EXTRA]
 inputs.append((0,'round53_four_atom_theta1',sorted(obs.ATOMS),F(1),{}));results=[]
 for batch,label,atoms,alpha,params in inputs:
  saved,radii,cells,_=flow.context(old,atoms,alpha);I=alpha/4;groups=defaultdict(F);natural=[F(0) for _ in atoms]
  for a,b,J,trace in cells:
   rec=record.intervals(trace,alpha/4,alpha/2);assert sum((hi-lo for lo,hi in rec.values()),F(0))==I
   x=(a+b)/2
   for k,(lo,hi) in rec.items():
    S=frozenset(i for i,(y,w) in enumerate(atoms) if abs(y-x)<radii[k]);q=sum((atoms[i][1] for i in S),F(0));g=q/(2*radii[k]);d=(b-a)*(hi-lo)/I*g
    assert d>0;groups[S]+=d
    # Strict full maximal witnesses must occur inside the selected ball.
    candidates=[dist for dist in {abs(y-x) for y,w in atoms} if 0<dist<radii[k]]
    assert candidates and max(sum((w for y,w in atoms if abs(y-x)<=dist),F(0))/(2*dist) for dist in candidates)>alpha
    assert max(sum((atoms[i][1] for i in S if abs(atoms[i][0]-x)<=dist),F(0))/(2*dist) for dist in candidates)>alpha
    for i in S:natural[i]+=d/q
  weights=[w for y,w in atoms];M=sum(groups.values(),F(0));X=alpha*F(saved['eligible_volume'])
  assert F(3,8)*X<=M<=F(3,4)*X
  C,U,history,allocation,loads=solve(groups,weights)
  assert C<=max(natural) and sum((weights[i]*natural[i] for i in range(len(weights))),F(0))==M
  # Independent exhaustive Hall optimum for small source counts.
  exhaustive=None
  if len(weights)<=10:
   exhaustive=max(sum((d for S,d in groups.items() if S<=frozenset(i for i in range(len(weights)) if mask>>i&1)),F(0))/sum((w for i,w in enumerate(weights) if mask>>i&1),F(0)) for mask in range(1,2**len(weights)))
   assert exhaustive==C
  chains=[]
  for S in sorted(groups,key=lambda S:(-len(S),tuple(sorted(S)))):
   for chain in chains:
    if S<=chain[-1]:chain.append(S);break
   else:chains.append([S])
  root_load=[sum(i in chain[0] for chain in chains) for i in range(len(weights))]
  kappa=max(root_load);chain_ratios=[]
  for chain in chains:
   for S in chain:
    tail=sum((groups[T] for T in chain if T<=S),F(0));mass=sum((weights[i] for i in S),F(0))
    assert tail<=F(3,2)*mass
    chain_ratios.append(tail/mass)
  assert C<=F(3,2)*kappa
  case=dict(batch=batch,label=label,N=len(atoms),D=saved['D'],groups=len(groups),X=str(X),M=str(M),M_over_X=str(M/X),optimal_C=str(C),natural_max=str(max(natural)),natural_over_optimal=str(max(natural)/C),tight_source_set=sorted(U),tight_source_mass=str(sum((weights[i] for i in U),F(0))),chain_cover=[list(map(sorted,chain)) for chain in chains],chain_root_congestion=kappa,chain_max_tail_ratio=str(max(chain_ratios)),iterations=len(history),history=history,exhaustive_crosscheck=exhaustive is not None,primal_rows=[dict(sources=sorted(S),demand=str(groups[S]),allocations={str(i):str(v) for i,v in row.items()}) for S,row in allocation.items()])
  results.append(case);print(label,float(C),float(max(natural)),len(history),flush=True)
 out=dict(status='passed',cases=results,batch_counts=[sum(c['batch']==i for c in results) for i in range(3)],scope='Finite full-E atomic shared-threshold allocation, exact rational primal and cut dual certificate; no all-input bound.')
 args.output.write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
