#!/usr/bin/env python3
"""Shared-threshold stability, original-task tree gluing, and band instability."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import importlib.util,sys,json,random,tempfile,hashlib
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def module(p,n):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def trace_check(p,q,a=F(1,4),b=F(1,2)):
 assert len(p)==len(q) and p[0]<=a and q[0]<=a and max(q)>b
 assert all(0<=y<=x for x,y in zip(p,q))
 assert all(p[j]<=2*p[j-1] and q[j]<=2*q[j-1] for j in range(1,len(p)))
 delta=max(x-y for x,y in zip(p,q));r=[F(0)];s=[F(0)]
 for x,y in zip(p,q):r.append(max(r[-1],x));s.append(max(s[-1],y))
 clip=lambda x:max(a,min(b,x))
 c=list(map(clip,r));d=list(map(clip,s));e=[x-y for x,y in zip(c,d)]
 first=sum(((r[j]-s[j])*(c[j]-c[j-1]) for j in range(1,len(r))),F(0))
 second=sum((e[j]*(s[j+1]-s[j]) for j in range(len(r)-1)),F(0))
 cuts=sorted({a,b}|{x for x in p+q if a<x<b});signed=absolute=removed=F(0);jump=0
 for lo,hi in zip(cuts,cuts[1:]):
  beta=(lo+hi)/2;j=next(j for j,x in enumerate(p) if x>beta);k=next(k for k,y in enumerate(q) if y>beta)
  assert j<=k and p[j]-q[k]<=p[j]-q[j]
  signed+=(hi-lo)*(p[j]-q[k]);absolute+=(hi-lo)*abs(p[j]-q[k]);removed+=(hi-lo)*(p[j]-q[j]);jump=max(jump,k-j)
 assert signed==first-second and 0<=first<=delta*(b-a) and 0<=second<=2*b*delta
 assert signed<=removed and absolute<=((b-a)+2*b)*delta
 return dict(D=len(p),delta=str(delta),signed=str(signed),absolute=str(absolute),original_selected_removed=str(removed),largest_jump=jump)

def glue(groups,weights,alloc,mode):
 N=len(weights);operations=0;nodes=0
 def recurse(U):
  nonlocal operations,nodes
  nodes+=1;local={S:d for S,d in groups.items() if S<=U}
  if len(U)==1:
   i=next(iter(U));d=local.get(U,F(0));return d/weights[i],({U:{i:d}} if d else {})
  ordered=sorted(U)
  if mode=='index':k=len(ordered)//2
  elif mode=='mass':
   total=sum((weights[i] for i in ordered),F(0));k=min(range(1,len(ordered)),key=lambda k:abs(sum((weights[i] for i in ordered[:k]),F(0))-total/2))
  else:
   ordered.sort(key=lambda i:(-weights[i],i));k=1
  children=[frozenset(ordered[:k]),frozenset(ordered[k:])]
  data=[recurse(V) for V in children];cross={S:d for S,d in local.items() if not any(S<=V for V in children)}
  ordered=sorted(U);index={i:j for j,i in enumerate(ordered)};reserve={}
  augmented=defaultdict(F)
  for V,(cv,_) in zip(children,data):
   for i in V:reserve[i]=cv;augmented[frozenset({index[i]})]+=cv*weights[i]
  for S,d in cross.items():augmented[frozenset(index[i] for i in S)]+=d
  augmented={S:d for S,d in augmented.items() if d}
  if augmented:
   cv,cut,_,routes,_=alloc.solve(augmented,[weights[i] for i in ordered]);operations+=1
  else:cv=F(0);routes={}
  output={S:row.copy() for _,rows in data for S,row in rows.items()}
  for S,d in cross.items():
   row=routes[frozenset(index[i] for i in S)]
   output[S]={i:row.get(index[i],F(0)) for i in S}
  for i in U:
   incoming=sum((output[S].get(i,F(0)) for S in cross),F(0))
   assert incoming<=(cv-reserve[i])*weights[i]
   assert sum((row.get(i,F(0)) for row in output.values()),F(0))<=cv*weights[i]
  assert set(output)==set(local)
  assert all(sum(output[S].values(),F(0))==d for S,d in local.items())
  return cv,output
 C,rows=recurse(frozenset(range(N)))
 return C,operations,nodes

def main():
 flow=module(ROOT/'rounds/047/flow_probe.py','flow59');rec=module(ROOT/'rounds/046/threshold_probe.py','rec59')
 prefix=module(ROOT/'rounds/046/prefix_probe.py','pre59');lag=module(ROOT/'rounds/054/lag_probe.py','lag59')
 obs=module(ROOT/'rounds/053/obstruction_probe.py','obs59');alloc=module(ROOT/'rounds/056/allocation_probe.py','alloc59')
 previous=module(ROOT/'rounds/058/verify.py','v58')
 trace_rows=[]
 for batch,Ds in enumerate(((8,12,16),(24,40,64),(96,160,256))):
  for case,D in enumerate(Ds):
   rng=random.Random(590100+10*batch+case);eps=F(1,2**(batch+case+5));q=[F(1,8)];loss=[eps]
   for j in range(1,D-4):
    top=min(512,int(2*q[-1]*1024));q.append(F(rng.randint(128,top),1024))
    loss.append(eps*F(rng.randint(128,256),256))
   while len(q)<D:q.append(min(F(3,4),2*q[-1]));loss.append(eps)
   p=[a+b for a,b in zip(q,loss)]
   result=trace_check(p,q);result.update(batch=batch,seed=590100+10*batch+case);trace_rows.append(result)
 # Explicit arbitrarily long clock jump with small integrated cost.
 long_rows=[]
 for batch,Ls in enumerate(((4,8,16),(32,64,128),(256,512,1024))):
  for L in Ls:
   eps=F(1,2**(batch+8));q=[F(1,8),F(1,4)]+[F(3,8)-eps]*L+[F(3,4)-2*eps]
   p=[x+eps for x in q];row=trace_check(p,q);assert row['largest_jump']==L
   row.update(batch=batch,L=L);long_rows.append(row)
 with tempfile.TemporaryDirectory(prefix='euclidean59_') as tmp:
  old=rec.load(ROOT,Path(tmp)/'scratch.json')
  inputs=[c for c in flow.cases(old,prefix) if c[1] in lag.ORIGINAL|lag.EXTRA]
  inputs.append((0,'round53_four_atom_theta1',sorted(obs.ATOMS),F(1),{}))
  for L in (4,6,10,12):inputs.append((2,f'chain_L{L}_alpha3/2',old.old.old.chain(L,den=4,position=F(9,8),weight=F(27,16)),F(3,2),{}))
  cases=[]
  for batch,label,atoms,alpha,_ in inputs:
   saved,radii,groups,tasks=previous.capture_data(flow,rec,old,atoms,alpha)
   weights=[w for y,w in atoms];exact,_,_,_,_=alloc.solve(groups,weights)
   trees={}
   for mode in ('index','mass','peel'):
    c,operations,nodes=glue(groups,weights,alloc,mode);assert c>=exact
    trees[mode]=dict(C=str(c),over_optimal=str(c/exact),flow_solves=operations,nodes=nodes)
   cases.append(dict(batch=batch,label=label,N=len(atoms),D=saved['D'],optimal_C=str(exact),trees=trees))
 plateau=[];pair=[]
 for batch,powers in enumerate(((6,8,10),(12,16,20),(24,32,48))):
  for power in powers:
   eps=F(1,2**power);p=[(1+eps)/16,(1+eps)/8,(1+eps)/4,(1+eps)/2]
   intervals=rec.intervals(p,F(1,4),F(1,2))
   demand=F(1,4)*sum(((hi-lo)*p[k-1] for k,(lo,hi) in intervals.items()),F(0))/F(1,4)
   assert demand==(1+eps)*(2-eps)/16 and next(i for i,g in enumerate(p,1) if g>F(1,2))==4
   deleted=eps/4;fraction=eps/(1+eps)
   plateau.append(dict(batch=batch,epsilon=str(eps),original_subtask_demand=str(demand),removed_mass=str(deleted),loss_over_removed=str(demand/deleted),selected_removed_fraction=str(fraction)))
   groups={frozenset({0}):F(1,32)-eps/2,frozenset({1}):F(1,32)-eps/2,frozenset({0,1}):eps}
   c,_,_=glue(groups,[F(1,8),F(1,8)],alloc,'index');assert c==F(1,4)
   pair.append(dict(batch=batch,epsilon=str(eps),leaf_C=str(F(1,4)-4*eps),root_C=str(c),increment=str(4*eps)))
 # Abstract unequal-slack example, no Euclidean realization claimed.
 gs={frozenset({0}):F(1,2),frozenset({0,1}):F(1,2)}
 c,_,_=glue(gs,[F(1,2),F(1,2)],alloc,'index');assert c==1
 abstract_growth=[]
 for batch,Ls in enumerate(((1,2,4),(8,16,32),(64,128,256))):
  for L in Ls:
   N=(L+1)*2**L
   # Every height contributes total demand 1/(L+1); disjoint source colors.
   assert sum((F(1,L+1) for h in range(L+1)),F(0))==1
   abstract_growth.append(dict(batch=batch,L=L,N=N,optimal_C='1',uniform_reservation_C=str(L+1),scope='Abstract task system only, no Euclidean realization'))
 concrete_growth=[]
 for L in (1,2,3):
  N=(L+1)*2**L;weights=[F(1,N)]*N;gs={}
  for k in range(2**L):gs[frozenset({k*(L+1)})]=F(1,N)
  for h in range(1,L+1):
   for start in range(0,2**L,2**h):
    S=frozenset(k*(L+1)+h for k in range(start,start+2**h));gs[S]=F(len(S),N)
  assert all(sum(i in S for S in gs)==1 for i in range(N))
  exact,_,_,_,_=alloc.solve(gs,weights);bound,ops,nodes=glue(gs,weights,alloc,'index')
  assert exact==1 and bound==L+1
  concrete_growth.append(dict(L=L,N=N,optimal_C=str(exact),uniform_reservation_C=str(bound),flow_solves=ops))
 deps=sorted((ROOT/'runtime').rglob('*.py'))+[ROOT/'rounds'/s for s in
  ('041/square_probe.py','042/row_probe.py','043/global_probe.py','044/tail_probe.py','046/prefix_probe.py','046/threshold_probe.py','047/flow_probe.py','053/obstruction_probe.py','054/lag_probe.py','056/allocation_probe.py','058/verify.py')]+[Path(__file__).resolve()]
 print(json.dumps(dict(status='passed',round=59,date='2026-10-08',python=sys.version.split()[0],
  trace_tests=trace_rows,long_jump_tests=long_rows,actual_tree_cases=cases,
  plateau_counterexample=plateau,weak_connection_tree_checks=pair,
  abstract_unequal_slack=dict(root_C='1',crude_bound='3/2'),abstract_reservation_growth=abstract_growth,concrete_abstract_growth=concrete_growth,
  scope='Exact rational scalar diagnostics and full-E n=1 atomic original-task routing. Plateau subtask demand exact; smooth analogue proved analytically only. No uniform root potential or dimension-free weak theorem established.',
  main_weak_type_theorem_proved=False,
  sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in deps}),ensure_ascii=False,indent=2))
if __name__=='__main__':main()
