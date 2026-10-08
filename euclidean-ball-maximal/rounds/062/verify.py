#!/usr/bin/env python3
"""Exact capacity-slack rerouting certificates on actual first-exit demands."""
from pathlib import Path
from fractions import Fraction as F
import importlib.util,json,sys,tempfile,random,hashlib
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[2]
def module(path,name):
 s=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def linear(A,b):
 n=len(b);mat=[list(row)+[v] for row,v in zip(A,b)]
 for col in range(n):
  pivot=next(i for i in range(col,n) if mat[i][col]);mat[col],mat[pivot]=mat[pivot],mat[col]
  scale=mat[col][col];mat[col]=[v/scale for v in mat[col]]
  for i in range(n):
   if i!=col and mat[i][col]:
    factor=mat[i][col];mat[i]=[x-factor*y for x,y in zip(mat[i],mat[col])]
 return [mat[i][-1] for i in range(n)]

def reroute(old,new,B,T,weights,C,eps,trace=False):
 n=len(weights);theta=1-eps;b=[sum((row.get(i,F(0)) for row in B.values()),F(0)) for i in range(n)]
 target_old={S:{i:v*old.get(S,F(0))/(old.get(S,F(0))+new.get(S,F(0))) for i,v in row.items()} for S,row in T.items() if old.get(S,F(0))}
 target_new={S:{i:v*new.get(S,F(0))/(old.get(S,F(0))+new.get(S,F(0))) for i,v in row.items()} for S,row in T.items() if new.get(S,F(0))}
 a=[sum((row.get(i,F(0)) for row in target_new.values()),F(0)) for i in range(n)]
 K=[[F(0)]*n for _ in range(n)]
 for S,row in B.items():
  for i,v in row.items():
   if not v:continue
   for j,t in target_old[S].items():K[i][j]+=v/b[i]*t/old[S]
 def transfer(r):return [a[j]+sum((r[i]*K[i][j] for i in range(n)),F(0)) for j in range(n)]
 assert all(sum(K[i])==int(b[i]>0) for i in range(n))
 assert all(x<=C*w for x,w in zip(b,weights)) and all(x<=C*w for x,w in zip(transfer(b),weights))
 # Start with the upper obstacle b. Release violated capped coordinates.
 # Each exact M-matrix solve decreases r; released coordinates stay below b.
 capped=set(range(n));r=b.copy();policies=0
 while True:
  g=transfer(r);release={i for i in capped if theta*g[i]<b[i]}
  if not release:break
  capped-=release;U=sorted(set(range(n))-capped)
  mat=[[int(i==j)-theta*K[j][i] for j in U] for i in U]
  rhs=[theta*(a[i]+sum((b[j]*K[j][i] for j in capped),F(0))) for i in U]
  solution=linear(mat,rhs);previous=r.copy()
  r=b.copy()
  for i,v in zip(U,solution):r[i]=v
  assert all(0<=x<=y for x,y in zip(r,previous));policies+=1;assert policies<=n
 g=transfer(r);assert all(r[i]==min(b[i],theta*g[i]) for i in range(n))
 output={}
 for S,row in B.items():
  withdrawn=sum((v*r[i]/b[i] for i,v in row.items() if v),F(0))
  output[S]={i:row.get(i,F(0))*(1-r[i]/b[i] if b[i] else 1)+target_old[S].get(i,F(0))*withdrawn/old[S] for i in S}
  assert sum(output[S].values())==old[S]
 cols=[sum((row.get(i,F(0)) for row in output.values()),F(0))+a[i] for i in range(n)]
 assert cols==[b[i]-r[i]+g[i] for i in range(n)]
 assert all(0<=v<=(1+eps)*C*w for v,w in zip(cols,weights))
 withdrawn=sum(r);delta=sum(new.values(),F(0));movement=sum((abs(output[S].get(i,F(0))-row.get(i,F(0))) for S,row in B.items() for i in S),F(0))/2
 assert movement<=withdrawn<=theta/eps*delta
 result=dict(epsilon=str(eps),added=str(delta),withdrawn=str(withdrawn),actual_movement=str(movement),bound=str(theta/eps*delta),maximum_density=str(max((v/w for v,w in zip(cols,weights)),default=F(0))),C=str(C),policies=policies)
 return (result,r,output) if trace else result

def main():
 p=module('rounds/060/verify.py','prof62');alloc=module('rounds/056/allocation_probe.py','alloc62');prev=module('rounds/058/verify.py','prev62');flow=module('rounds/047/flow_probe.py','flow62');rec=module('rounds/046/threshold_probe.py','rec62')
 rows=[]
 with tempfile.TemporaryDirectory(prefix='euclidean62_') as tmp:
  oldruntime=rec.load(ROOT,Path(tmp)/'scratch.json')
  inputs=[]
  witness=json.loads((ROOT/'rounds/060/verification.json').read_text())['actual_frozen_counterexample']
  inputs.append((0,'round60_witness',[(F(z['x']),F(z['w'])) for z in witness['atoms']],F(1)))
  for batch,N in enumerate((3,5,8)):
   for serial in range(10):
    seed=620100+100*batch+serial;rng=random.Random(seed);xs=sorted(rng.sample(range(-512,513),N));ws=[rng.randint(1,64) for _ in xs]
    atoms=[(F(x,1024),F(3,8)*F(w,sum(ws))) for x,w in zip(xs,ws)]+[(F(10),F(5,8))]
    inputs.append((batch,f'seed{seed}',atoms,rng.choice((F(1,2),F(1),F(2),F(4)))))
  for batch,label,atoms,alpha in inputs:
   saved,r,groups,tasks=prev.capture_data(flow,rec,oldruntime,atoms,alpha);weights=[w for x,w in atoms]
   phi,T,_=p.profile(groups,weights,alloc);C=max(phi)
   old={S:d*F(1+i%3,4) for i,(S,d) in enumerate(sorted(groups.items(),key=lambda z:(len(z[0]),tuple(sorted(z[0])))))}
   new={S:groups[S]-v for S,v in old.items()};bp,B,_=p.profile(old,weights,alloc)
   for eps in (F(1),F(1,2),F(1,4),F(1,16)):
    previous_r=[F(0)]*len(weights);previous_output=B;path_movement=F(0)
    for fraction in (F(1,3),F(2,3),F(1)):
     incremental={S:fraction*v for S,v in new.items()}
     Tp={S:{i:v*(old[S]+incremental[S])/groups[S] for i,v in row.items()} for S,row in T.items()}
     result,rr,output=reroute(old,incremental,B,Tp,weights,C,eps,trace=True)
     assert all(x>=y for x,y in zip(rr,previous_r))
     step=sum((abs(output[S].get(i,F(0))-previous_output[S].get(i,F(0))) for S in old for i in S),F(0))/2
     assert step<=sum(rr)-sum(previous_r)
     path_movement+=step;assert path_movement<=F(result['withdrawn'])
     previous_r=rr;previous_output=output
     rows.append(dict(batch=batch,label=label,N=len(weights),D=saved['D'],prefix=str(fraction),path_movement=str(path_movement),**result))
 chain=[]
 for batch,Ls in enumerate(((2,3,4),(5,6,7),(8,9,10))):
  for L in Ls:
   N=2**L
   for eps in (F(1,2),F(1,8),F(1,32)):
    # Chain old left mass b_i=(N-i)/N, new singleton demand1.
    # Allowed column capacity1+eps. Prefix balance forces moved edge amount.
    moves=[max(F(0),F(N-i,N)-i*eps) for i in range(1,N)]
    left=[F(N-i,N)-z for i,z in enumerate(moves,1)];cols=[F(0)]*N;cols[0]=1
    for i,v in enumerate(left):assert 0<=v<=1;cols[i]+=v;cols[i+1]+=1-v
    assert all(v<=1+eps for v in cols)
    value=sum(moves);assert value<=1/(2*(eps+F(1,N)))
    chain.append(dict(batch=batch,N=N,epsilon=str(eps),optimal_movement_per_added=str(value),upper_bound=str((1-eps)/eps)))
   if N in (4,32,256):
    weights=[F(1,N)]*N;old={frozenset({i,i+1}):F(1,N) for i in range(N-1)};new={frozenset({0}):F(1,N)}
    B={S:{i:F(N-1-i,N*N),i+1:F(i+1,N*N)} for i,S in enumerate(old)}
    T={S:{i+1:F(1,N)} for i,S in enumerate(old)};T[frozenset({0})]={0:F(1,N)}
    for eps in (F(1,2),F(1,8)):
     result=reroute(old,new,B,T,weights,F(1),eps);rows.append(dict(batch=batch,label=f'chain{N}',N=N,D=None,**result))
 asymptotic=[]
 for batch,qs in enumerate(((2,3,4),(6,8,10),(12,16,24))):
  for q in qs:
   N=2**(2*q);eps=F(1,2**q);tau=eps+F(1,N)
   reciprocal=1/tau;m=(reciprocal.numerator+reciprocal.denominator-1)//reciprocal.denominator-1
   value=m-tau*m*(m+1)/2
   assert 0<eps*value<F(1,2) and F(1,2)-eps*value<=2*eps
   asymptotic.append(dict(batch=batch,q=q,N=N,epsilon=str(eps),epsilon_times_movement=str(eps*value),half_minus_value=str(F(1,2)-eps*value)))
 deps=json.loads((ROOT/'rounds/061/verification.json').read_text())['sha256'];deps['rounds/062/verify.py']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
 for path,h in deps.items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h
 print(json.dumps(dict(round=62,status='passed',python=sys.version.split()[0],rerouting=rows,chain=chain,asymptotic=asymptotic,sha256=deps,main_theorem_proved=False,scope='Exact finite certificates; general measure lemma has a separate analytic proof.'),ensure_ascii=False,indent=2))
if __name__=='__main__':main()
