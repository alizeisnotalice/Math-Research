import json,random,itertools,time
from fractions import Fraction as F
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
HERE=Path(__file__).resolve().parent

def average(x,R,w,ends):
 p=[F(1)];den=F(1)
 for t,c in zip(x,ends):
  inner=min(t,F(1))+min(c,F(1));support=min(t,R)+min(c,R);outer=support-inner
  z=[F(0)]*(len(p)+1)
  for j,a in enumerate(p):z[j]+=a*outer;z[j+1]+=a*inner
  p=z;den*=t+c
 return sum(a*b for a,b in zip(p,w))/den

def one(task):
 n,R=task;rng=random.Random(739000+100*n+R);tests=0
 for sample in range(30):
  x=[rng.choice([F(1,4),F(3,4),F(1),F(R+1,2),F(R),F(R+1),F(2*R),F(3*R)]) for i in range(n)]
  w=sorted(F(rng.randrange(0,11),5) for i in range(n+1));outs=[i for i,t in enumerate(x) if t>R]
  brute=max(average(x,F(R),w,c) for c in itertools.product([F(1),F(R)],repeat=n))
  order=sorted(outs,key=lambda i:x[i],reverse=True);vals=[]
  for r in range(len(outs)+1):
   selected=set(order[:r]);c=[F(R) if i in selected else F(1) for i in range(n)];vals.append(average(x,F(R),w,c))
  assert max(vals)==brute,(n,R,x,w,vals,brute)
  tests+=1
 return dict(n=n,R=R,exact_tests=tests,pass_all=True)
if __name__=='__main__':
 start=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:res=list(pool.map(one,[(n,R) for n in range(1,9) for R in [2,4,8]]))
 out=dict(status='complete',exact_fraction_arithmetic=True,seed_rule='739000+100*n+R',cases=sum(r['exact_tests'] for r in res),workers=12,scope='endpoint subset exhaustive equality versus sorted candidates, monotone soft count profiles',results=res,elapsed=time.monotonic()-start)
 (HERE/'check_endpoints.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['status','cases','elapsed']}))
