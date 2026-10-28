from pathlib import Path
from fractions import Fraction as F
from bisect import bisect_left,bisect_right
import numpy as np,json,math
P=Path(__file__).resolve().parent
source=P.parent/'021/budget.py';exec(source.read_text().split('rows=[]')[0],globals());P=Path(__file__).resolve().parent
results=[]
for name,L in [('near',[F(29,25)+F(j,16*K) for j in [0,1,4]]),('far',[F(113,100),F(117,100),F(121,100)])]:
 edges=sorted({F(0),F(1)}|{(atom+edge+shift)%1 for ll in L for atom in atoms for edge in [-delta/2,delta/2] for shift in [-ll/2,ll/2]});V=F(0);eps=F(0);bias=F(0)
 for lo,hi in zip(edges[:-1],edges[1:]):
  mid=(lo+hi)/2;slopes=[]
  for ll in L:
   counts=[bisect_right(atoms,t+delta/2)-bisect_left(atoms,t-delta/2) for t in [(mid+ll/2)%1,(mid-ll/2)%1]]
   slopes.append(abs(counts[0]-counts[1])/delta)
  s=max(slopes);length=hi-lo;lip=2*a*F(9,10)*(s+F(3520,63));e=length*lip;eps=max(eps,e);V+=length*e*e
  r=F(1,64) if s==0 else F(1,256);q=length/(2*r);bias+=2*length*q**16/(1-q)
 ell=96.0;shift=float(n*bias)+math.sqrt(2*n*float(V)*ell)+4*float(eps)*ell/3
 row=dict(log_union_budget=96,name=name,eps=str(eps),second_moment_bound=str(V),n_mean_error_bound=float(n*bias),deterministic_shift=float(n*eps),n_variance_bound=float(n*V),failure_probability=1e-40,bernstein_shift=shift)
 row["subdivision_shifts"]={str(S):float(n*bias)/S**16+(shift-float(n*bias))/S for S in [1,8,32,64,128,256]}
 row['certified_integer_shifts']={}
 for S,display in row['subdivision_shifts'].items():
  SS=int(S);d=math.ceil(display);remainder=F(d)-n*bias/SS**16-128*eps/SS
  while remainder<=0 or remainder**2<192*n*V/SS**2:
   d+=1;remainder=F(d)-n*bias/SS**16-128*eps/SS
  row['certified_integer_shifts'][S]=d
 results.append(row);print(json.dumps(row),flush=True)
(P/'results.json').write_text(json.dumps(results,indent=2))
