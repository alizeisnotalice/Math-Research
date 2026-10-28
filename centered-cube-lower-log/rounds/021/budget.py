from pathlib import Path
from fractions import Fraction as F
from bisect import bisect_left,bisect_right
from collections import Counter
import numpy as np,json
P=Path(__file__).resolve().parent;K=1023;delta=F(1,100*K);a=F(1,2*K);n=900000000
z=1;ex=[]
for _ in range(K):
 ex.append(z);z<<=1
 if z&1024:z^=0x409
ex=np.array(ex);lg=np.full(1024,-1);lg[ex]=np.arange(K)
def tr(v):
 c=v.copy();out=v.copy()
 for _ in range(9):c=ex[2*lg[c]%K];out^=c
 return out
y=tr(ex)+2*tr(ex[(lg[ex]+lg[6])%K]);atoms=[F(8*j+2*int(v)+1,8*K) for j,v in enumerate(y)]
rows=[]
def up(q,digits=60):
 unit=10**digits;num=(q.numerator*unit+q.denominator-1)//q.denominator
 return str(num)+'/'+str(unit)
for name,L in [('near',[F(29,25)+F(j,16*K) for j in [0,1,4]]),('far',[F(113,100),F(117,100),F(121,100)])]:
 edges=sorted({F(0),F(1)}|{(atom+edge+shift)%1 for ll in L for atom in atoms for edge in [-delta/2,delta/2] for shift in [-ll/2,ll/2]});groups=Counter();qual=F(0);maxq=F(0)
 for lo,hi in zip(edges[:-1],edges[1:]):
  mid=(lo+hi)/2;slopes=[]
  for ll in L:
   counts=[]
   for t in [(mid+ll/2)%1,(mid-ll/2)%1]:counts.append(bisect_right(atoms,t+delta/2)-bisect_left(atoms,t-delta/2))
   slopes.append(abs(counts[0]-counts[1])/delta)
  slope=max(slopes);r=F(1,64) if slope==0 else F(1,256)
  bound=a*F(9,10)*(F(154,27)+r*(slope+F(3520,7)));qual=max(qual,bound);assert bound<F(1,2)
  q=(hi-lo)/(2*r);maxq=max(maxq,q);assert q<1
  groups[(hi-lo,r)]+=1
 for order in [2,4,8,12]:
  E=sum((count*F(9,2)*length*((length/(2*r))**(2*order))/(1-length/(2*r)) for (length,r),count in groups.items()),F(0))
  rows.append(dict(name=name,order=order,segments=len(edges)-1,width_classes=len(groups),max_disk_bound=str(qual),max_radius_ratio=str(maxq),integral_error_upper=up(E),n_log_error_upper=up(4*n*E),display_n_log_bound=float(4*n*E)))
(P/'results.json').write_text(json.dumps(rows,indent=2));print(json.dumps([{k:v for k,v in r.items() if k in ['name','order','display_n_log_bound','segments','width_classes','max_radius_ratio']} for r in rows],indent=2))
