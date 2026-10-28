from pathlib import Path
import numpy as np,json
P=Path(__file__).resolve().parent;records=[]
def rank(vs):
 basis={}
 for x in vs:
  v=int(x)
  while v:
   bit=v.bit_length()-1
   if bit in basis:v^=basis[bit]
   else:basis[bit]=v;break
 return len(basis)
for d,poly in [(10,0x409),(12,0x1053)]:
 K=(1<<d)-1;ex=[];z=1
 for _ in range(K):
  ex.append(z);z<<=1
  if z&(1<<d):z^=poly
 ex=np.array(ex);log=np.full(K+1,-1);log[ex]=np.arange(K)
 def sub(k):return [int(ex[k%K]),int(ex[(k+log[6])%K])]
 def trace(a):
  y=a.copy();cur=a.copy()
  for _ in range(d-1):cur=ex[2*log[cur]%K];y^=cur
  return y
 labels=trace(ex)+2*trace(ex[(log[ex]+log[6])%K]);second=12
 assert rank(sub(0)+sub(second))==4
 total=np.zeros(64,dtype=np.int64);eligible=0;defective=0
 for c in range(K):
  if rank(sub(0)+sub(c))==4 and rank(sub(second)+sub(c))==4:
   eligible+=1;defective+=rank(sub(0)+sub(second)+sub(c))<6
   total+=np.bincount(labels*16+np.roll(labels,-second)*4+np.roll(labels,-c),minlength=64)
 ratio=total*64/(eligible*K)
 bound=(K+1)/K*(1+3*45/(K-18))
 assert ratio.max()<=bound
 records.append(dict(d=d,K=K,second=second,eligible=eligible,defective=defective,min_ratio=float(ratio.min()),max_ratio=float(ratio.max()),theoretical_upper=bound,total_counts=total.tolist()))
(P/'average_check.json').write_text(json.dumps(records,indent=2));print(json.dumps([{k:v for k,v in r.items() if k!='total_counts'} for r in records],indent=2))
