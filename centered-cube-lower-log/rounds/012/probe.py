from pathlib import Path
import numpy as np,json
P=Path(__file__).resolve().parent

def labels(d,poly):
 K=(1<<d)-1;ex=[];z=1
 for _ in range(K):
  ex.append(z);z<<=1
  if z&(1<<d):z^=poly
 assert z==1 and len(set(ex))==K
 ex=np.array(ex);log=np.full(K+1,-1);log[ex]=np.arange(K)
 def mul(a,b):return ex[(log[a]+log[b])%K]
 for beta in range(2,K+1):
  W=[1,beta,1^beta];prod=[int(mul(W[i],W[j])) for i in range(3) for j in range(i,3)]
  D={(int(log[a])-int(log[b]))%K for a in W for b in W if a!=b}
  if len(set(prod))==6 and len(D)==6 and not D.intersection({1,2,K-1,K-2}) and not D.intersection({(z+1)%K for z in D}):break
 else:raise RuntimeError('no label')
 def trace(a):
  out=a.copy();cur=a.copy()
  for _ in range(d-1):cur=ex[(2*log[cur])%K];out^=cur
  assert np.all((out==0)|(out==1));return out
 U=(trace(ex)+2*trace(mul(ex,beta))+.5)/4
 y=(U*4-.5).astype(int);expected=np.full(16,(K+1)//16);expected[0]-=1
 audited=0
 for lag in range(1,K):
  if lag not in D:
   assert np.array_equal(np.bincount(y*4+np.roll(y,-lag),minlength=16),expected)
   audited+=1
 atoms=(np.arange(K)+U)/K;delta=1/(100*K)
 def primitive(t):
  y=np.mod(t,1);lo=np.searchsorted(atoms,y-delta/2);hi=np.searchsorted(atoms,y+delta/2)
  assert np.max(hi-lo)<=1
  partial=np.where(hi>lo,(y-atoms[np.minimum(lo,K-1)]+delta/2)/delta,0)
  return lo+partial-K*y
 return K,primitive,dict(d=d,audited_safe_integer_lags=audited,polynomial=poly,beta=beta,defects=sorted(D),label_counts=np.bincount((U*4-.5).astype(int),minlength=4).tolist())
K0,f0,c0=labels(10,0x409);K1,f1,c1=labels(12,0x1053)
m=np.arange(1,17);alpha=729/29;a=10/27;s=1/9
b=4*alpha*((np.sin(np.pi*m*a)/(np.pi*m))**4+3*(np.sin(np.pi*m*s)/(np.pi*m))**4)
rows=[]
for G in [32768,131072]:
 x=(np.arange(G)+.5)/G;cos=np.cos(2*np.pi*x[:,None]*m)
 def V(L):
  B=cos@(np.sqrt(2*b/3)*np.sin(np.pi*m*L))
  return (f0(x+L/2)-f0(x-L/2)+B)/L,(K0/K1)*(f1(x+L/2)-f1(x-L/2)+B)/L
 base=V(1.16)
 for w in [0.,.25,.5,1.]:
  v=(1-w)*base[0]+w*base[1];var=float(np.mean((v-v.mean())**2));items=[]
  for tau in [1/64,1/16,1/4,1]:
   other=V(1.16+tau/K0);y=(1-w)*other[0]+w*other[1]
   items.append(dict(tau=tau,increment=float(np.mean((v-y)**2)),correlation=float(np.corrcoef(v,y)[0,1])))
  rows.append(dict(G=G,fine_weight=w,scaled_variance=var,increments=items))
(P/'results.json').write_text(json.dumps(dict(labels=[c0,c1],K=[K0,K1],M=16,rows=rows),indent=2))
for r in rows:print(r['G'],r['fine_weight'],r['scaled_variance'],[(q['tau'],round(q['correlation'],6)) for q in r['increments']])
