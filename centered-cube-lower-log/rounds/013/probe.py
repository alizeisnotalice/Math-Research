from pathlib import Path
import numpy as np,json
P=Path(__file__).resolve().parent;K=1023;poly=0x409;ex=[];z=1
for _ in range(K):
 ex.append(z);z<<=1
 if z&1024:z^=poly
ex=np.array(ex);log=np.full(1024,-1);log[ex]=np.arange(K)
def trace(a):
 out=a.copy();cur=a.copy()
 for _ in range(9):cur=ex[(2*log[cur])%K];out^=cur
 return out
U=(trace(ex)+2*trace(ex[(log[ex]+log[6])%K])+.5)/4
rng=np.random.default_rng(2026102808)
configs={'field':U,'permuted':rng.permutation(U),'sorted':np.sort(U),'reversed':U[::-1]}
rows=[]
for name,u in configs.items():
 atoms=(np.arange(K)+u)/K;gap=np.min(np.diff(np.r_[atoms,atoms[0]+1]))
 for tau in [1/64,1/16,1/4]:
  delta=tau/K;t=delta/2;assert t<gap
  for typ,L in [('generic',1.16),('aligned',(K+164+.25-tau/2)/K)]:
   c1=np.mod(atoms-L/2-t/2,1);c2=np.sort(np.mod(atoms+L/2+t/2,1));ext=np.r_[c2-1,c2,c2+1]
   idx=np.searchsorted(ext,c1);overlap=0.
   for off in [-1,0]:overlap+=np.maximum(0,t-np.abs(c1-ext[idx+off])).sum()
   exact_float=K*delta+2*overlap-(K*delta)**2
   rows.append(dict(name=name,tau=tau,kind=typ,L=L,overlap=float(overlap),increment=float(exact_float),ratio=float(exact_float/(K*delta))))
(P/'results.json').write_text(json.dumps(rows,indent=2))
for r in rows:print(r['name'],r['tau'],r['kind'],r['ratio'])
