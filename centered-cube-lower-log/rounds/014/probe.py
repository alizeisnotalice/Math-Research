from pathlib import Path
import numpy as np,json
P=Path(__file__).resolve().parent;K=1023;ex=[];z=1
for _ in range(K):
 ex.append(z);z<<=1
 if z&1024:z^=0x409
ex=np.array(ex);log=np.full(1024,-1);log[ex]=np.arange(K);beta=6
D={67,78,145,878,945,956}
def rank(vs):
 basis={}
 for value in vs:
  v=int(value)
  while v:
   bit=v.bit_length()-1
   if bit in basis:v^=basis[bit]
   else:basis[bit]=v;break
 return len(basis)
def rows(a):return [int(ex[a%K]),int(ex[(a+log[beta])%K])]
def safe(a,b):return a!=b and ((a-b)%K) not in D
rng=np.random.default_rng(2026102809);counts={};accepted=0
for a,b in rng.integers(1,K,size=(20000,2)):
 if safe(0,int(a)) and safe(0,int(b)) and safe(int(a),int(b)):
  r=rank(rows(0)+rows(int(a))+rows(int(b)));counts[r]=counts.get(r,0)+1;accepted+=1
candidates=[]
for aa in range(8,K):
 bb=int(log[1^int(ex[aa])])
 if bb>=8 and abs(aa-bb)>=8 and safe(0,aa) and safe(0,bb) and safe(aa,bb):candidates.append((max(aa,bb),aa,bb))
_,a,b=min(candidates)
witness=dict(a=a,b=b,ranks=[rank(rows(0)+rows(a)),rank(rows(0)+rows(b)),rank(rows(a)+rows(b))],triple_rank=rank(rows(0)+rows(a)+rows(b)))
def trace(v):
 cur=v.copy();res=v.copy()
 for _ in range(9):cur=ex[2*log[cur]%K];res^=cur
 return res
y=trace(ex)+2*trace(ex[(log[ex]+log[6])%K]);joint=np.bincount(y*16+np.roll(y,-a)*4+np.roll(y,-b),minlength=64)
witness.update(nonzero_cells=int(np.count_nonzero(joint)),counts=joint.tolist())
# Complete H+B Gram for fixed spatial sample grid, three label arrangements.
m=np.arange(1,17);bm=4*(729/29)*((np.sin(np.pi*m*10/27)/(np.pi*m))**4+3*(np.sin(np.pi*m/9)/(np.pi*m))**4)
L=np.linspace(10/9,11/9,65);summary=[]
for G in [8192,32768]:
 x=(np.arange(G)+.5)/G;B=np.cos(2*np.pi*x[:,None]*m)@(np.sqrt(2*bm/3)[:,None]*np.sin(np.pi*m[:,None]*L))
 for name,yy in [('field',y),('random',np.random.default_rng(2026102809).permutation(y)),('sorted',np.sort(y))]:
  U=(yy+.5)/4
  def prim(t):
   kx=K*np.mod(t,1);j=np.floor(kx).astype(int);r=kx-j;return (r>=U[j])-r
  H=prim(x[:,None]+L/2)-prim(x[:,None]-L/2);A=(H+B)/L;A-=A.mean(axis=0);gram=A.T@A/G;var=np.diag(gram);rho=gram/np.sqrt(var[:,None]*var[None,:]);mask=(L[:,None]-L[None,:])>4/K;rv=rho[mask]
  summary.append(dict(G=G,name=name,var_min=float(var.min()),var_max=float(var.max()),far_rho_median=float(np.median(rv)),far_rho95=float(np.quantile(rv,.95)),far_rho_max=float(rv.max())))
out=dict(pair_safe_samples=accepted,rank_counts=counts,witness=witness,gram=summary)
(P/'results.json').write_text(json.dumps(out,indent=2));print('triple',accepted,counts,witness);print('gram',summary)
