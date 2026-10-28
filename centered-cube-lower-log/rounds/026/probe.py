import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['VECLIB_MAXIMUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import numpy as np,json,time,multiprocessing as mp,argparse
from fractions import Fraction as Q
from math import isqrt
P=Path(__file__).resolve().parent
source=(P.parent/'017/probe.py').read_text();exec(source.split('rows=[]')[0],globals());P=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=4);parser.add_argument('--per-worker',type=int,default=512);parser.add_argument('--dimension',type=int,choices=[10,12,14],default=10);parser.add_argument('--grid-size',type=int,choices=[262144,1048576,4194304],default=262144);parser.add_argument('--seed',type=int,default=2026102909);parser.add_argument('--output',default='results.json');args=parser.parse_args()
# Build and verify each finite-field source before any Monte Carlo output.
d=args.dimension;polynomial={10:0x409,12:0x1053,14:0x4443}[d];K=(1<<d)-1
z=1;exlist=[]
for _ in range(K):
 exlist.append(z);z<<=1
 if z&(1<<d):z^=polynomial
assert z==1 and len(set(exlist))==K
ex=np.array(exlist);log=np.full(K+1,-1);log[ex]=np.arange(K)
for beta in range(2,K+1):
 W=[1,beta,1^beta];products=[int(ex[(log[W[i]]+log[W[j]])%K]) for i in range(3) for j in range(i,3)]
 defects={(int(log[u])-int(log[v]))%K for u in W for v in W if u!=v}
 if len(set(products))==6 and len(defects)==6 and not defects.intersection({1,2,K-1,K-2}) and not defects.intersection({(v+1)%K for v in defects}):break
else:raise RuntimeError('no valid labels')
def trace_field(v):
 cur=v.copy();out=v.copy()
 for _ in range(d-1):cur=ex[(2*log[cur])%K];out^=cur
 assert np.all((out==0)|(out==1));return out
labels=trace_field(ex)+2*trace_field(ex[(log[ex]+log[beta])%K]);label_counts=np.bincount(labels,minlength=4).tolist();assert label_counts==[(K+1)//4-1]+[(K+1)//4]*3
U=(labels+.5)/4;atoms=(np.arange(K)+U)/K;delta=1/(100*K)
field=dict(dimension=d,K=K,polynomial=polynomial,beta=beta,defects=sorted(defects),label_counts=label_counts)

start=time.perf_counter();G=args.grid_size;n=(900000000*K*K+1023**2//2)//1023**2;a=.5/K;h=30.;eta=2.;x=(np.arange(G)+.5)/G;pp=prim(x);fp=np.fft.fft(pp);auto=np.fft.ifft(abs(fp)**2).real/G
kk=np.arange(int(np.ceil((10/9)*G/2)),int(np.floor((11/9)*G/2))+1);L=2*kk/G;r=len(kk);c=np.sqrt(2*b/3);sins=np.sin(np.pi*L[:,None]*m);cosx=np.cos(2*np.pi*x[:,None]*m);psin=pp@np.sin(2*np.pi*x[:,None]*m)/G
variance=(2*(np.mean(pp*pp)-auto[(2*kk)%G])+4*(sins*sins)@(c*psin)+.5*(sins*sins)@(c*c))/(L*L)
sizes=[s for s in [65,129,257,513,1025,2049,4097,8193,16385,32769,65537,131073] if s<r]+[r];ids=[np.unique(np.rint(np.linspace(0,r-1,size)).astype(int)) for size in sizes]
assert all(np.isin(ids[j],ids[j+1]).all() for j in range(len(ids)-1))
EP2=sum(int(count)*(Q(2*j+1,8)**2-Q(2*j+1,8)+Q(1,3)) for j,count in enumerate(label_counts))/K
Sbound=sum((4*Q(729,29)*(min(Q(10,27),Q(1,3*j))**4+3*min(Q(1,9),Q(1,3*j))**4) for j in range(1,17)),Q(0))
q_bound=(8*(EP2+Q(2*K,G))+2*Sbound/3)*Q(81,100);assert q_bound<2
aa=Q(1,2*K);AA=Q(77,15);rad=2*n*AA*AA*2*1100;root=Q(isqrt((rad.numerator*10**24)//rad.denominator)+1,10**12);dev=root+2*AA*AA*1100/3
eta_bound=aa*aa*dev/2+aa**3*AA*(n*2+dev)/(3*(1-aa*AA));assert eta_bound<2 and n*aa*aa*2<431
qualification=dict(field=field,n=n,G=G,r=r,na2=float(n*aa*aa),q_upper_rational=str(q_bound),eta_upper_rational=str(eta_bound),eta_upper_display=float(eta_bound),variance_range=[float(variance.min()),float(variance.max())])
if __name__=='__main__':(P/(args.output+'.qualification.json')).write_text(json.dumps(qualification,indent=2))
def worker(index):
 rng=np.random.default_rng(np.random.SeedSequence([args.seed,G,index]));records=[]
 for rep in range(args.per_worker):
  winner=int(rng.integers(r));k=int(kk[winner]);ll=L[winner];X=(np.roll(pp,-k)-np.roll(pp,k)+cosx@(c*np.sin(np.pi*m*ll)))/ll
  q=(1+a*X)/G;q/=q.sum();counts=rng.multinomial(n,q)
  corr=np.fft.ifft(np.conj(np.fft.fft(counts))*fp).real;linear=(corr[kk%G]-corr[(-kk)%G]+sins@(c*(counts@cosx)))/L
  T=a*linear-n*a*a*variance/2;mx=T.max();logD=mx+np.log(np.exp(T-mx).mean());w=float(np.exp(h-logD));maxima=np.array([T[jj].max() for jj in ids]);records.append(np.stack([w*(maxima>h+eta),w*(maxima>h),w*(maxima>h-eta)]).tolist())
 return records
if __name__=='__main__':
 print('grid',G,'scales',r,'workers',args.workers,flush=True)
 with mp.get_context('fork').Pool(args.workers) as pool:parts=pool.map(worker,range(args.workers))
 values=np.array([rec for part in parts for rec in part]);rows=[]
 for j,size in enumerate(sizes):
  z=values[:,:,j];means=z.mean(axis=0);se=z.std(axis=0,ddof=1)/np.sqrt(len(z));rows.append(dict(normalized_center=float(means[1]/np.sqrt(n)),normalized_se=float(se[1]/np.sqrt(n)),scales=len(ids[j]),N=len(z),A_plus=float(means[0]),A_center=float(means[1]),A_minus=float(means[2]),se_plus=float(se[0]),se_center=float(se[1]),se_minus=float(se[2]),lower_transfer_score=float(np.exp(-eta)*means[0]),upper_transfer_score=float(np.exp(eta)*means[2])))
 np.save(P/(args.output+'.raw.npy'),values)
 result=dict(qualification=qualification,seconds=time.perf_counter()-start,seed=args.seed,workers=args.workers,per_worker=args.per_worker,G=G,r=r,n=n,h=h,eta=eta,n_a2_q=float(n*a*a*2),log_P_bad=float(np.log(2*r)-1100),log_Q_bad_bound=float((np.log(2*r)-1100+n*a*a*2)/2),rows=rows)
 (P/args.output).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
