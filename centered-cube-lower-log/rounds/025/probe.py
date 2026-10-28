import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['VECLIB_MAXIMUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import numpy as np,json,time,multiprocessing as mp,argparse
P=Path(__file__).resolve().parent
source=(P.parent/'017/probe.py').read_text();exec(source.split('rows=[]')[0],globals());P=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=4);parser.add_argument('--per-worker',type=int,default=512);parser.add_argument('--grid-size',type=int,choices=[262144,1048576],default=262144);parser.add_argument('--seed',type=int,default=2026102909);parser.add_argument('--output',default='results.json');args=parser.parse_args()
start=time.perf_counter();G=args.grid_size;n=900000000;a=.5/K;h=30.;eta=2.;x=(np.arange(G)+.5)/G;pp=prim(x);fp=np.fft.fft(pp);auto=np.fft.ifft(abs(fp)**2).real/G
kk=np.arange(int(np.ceil((10/9)*G/2)),int(np.floor((11/9)*G/2))+1);L=2*kk/G;r=len(kk);c=np.sqrt(2*b/3);sins=np.sin(np.pi*L[:,None]*m);cosx=np.cos(2*np.pi*x[:,None]*m);psin=pp@np.sin(2*np.pi*x[:,None]*m)/G
variance=(2*(np.mean(pp*pp)-auto[(2*kk)%G])+4*(sins*sins)@(c*psin)+.5*(sins*sins)@(c*c))/(L*L)
sizes=[65,129,257,513,1025,2049,4097,8193,r];ids=[np.unique(np.rint(np.linspace(0,r-1,size)).astype(int)) for size in sizes]
assert all(np.isin(ids[j],ids[j+1]).all() for j in range(len(ids)-1))
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
  z=values[:,:,j];means=z.mean(axis=0);se=z.std(axis=0,ddof=1)/np.sqrt(len(z));rows.append(dict(scales=len(ids[j]),N=len(z),A_plus=float(means[0]),A_center=float(means[1]),A_minus=float(means[2]),se_plus=float(se[0]),se_center=float(se[1]),se_minus=float(se[2]),lower_transfer_score=float(np.exp(-eta)*means[0]),upper_transfer_score=float(np.exp(eta)*means[2])))
 np.save(P/(args.output+'.raw.npy'),values)
 gain=values[:,1,-1]-values[:,1,-2]
 paired=dict(mean=float(gain.mean()),se=float(gain.std(ddof=1)/np.sqrt(len(gain))))
 result=dict(paired_last_gain=paired,seconds=time.perf_counter()-start,seed=args.seed,workers=args.workers,per_worker=args.per_worker,G=G,r=r,n=n,h=h,eta=eta,n_a2_q=float(n*a*a*2),log_P_bad=float(np.log(2*r)-1100),log_Q_bad_bound=float((np.log(2*r)-1100+n*a*a*2)/2),rows=rows)
 (P/args.output).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
