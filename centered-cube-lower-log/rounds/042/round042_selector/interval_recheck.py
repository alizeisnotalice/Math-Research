"""Outward binary64 enclosure of the intended rational source at frozen dyadic queries.
All arithmetic in kernels, DP, positive sums/products uses adjacent-float expansion.
Does not certify sampling or a continuum query integral.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,math,sys,time
import numpy as np
import probe as core
HERE=Path(__file__).resolve().parent
DOWN=-np.inf;UP=np.inf

def dn(a):return np.nextafter(a,DOWN)
def up(a):return np.nextafter(a,UP)
def rational(a):
 a=F(a);v=float(a);f=F.from_float(v)
 return (float(np.nextafter(v,DOWN)) if f>a else v,
         float(np.nextafter(v,UP)) if f<a else v)
def add(a,b):return (np.maximum(0,dn(a[0]+b[0])),up(a[1]+b[1]))
def mul(a,b):return (np.maximum(0,dn(a[0]*b[0])),up(a[1]*b[1]))
def div(a,b):
 assert np.all(np.asarray(b[0])>0)
 return (np.maximum(0,dn(a[0]/b[1])),up(a[1]/b[0]))
def isum(a):
 lo=np.zeros(a[0].shape[:-1]);hi=lo.copy()
 for j in range(a[0].shape[-1]):lo,hi=add((lo,hi),(a[0][...,j],a[1][...,j]))
 return lo,hi

def overlap(x,L,left,right):
 h=rational(L/2);wl=(dn(x-h[1]),up(x-h[0]));wr=(dn(x+h[0]),up(x+h[1]))
 lo=np.zeros_like(x);hi=lo.copy()
 for shift in range(-2,3):
  al,au=rational(left+shift);bl,bu=rational(right+shift)
  ol=np.maximum(0.,dn(np.minimum(wr[0],bl)-np.maximum(wl[1],au)))
  oh=np.maximum(0.,up(np.minimum(wr[1],bu)-np.maximum(wl[0],al)))
  lo,hi=add((lo,hi),(ol,oh))
 return lo,hi

def kernels(x,c):
 N,n=x.shape;K=c['K'];R=c['R'];lengths=[F(str(v))for v in c['lengths']]
 assert (K,R,F(str(c['epsilon'])),F(str(c['delta'])))==(4,4,F(2,25),F(1,20))
 lo=np.empty((R,n,N,len(lengths),K));hi=np.empty_like(lo)
 for b in range(R):
  for j in range(K):
   phase=F(2*b+1,2*R) if j%2==0 else F(2*R-2*b-1,2*R)
   center=(j+phase)/K;delta=F(1,20)
   assert F(j,K)<=center-delta/2<center+delta/2<=F(j+1,K)
   for li,L in enumerate(lengths):
    a=mul(overlap(x,L,F(j,K),F(j+1,K)),rational(F(2,25)))
    z=mul(overlap(x,L,center-delta/2,center+delta/2),rational(F(23,5)))
    kl,kh=div(add(a,z),rational(L));lo[b,:,:,li,j]=kl.T;hi[b,:,:,li,j]=kh.T
 assert np.isfinite(hi).all()and np.all(lo>=0)and np.all(lo<=hi)
 return lo,hi

def prior_mass(types,s,c):
 pi=[F(0)for _ in range(c['R'])]
 for a,b in zip(types,s):pi[int(b)]+=F(math.factorial(c['n']),math.prod(math.factorial(v)for v in a)*c['K']**c['n'])
 assert sum(pi)==1
 return pi

def bounds(x,c,s):
 zlo,zhi=kernels(x,c);layers,tr=core.transitions(c['n'],c['K']);N=len(x);J=len(c['lengths'])
 gate=(np.zeros((N,J)),np.zeros((N,J)));fixed=[]
 for b in range(c['R']):
  lo=np.ones((N,J,1));hi=lo.copy();f=(np.ones((N,J)),np.ones((N,J)))
  for i,tab in enumerate(tr):
   nl=np.zeros((N,J,len(layers[i+1])));nh=nl.copy()
   for j,(dst,src) in enumerate(tab):
    q=mul((lo[:,:,src],hi[:,:,src]),(zlo[b,i,:,:,j,None],zhi[b,i,:,:,j,None]))
    nl[:,:,dst],nh[:,:,dst]=add((nl[:,:,dst],nh[:,:,dst]),q)
   lo,hi=nl,nh
   f=mul(f,isum((zlo[b,i],zhi[b,i])))
  selected=isum((lo[...,s==b],hi[...,s==b]));gate=add(gate,selected);fixed.append(f)
 pi=prior_mass(layers[-1],s,c);prior=(np.zeros((N,J)),np.zeros((N,J)))
 for b,f in enumerate(fixed):prior=add(prior,mul(f,rational(pi[b])))
 lows=np.array([gate[0],prior[0]]+[f[0]for f in fixed]);highs=np.array([gate[1],prior[1]]+[f[1]for f in fixed])
 # rho_b(t+1/2)=rho_b(t); length3/2 contains exactly3 periods.
 # Every fixed-table product average and their actual-prior mixture equal1.
 for j,L in enumerate(c['lengths']):
  if F(str(L))==F(3,2):lows[1:,:,j]=1.;highs[1:,:,j]=1.
 assert np.all(lows<=highs)and np.isfinite(highs).all()
 return lows,highs

def run_case(folder,batch=256):
 plan=json.loads((folder/'PLAN.json').read_text());c=plan.get('config');frozen=json.loads((folder/'FROZEN_SELECTOR.json').read_text())
 assert core.sha(folder/'PLAN.json')==frozen['plan_SHA256']
 s=np.array(frozen['selector']);lam=F(str(frozen['lambda_']));assert lam==1
 alpha=.05/36 # global36 directions: three frozen cases × six models × two.
 rng=np.random.default_rng(c['test_seed']);lower=np.zeros(6,dtype=int);upper=lower.copy();maxwidth=0.;max_by_model=np.zeros(6)
 for start in range(0,c['test_N'],batch):
  N=min(batch,c['test_N']-start);x=rng.random((N,c['n']));lo,hi=bounds(x,c,s)
  lower+=(lo.max(axis=2)>float(lam)).sum(axis=1);upper+=(hi.max(axis=2)>float(lam)).sum(axis=1)
  width=(hi-lo).max(axis=(1,2));max_by_model=np.maximum(max_by_model,width)
  maxwidth=max(maxwidth,float(width.max()))
 rows=[]
 for name,l,u in zip(['gated','prior','fixed0','fixed1','fixed2','fixed3'],lower,upper):
  lb=core.ci(int(l),c['test_N'],alpha)[0];ub=core.ci(int(u),c['test_N'],alpha)[1]
  rows.append(dict(model=name,lower_successes=int(l),upper_successes=int(u),ambiguous=int(u-l),N=c['test_N'],weak_CI=[float(lam)*lb,float(lam)*ub]))
 return dict(n=c['n'],rows=rows,plan_SHA256=core.sha(folder/'PLAN.json'),selector_SHA256=core.sha(folder/'FROZEN_SELECTOR.json'),
   max_per_scale_interval_width=maxwidth,per_model_max_width=max_by_model.tolist(),
   simultaneous_gated_lower_minus_prior_upper=rows[0]['weak_CI'][0]-rows[1]['weak_CI'][1],
   simultaneous_gated_lower_minus_best_fixed_upper=rows[0]['weak_CI'][0]-max(row['weak_CI'][1]for row in rows[2:]))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',default=str(HERE/'INTERVAL_RESULTS.json'));ap.add_argument('--batch',type=int,default=256);args=ap.parse_args()
 outpath=Path(args.output);assert not outpath.exists();start=time.time()
 planfile=outpath.with_name(outpath.stem+'_PLAN.json');assert not planfile.exists()
 plan=dict(script_SHA256=core.sha(__file__),core_SHA256=core.sha(core.__file__),alpha_each=.05/36,
   exact_source={'epsilon':'2/25','delta':'1/20','K':4,'R':4},
   cases=[dict(folder=str(f.relative_to(HERE)),plan_SHA256=core.sha(f/'PLAN.json'),selector_SHA256=core.sha(f/'FROZEN_SELECTOR.json'))for f in [HERE,HERE/'balanced/n7',HERE/'balanced/n15']],
   scope='Repair original frozen observations only; no tuning/search/new samples; all positive arithmetic outward, exact period identity; inference remains Monte Carlo, not continuum deterministic certificate')
 planfile.write_text(json.dumps(plan,indent=2)+'\n');cases=[]
 for folder in [HERE,HERE/'balanced/n7',HERE/'balanced/n15']:
  cases.append(run_case(folder,args.batch));print(json.dumps(cases[-1]),flush=True)
 out=dict(cases=cases,plan_SHA256=core.sha(planfile),alpha_each=.05/36,
  source_model='Exact rational source, dyadic query coordinates as generated',
  local_arithmetic_enclosure='Binary64 correctly rounded basic operations expanded one adjacent float per operation; no libm kernels, reductions, BLAS, or chosen tolerance',
  exact_period_identity=True,original_strict_boundary_counts_replaced=True,
  deterministic_continuum_integral_certificate=False,Cn_log_gain_proved=False,
  python=sys.version,numpy=np.__version__,seconds=time.time()-start)
 outpath.write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
