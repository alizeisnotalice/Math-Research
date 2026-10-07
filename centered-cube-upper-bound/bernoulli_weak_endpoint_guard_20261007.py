"""Original kernel weak-level pressure, registered inputs; numerical screens are not intervals."""
from pathlib import Path
from fractions import Fraction as F
from math import comb, isqrt, ceil, log, pi
from hashlib import sha256
import json, ctypes, subprocess, time
import numpy as np

HERE=Path(__file__).resolve().parent
REG=HERE/'bernoulli_weak_endpoint_registration_20261007.json'
OUT=HERE/'bernoulli_weak_endpoint_guard_results_20261007.json'
CPP=HERE/'bernoulli_weak_endpoint_dp_20261007.cpp'
LIB=HERE/'bernoulli_weak_endpoint_dp_20261007.dylib'
assert REG.exists() and not OUT.exists(), 'Preserve all prior results'
started=time.time()
subprocess.run(['clang++','-O3','-std=c++17','-shared','-fPIC',str(CPP),'-o',str(LIB)],check=True)
lib=ctypes.CDLL(str(LIB))
ptr=ctypes.POINTER(ctypes.c_double)
lib.profiles.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.c_int,ptr,ptr,ptr,
                      ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double,ptr]
def pointer(a): return a.ctypes.data_as(ptr)
K=16384; table_n=262144; C=512/511
tables={}; errors={}
for c in (1.0,0.5):
 for period in (1.0,2.0,4.0):
  spec=np.zeros(table_n//2+1,dtype=np.complex128); spec[0]=0.5*table_n
  ks=np.arange(1,K+1,2); v=(2*pi*ks/period)**2
  gc=1/(1+c*(v/np.log1p(v)-1))
  spec[ks]=-1j*table_n*gc/(pi*ks)
  table=np.fft.irfft(spec,n=table_n)
  aa=(2*pi/period)**2
  tail=4/(pi*c*aa)*(log(K)+0.5+0.5*log(1+aa))/(K*K)
  # G_c<=2w for c>=1/2, periodic peak <=2*coth(period/2).
  lip=2*(1+2/np.expm1(period))
  interpolation=lip*period/(2*table_n)
  tables[c,period]=table
  errors[c,period]={'fourier_tail':tail,'interpolation':interpolation,
                    'total_analytic_error':tail+interpolation,
                    'floating_screen_allowance':1e-10,
                    'not_interval':True,
                    'range_unclipped':[float(table.min()),float(table.max())]}

rounds=[]; exact_checks=0; diagnostics=0
for n,M,ns in ((8,128,1024),(32,256,512),(128,512,512)):
 rrq=[F(j*j,j*j+(M-j)*(M-j)) for j in range(M+1)]
 rr=np.asarray([float(q) for q in rrq],dtype=np.float64)
 rs=F(0)
 for j in range(M):
  if j==0 or j==M-1: rj=F(1,(M-1)**2)
  else:
   l,u=rrq[j],rrq[j+1]
   sl=F(j*(M-j),j*j+(M-j)*(M-j))
   sj=j+1; su=F(sj*(M-sj),sj*sj+(M-sj)*(M-sj))
   rj=((u-l)/(sl+su))**2
  rs=max(rs,rj)
 assert rs==F(4,M*M) and n*rs==F(1,512)
 exact_checks+=M+2
 shift=isqrt(n)//2; mid=n//2
 central_q=F(comb(n,mid),2**n)
 shifted_q=F(comb(n,mid-shift)+comb(n,mid+shift),2**n)
 band_q=F(sum(comb(n,j) for j in (mid-1,mid,mid+1)),2**n)
 central=float(central_q); shifted=float(shifted_q); band=float(band_q)
 # Equal mixture of separately normalized symmetric shells = their normalized union.
 heights=np.asarray([1/central,1/shifted,1/band,
                     (1/band+1/shifted+1/central)/3])
 aggregate={c:[] for c in (1.0,0.5)}
 hashes=[]
 for rep in range(4):
  seed=202610070900+10000*n+rep
  rng=np.random.default_rng(seed)
  x=rng.uniform(0,4,size=(ns,n))
  bits=np.stack([((x%p)<p/2).astype(np.float64) for p in (1.,2.,4.)])
  for c in (1.0,0.5):
   gs=[]
   for p in (1.,2.,4.):
    pos=(x%p)*table_n/p; ind=np.floor(pos).astype(np.int64); a=pos-ind
    tab=tables[c,p]
    vals=(1-a)*tab[ind]+a*tab[(ind+1)%table_n]
    gs.append(np.clip(vals,0,1))
   gg=np.ascontiguousarray(np.stack(gs))
   full=np.empty((4,ns,M+1),dtype=np.float64)
   lib.profiles(ns,n,M+1,pointer(rr),pointer(bits),pointer(gg),shift,
                central,shifted,band,pointer(full))
   assert np.isfinite(full).all()
   assert float(full.min())>=-1e-12 and (full<=heights[:,None,None]+1e-9).all()
   diagnostics+=2
   peak=full.max(axis=2); winners=full.argmax(axis=2)
   name=f'bernoulli_weak_endpoint_n{n}_c{c:g}_rep{rep}_20261007.npz'
   path=HERE/name
   assert not path.exists()
   np.savez_compressed(path,receiver=x,profile=peak,winner_index=winners,
                       r_nodes=rr,heights=heights,seed=seed,c=c)
   hashes.append({'file':name,'sha256':sha256(path.read_bytes()).hexdigest()})
   aggregate[c].append(peak)
  print(json.dumps({'progress_n':n,'completed_batch':rep+1,'seconds':time.time()-started}),flush=True)
 reports=[]
 dkw=(log(2*24/0.01)/(2*ns*4))**0.5
 for c in (1.0,0.5):
  prof=np.concatenate(aggregate[c],axis=1)
  ee={p:errors[c,p]['total_analytic_error']+1e-10 for p in (1.,2.,4.)}
  additive=np.asarray([n*heights[0]*ee[4.],n*heights[1]*ee[4.],
                       n*heights[2]*ee[4.],
                       n*(ee[1.]/band+ee[2.]/shifted+ee[4.]/central)/3])
  for model in range(4):
   hi=heights[model]
   jj=np.arange(-8,ceil(8*np.log2(hi))+9)
   lambdas=2.0**(jj/8)
   empirical=np.asarray([(prof[model]>v).mean() for v in lambdas])
   low=np.asarray([(prof[model]-additive[model]>v).mean() for v in lambdas])
   high=np.asarray([(C*(prof[model]+additive[model])>v).mean() for v in lambdas])
   lower=lambdas*np.maximum(0,low-dkw)
   upper=lambdas*np.minimum(1,high+dkw)
   upper[lambdas>=hi]=0  # exact Markov pointwise cap by input height.
   best=int(np.argmax(lambdas*empirical))
   reports.append({'c':c,'model':model,'input_height':float(hi),
                   'additive_profile_error_analytic_plus_float_screen':float(additive[model]),
                   'empirical_grid_weak_ratio':float((lambdas*empirical)[best]),
                   'empirical_best_threshold':float(lambdas[best]),
                   'simultaneous_lower':float(lower.max()),
                   'simultaneous_upper_including_between_threshold_factor':float(2**(1/8)*upper.max()),
                   'lambda_grid':lambdas.tolist(),'empirical_probabilities':empirical.tolist(),
                   'lower_ratios':lower.tolist(),'upper_ratios':upper.tolist(),
                   'confidence_is_conditional_on_float_screen_not_interval':True})
 rounds.append({'n':n,'M':M,'samples_per_batch':ns,'independent_batches':4,
                'total_receivers':ns*4,'dkw_epsilon':dkw,'r_grid_factor':'512/511',
                'source_component_masses':{'central':str(central_q),'shifted_union':str(shifted_q),
                                           'band':str(band_q)},
                'reports':reports,'saved_profiles':hashes})
result={'status':'PASS_EXACT_AND_NUMERICAL_SCREEN','exact_checks':exact_checks,
        'numerical_diagnostics':diagnostics,'elapsed_seconds':time.time()-started,
        'registration_sha256':sha256(REG.read_bytes()).hexdigest(),
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
        'cpp_sha256':sha256(CPP.read_bytes()).hexdigest(),
        'one_dimensional_errors':{f'c{c:g}_P{p:g}':v for (c,p),v in errors.items()},
        'rounds':rounds,
        'scope':'Original G_c finite-period source weak-level pressure; not a universal upper bound, not interval numerics, no actual geom history certification.'}
OUT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':result['status'],'exact_checks':exact_checks,
                  'diagnostics':diagnostics,'seconds':result['elapsed_seconds'],
                  'empirical_ratios':[[r['n'],[v['empirical_grid_weak_ratio'] for v in r['reports']]] for r in rounds]}),flush=True)
