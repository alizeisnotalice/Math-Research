"""Frozen macroscopic-window modular-source diagnostics, no weak-norm claim."""
from pathlib import Path
import numpy as np,json,hashlib
OUT=Path(__file__).resolve().parent
reg=OUT/'registration.json'
if reg.exists():raise RuntimeError('Refusing repeat registered run')
plan={'seed':202610290361,'K':8,'R':[2,4,8],'n':[3,7,15,63],
      'L':[10/9,7/6,11/9],'centers_per_case':128,'delta':'1/(64 K R)',
      'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
reg.write_text(json.dumps(plan,indent=2));rng=np.random.default_rng(plan['seed']);rows=[]
K=plan['K']
for R in plan['R']:
 labels=rng.integers(R,size=K);delta=1/(64*K*R)
 atoms=(np.arange(K)[None,:]+((labels[None,:]+np.arange(R)[:,None])%R+.5)/R)/K
 for n in plan['n']:
  x=rng.random((128,n))
  for L in plan['L']:
   terms=[];products=[];eta_products=[];extra_max=0.
   for b in range(R):
    coeff=np.zeros((128,R));coeff[:,0]=1;prod=np.ones(128);ep=np.ones(128)
    for i in range(n):
     w=np.zeros((128,K))
     for sh in [-1,0,1]:
      a=atoms[b][None,:]+sh
      w+=np.maximum(0,np.minimum(a+delta/2,x[:,i,None]+L/2)-np.maximum(a-delta/2,x[:,i,None]-L/2))/delta/(K*L)
     z=np.stack([w[:,r::R].sum(axis=1)for r in range(R)],axis=1)
     S=z.sum(axis=1);F=np.fft.fft(z,axis=1)
     ep*=np.abs(F[:,1:]).max(axis=1)/S;prod*=S
     coeff=sum(z[:,r,None]*np.roll(coeff,r,axis=1)for r in range(R))
     extra_max=max(extra_max,float((S*L-1).max()))
    terms.append(coeff[:,b]);products.append(prod);eta_products.append(ep)
   joint=np.sum(terms,axis=0);prior=np.mean(products,axis=0)
   rel=float(np.max(np.abs(joint/prior-1)))
   empirical_envelope=float((R-1)*np.max(eta_products))
   coarse=(R-1)*((L-1+2/K)/(L+2/K))**n
   rows.append({'K':K,'R':R,'n':n,'L':L,'max_relative_joint_prior_difference':rel,
     'Fourier_product_envelope':empirical_envelope,'candidate_coarse_envelope':coarse,
     'max_extra_period_mass':extra_max,'relative_error_to_Fourier_envelope_floor':rel/max(empirical_envelope,1e-14)})
res={'numpy':np.__version__,'registration_sha256':hashlib.sha256(reg.read_bytes()).hexdigest(),'rows':rows,
 'scope':'ordinary floating-point query diagnostics; apparent floor at 1e-14 is not mathematical residual'}
(OUT/'results.json').write_text(json.dumps(res,indent=2))
for n in plan['n']:
 a=[r for r in rows if r['n']==n];print('n',n,'max observed',max(r['max_relative_joint_prior_difference']for r in a),'max envelope',max(r['candidate_coarse_envelope']for r in a))
