from pathlib import Path
import numpy as np,json
P=Path(__file__).resolve().parent;rows=[];K=256
sig=np.array([1.]);q=sig.copy()
while len(sig)<K:sig,q=np.r_[sig,q],np.r_[sig,-q]
for R in [1,4,16,64]:
 m=R*K+np.arange(K)
 for c in [1.,1/R**2]:
  size=2**int(np.ceil(np.log2(64*(m[-1]+1))))
  coef=np.zeros(size,dtype=complex);coef[m]=np.pi*m*np.sqrt(2*c/(3*K))*sig
  d=(np.fft.ifft(coef)*size).real/K
  for p in [.1,.5]:
   need=np.maximum(0,-(1-p)/p-d);mass=float(need.mean())
   row=dict(K=K,R=R,c=c,p=p,grid=size,needed_mass=mass,quadrature_error_bound=float(np.pi**2/(np.sqrt(3)*K*size)*np.sqrt(np.sum(m.astype(float)**4)*(c/K))),L1=float(abs(d).mean()),L2_squared=float(np.mean(d*d)),sup=float(abs(d).max()),feasible_grid=mass<=1)
   if mass<=1:
    mu=need+1-mass;rho=1-p+p*mu+p*d
    row.update(mu_mean=float(mu.mean()),rho_mean=float(rho.mean()),rho_min=float(rho.min()))
   rows.append(row)
(P/'results.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
