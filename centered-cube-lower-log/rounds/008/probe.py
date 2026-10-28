from pathlib import Path
import numpy as np,json
P=Path(__file__).resolve().parent
rows=[]
for K in [64,256,1024]:
 signs=np.array([1.]);q=signs.copy()
 while len(signs)<K:signs,q=np.r_[signs,q],np.r_[signs,-q]
 for R in [1,4,16]:
  m=R*K+np.arange(K);c=1/R**2;b=np.full(K,c/K)
  # Oversampled full periodic grid; FFT evaluates real trigonometric sum.
  size=2**int(np.ceil(np.log2(64*(m[-1]+1))))
  coeff=np.zeros(size,dtype=complex);coeff[m]=np.pi*m*np.sqrt(2*b/3)*signs
  d=(np.fft.ifft(coeff)*size).real/K
  coh=np.pi*np.sum(m*np.sqrt(2*b/3))/K
  L=np.linspace(10/9,11/9,1001);J=L+1/(K*R)
  diff=(np.sin(np.pi*J[:,None]*m)/J[:,None]-np.sin(np.pi*L[:,None]*m)/L[:,None])*np.sqrt(b/3)
  rows.append(dict(K=K,R=R,c=c,grid_size=size,grid_sup=float(abs(d).max()),grid_min=float(d.min()),coherent_sup=float(coh),increment_sq_max=float(np.max(np.sum(diff**2,axis=1))),theory_sup=float(np.pi*np.sqrt(4*c/3)*(R+1))))
(P/'results.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
