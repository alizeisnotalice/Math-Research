"""Frozen shorter-window diagnostics, stable Fourier class selection."""
from pathlib import Path
import numpy as np,json,math,hashlib
OUT=Path(__file__).resolve().parent
reg=OUT/'short_registration.json'
if reg.exists():raise RuntimeError('Refusing repeated registered run')
plan={'seed':202610290362,'K':32,'R':32,'n':[15,127,511],
      'L':[.125,.25,.5], 'centers':16,'delta':'1/(64 K R)',
      'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
reg.write_text(json.dumps(plan,indent=2));rng=np.random.default_rng(plan['seed']);K=R=32
labels=rng.permutation(K);delta=1/(64*K*R);rows=[]
for n in plan['n']:
 x=rng.random((16,n))
 for L in plan['L']:
  logits=[];modprob=[];caps=[]
  for b in range(R):
   a=(np.arange(K)+((labels+b)%R+.5)/R)/K
   w=np.zeros((16,n,K))
   for sh in [-1,0,1]:
    aa=a[None,None,:]+sh
    w+=np.maximum(0,np.minimum(aa+delta/2,x[:,:,None]+L/2)-np.maximum(aa-delta/2,x[:,:,None]-L/2))/delta/(K*L)
   S=w.sum(axis=-1);assert (S>0).all()
   char=np.fft.fft(w,axis=-1)/S[:,:,None]
   cp=np.prod(char,axis=1)
   # np FFT convention: residue b is inverse DFT at b.
   pp=np.fft.ifft(cp,axis=-1)[:,b].real
   logits.append(np.log(S).sum(axis=1));modprob.append(pp)
   caps.append((R-1)*np.max(np.prod(abs(char[:,:,1:]),axis=1),axis=1))
  logits=np.array(logits);prob=np.array(modprob);cap=np.array(caps)
  aw=np.exp(logits-logits.max(axis=0));aw/=aw.sum(axis=0)
  residual=R*(aw*(prob-1/R)).sum(axis=0)
  rows.append({'K':K,'R':R,'n':n,'L':L,'max_relative_joint_prior_residual':float(abs(residual).max()),
   'weighted_character_envelope':float((aw*cap).sum(axis=0).max()),
   'candidate_exp_envelope':(R-1)*math.exp(-n*L*L/16)})
# Separate finite Dirichlet diagnostic, using all h and N in this grid.
worst=-1;checked=0
for rr in range(2,129):
 for N in range(2,2*rr+1):
  h=np.arange(1,rr);alpha=2*np.pi*h/rr
  q=np.abs(np.sin(N*alpha/2)/(N*np.sin(alpha/2)))
  bound=math.sqrt(1-min(.5,4*(N*N-1)/(3*rr*rr)))
  worst=max(worst,float((q-bound).max()));checked+=len(h)
res={'numpy':np.__version__,'registration_sha256':hashlib.sha256(reg.read_bytes()).hexdigest(),
 'rows':rows,'Dirichlet_diagnostic':{'checked':checked,'max_ratio_minus_candidate_bound':worst},
 'scope':'ordinary float finite controls, no weak coefficient or interval certificate'}
(OUT/'short_results.json').write_text(json.dumps(res,indent=2))
for row in rows:print(row['n'],row['L'],row['max_relative_joint_prior_residual'])
print('Dirichlet',checked,worst)
