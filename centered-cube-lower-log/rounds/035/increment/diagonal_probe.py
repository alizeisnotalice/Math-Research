"""Fixed diagnostic of unsmoothed full Gram diagonal; not interval certified."""
from pathlib import Path
import json,hashlib,math
import numpy as np
out=Path(__file__).resolve().parent
reg=out/'diagonal_registration.json'
if reg.exists():raise RuntimeError('Refusing repeated registered run')
plan={'seed':202610290358,'K':[63,255,1023], 'v':[8,32,128],
      'L_grid':np.linspace(10/9,11/9,129).tolist(),
      'smoothing':0,'scope':'unsmoothed Gram diagonal source diagnostic only',
      'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
reg.write_text(json.dumps(plan,indent=2))
rng=np.random.default_rng(plan['seed']);rows=[]
L=np.array(plan['L_grid']);t=L%1
for K in plan['K']:
 u=rng.random(K);a=(np.arange(K)+u)/K
 d=(a[:,None]-a[None,:])%1;d=d[d>0];d.sort();cs=np.r_[0.,np.cumsum(d)]
 ix=np.searchsorted(d,t,side='right');r0=float(np.mean(u*u-u+1/3))
 D=K*t/2-K*K*t*t/2+t*ix-cs[ix]
 hh=2*D
 for v in plan['v']:
  M=math.ceil(math.sqrt(v));m=np.arange(1,M+1)
  b=4*(729/29)*((np.sin(np.pi*m*10/27)/(np.pi*m))**4+3*(np.sin(np.pi*m/9)/(np.pi*m))**4)
  c=np.sqrt(2*b/3)
  phat=np.exp(-2j*np.pi*a[:,None]*m).sum(axis=0)/(2j*np.pi*m)
  s=np.sin(np.pi*L[:,None]*m)**2
  bb=(s*(c*c/2)).sum(axis=1)
  cross=(s*(-4*c*phat.imag)).sum(axis=1)
  gram=hh+bb+cross;target=9*L*L/29
  rows.append({'K':K,'v':v,'M':M,'R0':r0,'HH_diagonal_spread':float(np.ptp(hh)),
    'HH_max_error_from_one_third':float(abs(hh-1/3).max()),
    'cross_max_abs':float(abs(cross).max()),'full_Gram_max_error':float(abs(gram-target).max()),
    'v_times_full_Gram_max_error':float(v*abs(gram-target).max()),
    'min_Gram':float(gram.min()), 'finite_BB_tail_bound':48/(29*M**3)})
res={'numpy':np.__version__,'registration_sha256':hashlib.sha256(reg.read_bytes()).hexdigest(),
     'rows':rows,'caution':'ordinary float, unsmoothed finite L grid, no good-event certificate and no C_n lower bound'}
(out/'diagonal_results.json').write_text(json.dumps(res,indent=2))
for row in rows:print(row['K'],row['v'],'v_G_error',round(row['v_times_full_Gram_max_error'],6),'cross',round(row['cross_max_abs'],6))
