"""Finite mask/Fourier cancellation diagnostics, deterministic ordinary float."""
from pathlib import Path
import numpy as np,json,hashlib
OUT=Path(__file__).resolve().parent;rows=[]
for K in [12,63,255]:
 for R in [r for r in range(2,min(K,16)+1)if K%r==0]:
  m=np.arange(1,K+1);j=np.arange(K)
  cell=np.exp(2j*np.pi*m[:,None]*(j+.5)/K)*np.sin(np.pi*m[:,None]/K)/(np.pi*m[:,None])
  masses=np.stack([cell[:,r::R].sum(axis=1)for r in range(R)],axis=1)
  selected=m<K//R
  rows.append({'K':K,'R':R,'first_possible_nonzero_mode':K//R,
   'max_abs_below_threshold':float(abs(masses[selected]).max())if selected.any()else 0,
   'max_abs_at_threshold':float(abs(masses[K//R-1]).max())})
res={'numpy':np.__version__,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'rows':rows,'scope':'finite deterministic diagnostic, exact cancellation proved separately'}
(OUT/'class_mode_results.json').write_text(json.dumps(res,indent=2)+'\n')
print('cases',len(rows),'max below threshold',max(r['max_abs_below_threshold']for r in rows))
