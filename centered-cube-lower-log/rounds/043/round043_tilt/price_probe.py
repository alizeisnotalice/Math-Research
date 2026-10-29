"""Exact2x2 Gram and numerical leading-eigenvalue price before formal target."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,math
import numpy as np
import exact_tilt_example as e
HERE=Path(__file__).resolve().parent

def main():
 target=HERE/'PRICE.json';assert not target.exists();saved=json.loads((HERE/'EXACT_TILT_EXAMPLE.json').read_text());q=1;eta=F(1,2);m=F(saved['m_t_exact']);G=[[F(0)for _ in range(2)]for _ in range(2)]
 for part in saved['partitions']:
  h=F(part['right'])-F(part['left']);z=[[F(v)for v in row]for row in part['z_endpoints']]
  for i in range(2):
   for j in range(2):G[i][j]+=e.affine_product(h,z[0][i],z[1][i],z[0][j],z[1][j])/eta/m
 H=np.array([[float(v)for v in row]for row in G]);ev=np.linalg.eigvalsh(H);gamma=float(ev[-1]);small=float(ev[0]);mu=[F(v)for v in saved['mu_t_exact']]
 chi=sum(v*v/eta for v in mu)-1;kappa=F(1,q*(q+1));ratio=(gamma-1)/float(chi)
 assert 1-1e-12<=ratio<=1/float(1-kappa)+1e-12
 rows=[]
 for n in [1,2,4,8,32,1000,1000000]:
  r=small/gamma;logvalue=n*math.log(gamma)+math.log1p(-r**(n+1))-math.log1p(-r)
  rows.append(dict(n=n,log_wrong_prior_trace=logvalue,log_leading_price=n*math.log(gamma),normalized_prefactor=math.exp(logvalue-n*math.log(gamma))))
 out=dict(script_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),exact_source_sha=hashlib.sha256((HERE/'EXACT_TILT_EXAMPLE.json').read_bytes()).hexdigest(),
   original_prior_Gram_exact=[[str(v)for v in row]for row in G],original_prior_Gram_float=H.tolist(),
   eigenvalues=ev.tolist(),chi2_prior_shift_exact=str(chi),gamma_minus1_over_chi2=ratio,kappa_exact=str(kappa),rows=rows,
   scope='Auxiliary spectral price of original-prior normalization, not C_n growth')
 target.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(eigenvalues=ev.tolist(),ratio=ratio,chi2=float(chi))))
if __name__=='__main__':main()
