"""Exact parameter budgets and finite-atom weighted-mixture/max decomposition checks."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
import json,hashlib,math
import numpy as np
HERE=Path(__file__).resolve().parent

def weak(values):
 # Sup over strict thresholds equals the left limit at an attained positive value.
 return max(t*F(sum(x>=t for x in values),len(values)) for t in set(values))

def main():
 alpha=F(261,3388);assert alpha>F(1,16)
 p=F(99,100);n=10**32;K=2**47-1;B=2048;v=9*p*p*n/(29*K*K)
 e=alpha*v;eps=1/v;eta_phase=F(1,64*v)
 var=F(B*B-1,6*B*B)-eta_phase*F(B-1,6*B)
 assert var>=F(1,7)
 kappa=p/F(K)*F(2+F(40,9))/F(10,9);assert kappa<1
 J=math.ceil(K*v)+1
 with localcontext() as dc:
  dc.prec=65
  dv=Decimal(v.numerator)/v.denominator;de=Decimal(e.numerator)/e.denominator
  assert de>Decimal(n).ln()
  parameter_rows=[]
  for R in [4,math.floor(v),math.floor(v*v)]:
   # sqrt(v) is upper-bounded rationally; exp(-E) <= 2^-floor(E).
   sqrtv_upper=F(math.isqrt(v.numerator//v.denominator)+1)
   strict_fee=F(J*(R-1))*sqrtv_upper/F(2**math.floor(e))
   n_fee=F(J*(R-1))*sqrtv_upper/n
   H=1+R*R*(B-1)**2
   bad=min(J,H*(math.floor(4*v)+1))
   parameter_rows.append(dict(R=R,J_upper=J,E_lower_exact=str(e),floor_E_lower=math.floor(e),
     good_fee_power2_upper=str(strict_fee),good_fee_numeric=float(strict_fee),coarse_n_inverse_fee=float(n_fee),
     bad_grid_count_upper=bad,bad_over_sqrt_n=float(F(bad,10**16)),
     note='Arithmetic budget only; no polynomial-R family is constructed here.'))
 # Atom models: all columns have exact mass1; good scales use disjoint peak blocks;
 # bad scales deliberately coincide across banks and are paid only by their count.
 rng=np.random.default_rng(500001);m=36;R=3;Js=3;Jb=3;epsilon=F(1,16);peak=1000;S=12*peak*peak+24
 weights_cases=[[F(1,3)]*3,[F(1,9),F(4,9),F(4,9)]];checks=[]
 for trial in range(8):
  raw=[]
  for l in range(Js+Jb):
   perm=rng.permutation(m);cols=[]
   for b in range(R):
    z=[1]*m
    block=perm[12*b:12*(b+1)] if l<Js else perm[:12]
    for j in block:z[int(j)]=peak*peak
    cols.append(z)
   raw.append(cols)
  A=[[[F(m*z,S) for z in col] for col in cols] for cols in raw]
  affinities=[F(sum(math.isqrt(raw[l][b][j]*raw[l][c][j]) for j in range(m)),S)
               for l in range(Js) for b in range(R) for c in range(b+1,R)]
  eta=max(affinities)
  for weights in weights_cases:
   delta=4*eta*(2 if weights[0]==F(1,3) else F(16,9))
   G=[[sum(weights[b]*A[l][b][j] for b in range(R)) for j in range(m)] for l in range(Js+Jb)]
   Mgood=[max(G[l][j] for l in range(Js)) for j in range(m)]
   Mbad=[max(G[l][j] for l in range(Js,Js+Jb)) for j in range(m)]
   Mall=[max(Mgood[j],Mbad[j]) for j in range(m)]
   productweak=[weak([max(A[l][b][j] for l in range(Js+Jb)) for j in range(m)]) for b in range(R)]
   rhs=(1+epsilon)*sum(weights[b]*productweak[b] for b in range(R))+Js*delta+Jb
   assert weak(Mall)<=rhs and weak(Mbad)<=Jb
   # Direct same-threshold event decomposition, independent of weak-norm algebra.
   for lam in sorted(set(Mall+Mgood+Mbad+[F(1,2),F(1),F(2)])):
    left=lam*F(sum(x>lam for x in Mall),m)
    split=lam*F(sum(x>lam for x in Mgood)+sum(x>lam for x in Mbad),m)
    assert left<=split
   for l in range(Js):
    E=[]
    for b in range(R):E.append([sum(weights[c]*A[l][c][j] for c in range(R) if c!=b)<=epsilon*weights[b]*A[l][b][j] for j in range(m)])
    assert all(sum(E[b][j] for b in range(R))<=1 for j in range(m))
    bad=[not any(E[b][j] for b in range(R)) for j in range(m)]
    badmass=sum(G[l][j] for j in range(m) if bad[j])/m;assert badmass<=delta
   checks.append(dict(trial=trial,weights=list(map(str,weights)),affinity=str(eta),delta_per_good_scale=str(delta),
                       full_max_weak=str(weak(Mall)),bad_max_weak=str(weak(Mbad)),comparison_rhs=str(rhs),passed=True))
 out=dict(scope='Exact Fraction parameter and finite-atom checks; asymptotic theorem proved separately, no C_n improvement.',
          alpha=str(alpha),v=str(v),var=str(var),E_lower=str(e),epsilon=str(eps),parameter_rows=parameter_rows,
          finite_atom_models=checks,seed=500001,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 (HERE/'RESULTS.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'finite_models_passed':len(checks),'E_lower':float(e),'parameter_rows':parameter_rows}))
if __name__=='__main__':main()
