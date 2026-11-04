"""Exact finite certificates for the n-dependent parameter hierarchy.
The binomial argument is an analytic lemma; astronomical exponential powers
are not evaluated. No runtime limits and no floating arithmetic.
"""
from fractions import Fraction as F
from math import factorial
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json

def ceilf(x): return -(-x.numerator//x.denominator)
def check(n):
    c=F(1,2**n); c1=20**n
    c2=2*n*4**n*3**(n-1); boundary=2*c2+2
    alpha=c/(4*boundary)
    lower=F(2**(2*n+3)*factorial(n+1)*c1,1)/alpha**(n+1)
    K=max(2*(n+1),ceilf(lower))
    eta=F(1,1000*(K+1))
    D=ceilf(F(16*n*c1*(K+1)**n,1)/c)
    N=(3000*(K+1))**n  # (3/eta)^n, an integer half-box count majorant
    h=eta/(16*N)
    eps=(2*h)**n
    v=eps/(4*D)  # only certifies the explicit covering threshold
    checks={
        "alpha_range":0<alpha<1,
        "K_at_least_2n_plus2":K>=2*(n+1),
        "binomial_sufficient_bound":F(K,1)>=lower,
        "geometric_small_window":0<eta<F(1,40),
        "graph_depth_radius":20*eta*(K+1)<=F(1,4),
        "internal_volume_payment":F(2*n*c1*(K+1)**n,D)<=c/4,
        "thickening_step":h<eta/2 and 4*N*h<=eta,
        "covering_defect_threshold":D*v<eps,
        "v_smaller_than_eta":v<eta
    }
    assert all(checks.values())
    return {"n":n,"checks":checks,"K_bits":K.bit_length(),
       "D_bits":D.bit_length(),"eta_denominator_bits":eta.denominator.bit_length(),
       "v_denominator_bits":v.denominator.bit_length()}
if __name__=="__main__":
    with ProcessPoolExecutor(max_workers=15) as pool:
        rows=list(pool.map(check,range(2,17)))
    out={"scope":"n=2..16 exact hierarchy certificates; not a substitute for the all-n geometric proof",
         "formulas":"c=2^-n,C1=20^n,C2=2n*4^n*3^(n-1),C*=2C2+2,alpha=c/(4C*); K>=2^(2n+3)(n+1)!C1 alpha^(-n-1), K>=2(n+1); eta=1/[1000(K+1)]; D>=16nC1(K+1)^n/c",
         "all_pass":all(all(r["checks"].values()) for r in rows),"cases":rows}
    Path(__file__).with_suffix(".json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({"cases":len(rows),"all_pass":out["all_pass"],
        "max_K_bits":max(r["K_bits"] for r in rows),
        "max_v_denominator_bits":max(r["v_denominator_bits"] for r in rows)}))
