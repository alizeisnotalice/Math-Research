"""Exact capacity and one large-n planning example, not a constructed label list."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import json

HERE=Path(__file__).resolve().parent

def gates(B):
    N=B-1;T=N*(N+1)//2;R=N*(N-1)
    H=N+T*(T-1)+2*R*R+4*R
    return dict(B=B, basis_bits=B.bit_length()-1, nonzero_labels=N,
                distinct_products=T, distinct_ratios=R,
                cyclic_no_adjacency_necessary_K=2*R+3,
                sufficient_H=H, smallest_d_with_2d_gt_H=H.bit_length())

n=10**32;p=F(1,2);d=46;K=2**d-1
v=9*p*p*n/(29*K*K)
v_next=9*p*p*n/(29*(2*K+1)**2)
# log(10) in (2,7/3), hence 16 log(10^32) in (1024,3584/3).
assert v>F(3584,3) and v_next<1024
B=1
while B<v:B*=2
M=isqrt(v.numerator//v.denominator)+1
assert (M-1)**2<v<M*M
g=gates(B)
assert K+1>g['sufficient_H'] and M<K
assert F(20,3)*p*M/K<F(1,4) # crude H_M<=M background qualification
assert K>=47858724           # previously audited global-cross sufficient K
out=dict(B512=gates(512), concrete_example=dict(
    n=str(n),p=str(p),d=d,K=str(K),v=str(v),v_display=float(v),
    B=B,M=M,q_field=str(K+1),H_B=str(g['sufficient_H']),
    next_period_v_display=float(v_next),
    original_max_period_certified_by_log10_bounds=True,
    existence_lemma_qualifies=True,
    explicit_field_polynomial_or_label_basis_constructed=False,
    checked_scope='source-level existence/normalization/background/cross gates; not a new lower bound'))
(HERE/'parameter_gate_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
