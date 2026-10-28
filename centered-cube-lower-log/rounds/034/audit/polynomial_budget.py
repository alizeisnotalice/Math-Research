"""Exact integer/Fraction gates for v >= n**beta, n=2**m.

No source simulation. Rational upper bounds use pi<22/7, sum m^-2<2,
H_M<=1+ceil(log2 M), and sqrt(K+1)<=ceil_sqrt(K+1).
"""
from fractions import Fraction as F
from pathlib import Path
import json,math

BETAS=[F(1,16),F(1,12),F(1,10),F(1,9),F(1,8),F(1,6),F(1,5),F(1,4),F(1,3)]
def H(B):
    z=B-1;t=z*(z+1)//2;r=z*(z-1)
    return z+t*(t-1)+2*r*r+4*r
def l2(x):return math.log2(x.numerator)-math.log2(x.denominator)
rows=[]
for m in [720,1440,2880]:
    n=1<<m
    for beta in BETAS:
        assert (beta*m).denominator==1
        target=1<<int(beta*m)
        limit=math.isqrt(9*n//(116*target))
        d=(limit+1).bit_length()-1; K=(1<<d)-1
        v=F(9*n,116*K*K)
        assert v>=target and F(9*n,116*((1<<(d+1))-1)**2)<target
        B=1
        while B<v:B*=2
        M=math.isqrt(v.numerator//v.denominator)
        if M*M*v.denominator<v.numerator:M+=1
        hbound=1+(M-1).bit_length()
        sq=math.isqrt(K+1)
        if sq*sq<K+1:sq+=1
        cx=F(40,9)*(F(2*sq,K)+F(11*hbound,7*K))
        face=v/F(1<<(3*m//8))+F(M,1<<(m//2))+F(M,1<<(4*m))+F(M*math.isqrt(m+1)+M,1<<(m//8))
        budgets=dict(a_n_v=F(77,30)*v/K,v_squared_over_n=v*v/n,
                     v_Cx_upper=v*cx,safe_v_squared_over_K=v*v/K,
                     v_T_M=F(48,29)*v/M**3,
                     v_etaN=v*(F(2,3*B*B)+F(2,K)),
                     fixed_sigma_face_upper=face)
        ratios=(B-1)*(B-2)
        rows.append(dict(n_power2=m,beta=str(beta),target_v_power2=int(beta*m),
                         field_degree=d,K=str(K),v=str(v),B_power2=B.bit_length()-1,
                         M=str(M),q_gt_H=(K+1>H(B)),
                         necessary_no_adjacent_ratio_capacity=(K>=2*ratios+1),
                         M_lt_K=(M<K),v_ge_16logn_sufficient=(v>=16*m),
                         background_positive_sufficient=F(1,2)-F(10*hbound,3*K)>0,
                         log2_field_margin=l2(F(K+1,H(B))),
                         log2_ratio_capacity_margin=l2(F(K,2*ratios+1)),
                         log2_budgets={a:l2(b) for a,b in budgets.items()},
                         exact_budgets={a:str(b) for a,b in budgets.items()},
                         exact_cancellation_squared=str(F(K*K)*v/n)))
        assert F(K*K)*v/n==F(9,116)
out=dict(scope='Exact scalar parameter gates/upper budgets, not full finite-n proof or explicit field labels.',
         p='1/2',rows=rows)
Path(__file__).with_name('polynomial_budget_results.json').write_text(json.dumps(out,indent=2)+'\n')
for r in rows:
    if r['n_power2'] in [720,2880]:
        print(json.dumps({k:r[k] for k in ['n_power2','beta','field_degree','B_power2','q_gt_H','necessary_no_adjacent_ratio_capacity','log2_field_margin','log2_budgets']},sort_keys=True))
