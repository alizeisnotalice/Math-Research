"""Scalar planning/error budgets; no cube or Gaussian-process simulation.

Log comparisons and displayed Cx are float diagnostics.  Finite-field
capacity H(B) and v are computed with Python integers/Fraction.
"""
from pathlib import Path
from fractions import Fraction as F
import json,math

HERE=Path(__file__).resolve().parent
p=F(1,2);cp=3*float(p)/math.sqrt(29);w=4

def H(B):
    z=B-1;t=z*(z+1)//2;r=z*(z-1)
    return z+t*(t-1)+2*r*r+4*r

rows=[]
for power in [12,24,32,48,64,96,160]:
    n=10**power;ln=power*math.log(10)
    for exponent in [1,2]:
        threshold=16*ln**exponent
        d=1
        while float(9*p*p*n/(29*(2**(d+1)-1)**2))>=threshold:d+=1
        K=2**d-1;v=9*p*p*n/(29*K*K);vf=float(v)
        B=1
        while B<v:B*=2
        M=math.isqrt(v.numerator//v.denominator)+1
        hm=sum(1/j for j in range(1,M+1));sm2=sum(1/j**2 for j in range(1,M+1))
        Cx=(40/9)*((1-1/B)*math.sqrt(K+1)/K*sm2
                         +math.pi*(1-1/B)**2/(2*K)*hm)
        an=(77/15)*float(p)/K
        u=math.sqrt(vf)
        face=vf*math.exp(-3*ln/8)+u/math.sqrt(n)+u*math.exp(-4*ln)+u*math.exp(-ln/8)*math.sqrt(ln)
        rows.append(dict(log_power=exponent,n=str(n),n_power10=power,d=d,K=str(K),
            v=str(v),v_display=vf,B=B,M=M,H_B=str(H(B)),
            field_existence_q_gt_H=(K+1>H(B)),
            amplitude_times_v=an*vf,v_squared_over_n=vf*vf/n,
            v_Cx=vf*Cx,v_T_M=vf*48/(29*M**3),
            v_etaN=vf*(2/(3*B*B)+2/K),
            defect_count_over_K=((B-1)*(B-2))/K,
            v_squared_over_K=vf*vf/K,
            face_comparison_budget=face,
            log_far_pair_budget=math.log(K*vf)-vf/25,
            single_point_Mills_prefactor=1/u,
            natural_cells_Kv_over_sqrt_n=K*vf/math.sqrt(n),
            first_sum_Ksqrtv_over_sqrt_n=K*u/math.sqrt(n),
            exact_cancellation_reference=cp,
            h=vf/2-2,log_lambda=vf/2-2))

# Exact first-order cancellation across all plans (to float print accuracy).
assert max(abs(r['first_sum_Ksqrtv_over_sqrt_n']-cp) for r in rows)<1e-14
(HERE/'parameter_results.json').write_text(json.dumps(dict(
    p=str(p),w=w,note='planning diagnostics, not a finite-n lower-bound certificate',rows=rows),indent=2)+'\n')
for r in rows:
    print({k:r[k] for k in ['n_power10','log_power','d','v_display','B','field_existence_q_gt_H','amplitude_times_v','v_Cx','defect_count_over_K','face_comparison_budget','log_far_pair_budget','first_sum_Ksqrtv_over_sqrt_n']})
