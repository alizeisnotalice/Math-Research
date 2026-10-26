"""153 exact arithmetic interface. No H/gradient/log/corner calculation.

Caller must qualify the same-model, same-point H_(2eps) lower L and H_eps
upper U, and any old q upper, in the physical Euclidean-coordinate metric.
Arithmetic alone cannot establish those source/domain qualifications.
"""
from fractions import Fraction as F

def exact(x):
    if isinstance(x,(float,bool)):raise TypeError('Exact rational input required')
    return F(x)

def curvature_caps(epsilon,mass_lower,mass_upper,old_q_upper=None):
    e,L,U=map(exact,(epsilon,mass_lower,mass_upper))
    # Strict positivity of p_t implies H_2eps < H_eps, hence valid L < U.
    if not (e>0 and U>0 and 0<=L<U):raise ValueError('Require eps>0,U>0,0<=L<U')
    A=1/(16*e*e);r=L/U
    uniform=F(2,3)*A;split=A*(F(7,10)-F(8,15)*r)
    old=None if old_q_upper is None else exact(old_q_upper)
    if old is not None and old<=0:raise ValueError('Qualified old q must be positive')
    new=min(uniform,split) if old is None else min(old,uniform,split)
    effective=max(F(1,16),r)
    assert min(uniform,split)==A*(F(7,10)-F(8,15)*effective)
    assert 0<new<=uniform
    # No assertion q>=A/4, q>=256, or q>=A/6 is valid for the final min.
    return {'A':A,'mass_ratio_lower':r,'effective_mass_ratio_lower':effective,
        'uniform_time_cap':uniform,'split_time_mass_cap':split,
        'old_q_upper':old,'new_q_upper':new}
