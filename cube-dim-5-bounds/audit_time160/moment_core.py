"""Prepared exact time-moment factor core. No source values run on import.

Input J bounds must already be qualified for the same complete offset profile
and the inverse of the supplied physical time. No R/J/H evaluator is imported.
"""
from fractions import Fraction as F
from math import comb

def Q(x):
    if isinstance(x,(float,bool)):raise TypeError('Exact rational input required')
    return F(x)

def direct_weights(z_hi,z_lo):
    a,b=Q(z_hi),Q(z_lo);assert 0<=b<a
    E={j:F(4,j+4)*(a**(j+4)-b**(j+4))/(a**4-b**4) for j in (1,2,3)}
    alpha=(E[3]-b*E[2])/(a-b);beta=(a*E[2]-E[3])/(a-b)
    assert alpha>0 and beta>0 and alpha+beta==E[2]
    return dict(E1=E[1],E2=E[2],E3=E[3],alpha=alpha,beta=beta)

def positive_weights(z_hi,z_lo):
    """Independent affine-u integration,using positive binomial coefficients."""
    a,b=Q(z_hi),Q(z_lo);assert 0<=b<a;d=a-b
    integral_power=lambda n:sum((F(comb(n,k),k+1)*b**(n-k)*d**k for k in range(n+1)),F(0))
    D=integral_power(3);assert D>0
    E={j:integral_power(j+3)/D for j in (1,2,3)}
    alpha=sum((F(comb(5,k),k+2)*b**(5-k)*d**k for k in range(6)),F(0))/D
    beta=sum((F(comb(5,k),(k+1)*(k+2))*b**(5-k)*d**k for k in range(6)),F(0))/D
    assert alpha>0 and beta>0
    return dict(E1=E[1],E2=E[2],E3=E[3],alpha=alpha,beta=beta)

def ceil_grid(x,bits=192):
    x=Q(x);D=1<<bits;v=x*D
    return F(-((-v.numerator)//v.denominator),D)

def factor_from_exact_bounds(time_left,time_right,profile,J_left_upper,
                             J_right_upper,old_factor_upper,range_upper=None):
    """Future authorized finite case only. None righttime means the full tail.

    profile is the actual111/112 profile; J_left_upper=J(1/time_left)_U.
    old_factor_upper is an already qualified192-grid band factor. Source
    binding to geometry/offset hash,time keys,catalog hash is a caller gate.
    No mean or stored band mass is evaluated here.
    """
    tl=Q(time_left);assert tl>0
    tr=None if time_right is None else Q(time_right)
    assert tr is None or tr>tl
    a,b=1/tl,F(0) if tr is None else 1/tr
    weights=direct_weights(a,b);assert weights==positive_weights(a,b)
    mu,J0=Q(profile['mu_infinity']),Q(profile['J0']);assert mu>=0 and J0>=0
    A=Q(J_left_upper);assert A>=J0
    if tr is None:
        assert J_right_upper is None;B=J0
    else:
        B=Q(J_right_upper);assert B>=J0
    used_B=min(B,A) # Upper bound remains valid even if saved B>A.
    trial=1+mu*weights['E1']+F(25,24)*(weights['alpha']*A+weights['beta']*used_B)
    old=Q(old_factor_upper);D=1<<192
    assert old>=1 and (old*D).denominator==1 and trial>=1
    choices=[trial,old]
    if range_upper is not None:
        r=Q(range_upper);assert r>=1;choices.append(r)
    chosen=min(choices);rounded=ceil_grid(chosen)
    assert chosen<=rounded<=old
    return dict(**weights,J_left_upper=A,J_right_original_upper=B,J_right_used_upper=used_B,
                taylor_moment_upper=trial,selected_exact_upper=chosen,factor_upper=rounded,
                time_left=tl,time_right=tr,z_hi=a,z_lo=b,J0=J0,mu_infinity=mu)
