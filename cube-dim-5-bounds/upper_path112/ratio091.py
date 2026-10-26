"""Exact outward ratio bound for body16 Cauchy shifts with u1=0.

No H, heat-kernel, stencil, exp, log, asinh or old-ratio evaluation.
The caller supplies the true box h_infinity relative to its fixed anchor.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt


def _fraction(x):
    if isinstance(x,float):
        raise TypeError('Use an exact Fraction, rational string or integer')
    return F(x)


def chain_lift(tail):
    """Return five body columns and exact coefficients for (0,tail).

    Absolute values sort decreasingly, with original index breaking ties.
    The sign of zero is fixed as +1. Zero coefficients are retained.
    """
    if len(tail)!=4:
        raise ValueError('Expected four tail coordinates')
    u=list(map(_fraction,tail))
    order=sorted(range(4),key=lambda j:(-abs(u[j]),j))
    signs=[1 if z>=0 else -1 for z in u]
    weights=[abs(u[j]) for j in order]+[F(0)]
    radius=weights[0]
    direction=[1]+[-sign for sign in signs]
    columns=[tuple(direction)]
    coefficients=[-radius/2]
    gaps=[]
    for k,j in enumerate(order):
        direction[j+1]=signs[j]
        columns.append(tuple(direction))
        gap=weights[k]-weights[k+1]
        gaps.append(gap)
        coefficients.append(gap/2)
    assert all(g>=0 for g in gaps) and sum(gaps)==radius
    assert sum(coefficients)==0 and sum(map(abs,coefficients))==radius
    image=[sum((a*v[j] for a,v in zip(coefficients,columns)),F(0)) for j in range(5)]
    assert image==[F(0)]+u
    return {'tail':u,'order':order,'signs':signs,'radius':radius,
            'columns':columns,'coefficients':coefficients,'gaps':gaps}


def _sqrt_interval(q,bits):
    D=1<<bits
    scaled=q.numerator*D*D//q.denominator
    n=isqrt(scaled)
    lo,hi=F(n,D),F(n+1,D)
    assert lo*lo<=q<hi*hi
    return lo,hi


@lru_cache(maxsize=4096)
def _ratio_interval(x,bits):
    if x==0:
        return F(1),F(1)
    lo1,hi1=_sqrt_interval(1+x*x/16,bits)
    lo2,hi2=_sqrt_interval(1+x*x/256,bits)
    lower=(lo1+x/4)**2*(lo2+x/16)**8
    upper=(hi1+x/4)**2*(hi2+x/16)**8
    D=1<<bits
    lo=F((lower*D).numerator//(lower*D).denominator,D)
    hi=F((upper*D).numerator//(upper*D).denominator+1,D)
    assert 1<=lo<=lower<upper<hi
    return lo,hi


def ratio_interval(x,bits=192):
    """Enclose Rnew(x) on an outward dyadic grid, x=h_infinity/time_lower.

    Rnew(x)=(sqrt(1+x^2/16)+x/4)^2*(sqrt(1+x^2/256)+x/16)^8.
    Mathematical comparisons must use these endpoints, not floats.
    """
    x=_fraction(x)
    if x<0 or not isinstance(bits,int) or bits<32:
        raise ValueError('Need x>=0 and integer bits>=32')
    return _ratio_interval(x,bits)


def new_ratio_upper(h_infinity,time_lower,bits=192):
    """Valid for a fixed-anchor four-tail box and every t>=time_lower>0."""
    h,a=_fraction(h_infinity),_fraction(time_lower)
    if h<0 or a<=0:
        raise ValueError('Need h_infinity>=0 and time_lower>0')
    return ratio_interval(h/a,bits)[1]


def combine_ratio_upper(h_infinity,time_lower,existing_upper_bounds=(),bits=192):
    """Take min with existing legal upper bounds for the SAME box/anchor/time.

    Existing bounds are preserved and never recalculated. Their validity and
    matching model/geometry remain the caller's responsibility.
    """
    old=tuple(map(_fraction,existing_upper_bounds))
    if any(value<1 for value in old):
        raise ValueError('A supplied two-sided ratio cap must be >=1')
    return min((new_ratio_upper(h_infinity,time_lower,bits),)+old)
