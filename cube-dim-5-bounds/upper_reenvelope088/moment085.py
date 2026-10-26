"""Exact moments for product-uniform boxes, including fixed coordinates.

No heat kernel, band mass, gradient, transcendental evaluation or sampling.
All public numeric outputs are fractions.Fraction values.
"""
from fractions import Fraction as F
from itertools import combinations


def _multiply(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i+j] += a*b
    return out


def max_abs_square_moment(offsets, with_pieces=False):
    """Return E max_j |U_j|^2 for independent U_j uniform on [lo_j,hi_j].

    Zero-width coordinates are fixed point masses, not division by zero.
    With no coordinates, max is defined as zero. with_pieces returns a
    certificate of each exact CDF polynomial and its integrated contribution.
    """
    spans = [(F(lo), F(hi)) for lo, hi in offsets]
    if any(lo > hi for lo,hi in spans):
        raise ValueError('Reversed coordinate interval')
    knots = sorted({F(0)} | {abs(x) for span in spans for x in span})
    total = F(0)
    pieces = []
    for left, right in zip(knots, knots[1:]):
        mid = (left + right)/2
        cdf = [F(1)]
        for lo,hi in spans:
            if lo == hi:
                if mid < abs(lo):
                    cdf = [F(0)]
                    break
                continue
            # min(hi,r)-max(lo,-r), selected on the open knot interval.
            top = (hi, F(0)) if hi < mid else (F(0), F(1))
            bottom = (lo, F(0)) if lo > -mid else (F(0), F(-1))
            c0, c1 = top[0]-bottom[0], top[1]-bottom[1]
            if c0+c1*mid <= 0:
                cdf = [F(0)]
                break
            cdf = _multiply(cdf, [c0/(hi-lo), c1/(hi-lo)])
        contribution = right**2-left**2
        contribution -= sum((2*c*(right**(k+2)-left**(k+2))/F(k+2)
                             for k,c in enumerate(cdf)), F(0))
        assert contribution >= 0
        total += contribution
        pieces.append({'lo':left, 'hi':right, 'cdf_coefficients':cdf,
                       'moment_contribution':contribution})
    assert 0 <= total <= (knots[-1]**2 if knots else 0)
    return (total,pieces) if with_pieces else total


def box_statistics(bounds, anchor_tail):
    """Four-tail box statistics relative to an arbitrary fixed anchor.

    The anchor need not be the center or lie inside the box. 'volume' is
    ordinary four-dimensional volume; 'support_measure' is the product of
    positive widths for the lower-dimensional probability interpretation.
    A zero-volume box contributes zero to a four-dimensional integral.
    """
    if len(bounds) != 4 or len(anchor_tail) != 4:
        raise ValueError('Expected exactly four tail coordinates')
    bounds = [(F(lo),F(hi)) for lo,hi in bounds]
    anchor_tail = list(map(F,anchor_tail))
    if any(lo > hi for lo,hi in bounds):
        raise ValueError('Reversed coordinate interval')
    offsets = [(lo-c,hi-c) for (lo,hi),c in zip(bounds,anchor_tail)]
    mu = [(lo+hi)/2 for lo,hi in offsets]
    second = [(lo*lo+lo*hi+hi*hi)/3 for lo,hi in offsets]
    eta = sum(second,F(0))/4
    eta -= (sum(second,F(0))+2*sum((mu[i]*mu[j] for i,j in combinations(range(4),2)),F(0)))/36
    m2,pieces = max_abs_square_moment(offsets,with_pieces=True)
    volume = F(1)
    support_measure = F(1)
    support_dimension = 0
    for lo,hi in bounds:
        volume *= hi-lo
        if hi > lo:
            support_measure *= hi-lo
            support_dimension += 1
    old_coefficient = F(21,8)*eta
    proposed_coefficient = F(9,16)*m2
    face_coefficient = F(25,48)*m2
    assert face_coefficient <= F(20,21)*old_coefficient
    return {'offsets':offsets, 'mu':mu, 'mu_infinity':max(map(abs,mu)),
            'h_infinity':max(abs(z) for span in offsets for z in span),
            'E_max_abs_squared':m2, 'eta_old_Q':eta,
            'old_curvature_coefficient':old_coefficient,
            'proposed_9_over_8_coefficient':proposed_coefficient,
            'face_25_over_24_coefficient':face_coefficient,
            'new_to_old_coefficient_ratio':face_coefficient/old_coefficient if old_coefficient else F(0),
            'volume':volume,'support_dimension':support_dimension,
            'support_measure':support_measure,'cdf_pieces':pieces}


def band_mean_factor(time_lower, ratio_upper, stats):
    """Legal factor multiplying a positive band mass upper endpoint.

    ratio_upper must already bound band(c+s*u)/band(c) for every box offset
    and 0<=s<=1; this function does not construct or evaluate that ratio.
    The range, new face Taylor and old Taylor factors are all retained.
    """
    a, R = F(time_lower), F(ratio_upper)
    if a <= 0 or R < 1:
        raise ValueError('Need time_lower>0 and an outward ratio_upper>=1')
    linear = stats['mu_infinity']/a
    new = 1+linear+stats['face_25_over_24_coefficient']*R/a**2
    old = 1+linear+stats['old_curvature_coefficient']*R/a**2
    return min(R,new,old)
