"""Exact endpoint interface only. No roots, unions, grids or floating arithmetic.

Scalars a,s,r0 and rates may belong to a caller's exact ordered algebraic field.
The caller supplies the certified equality root r0=N^(1/5)/2. For finite parameters,
the caller supplies radius(mass), with its exact/inner/outer status recorded separately.
"""
from fractions import Fraction as Q
from itertools import product
from math import prod

STARTS = (-1, 0, 1, 2)


def parity(j):
    return 1 if j % 2 == 0 else -1


def tail_indices(ell):
    start = -((ell-1)//2)
    return tuple(range(start, start+ell))


def phases_extrema(start, ell, displacement):
    if ell == 1:
        z = parity(start)*displacement
        return z, z
    return -abs(displacement), abs(displacement)


def mass(lengths, start1, mu):
    N = prod(lengths)
    coefficient = parity(start1)*prod(lengths[1:]) if lengths[0] % 2 else 0
    return N+mu*coefficient, coefficient


def atom(j, a, s, eps, eta, mu):
    return ((a*j[0]+eps*parity(j[1]), s*j[1]+eta*parity(j[0]),
             s*j[2], s*j[3], s*j[4]), 1+mu*parity(j[0]))


def extrema(lengths, starts, a, s, eps, eta):
    l1,l2,*tail = lengths
    j1,j2 = starts
    e0,e1 = phases_extrema(j2,l2,eps)
    h0,h1 = phases_extrema(j1,l1,eta)
    ans = [(a*j1+e0,a*(j1+l1-1)+e1),
           (s*j2+h0,s*(j2+l2-1)+h1)]
    for ell in tail:
        ids = tail_indices(ell)
        ans.append((s*ids[0],s*ids[-1]))
    return tuple(ans)


def receiver_domain(a,s):
    return ((0,2*a),(0,s),(0,s/2),(0,s/2),(0,s/2))


def descriptors():
    """Exactly 4^5*4^2=16384 raw descriptions; caller accounts for every one."""
    for lengths in product(range(1,5),repeat=5):
        for starts in product(STARTS,repeat=2):
            yield lengths,starts


def finite_receiver(lengths, starts, a, s, eps, eta, mu, radius):
    m,_ = mass(lengths,starts[0],mu)
    assert m > 0
    r = radius(m)
    assert r > 0
    extremes = extrema(lengths,starts,a,s,eps,eta)
    box = tuple((max(dl,hi-r),min(du,lo+r))
                for (lo,hi),(dl,du) in zip(extremes,receiver_domain(a,s)))
    return {'lengths':lengths,'starts':starts,'mass':m,'radius':r,
            'extrema':extremes,'box':box,'positive_volume_possible':all(lo<hi for lo,hi in box)}


def raw_receiver_jet(lengths, starts, a, s, r0, rates):
    """Right jet for (eps,eta,mu)=t*(e,h,m), t down to 0 positive.

    Each endpoint is (constant,speed). Lexicographic clipping retains all ties.
    A zero constant width with positive speed is kept, including possible births.
    All raw descriptors must be visited before direction-dependent pruning.
    """
    e,h,m = rates
    N = prod(lengths)
    _,coefficient = mass(lengths,starts[0],0)
    r1 = r0*Q(coefficient,5*N)*m
    base = extrema(lengths,starts,a,s,0,0)
    e0,e1 = phases_extrema(starts[1],lengths[1],e)
    h0,h1 = phases_extrema(starts[0],lengths[0],h)
    speeds = ((e0,e1),(h0,h1),(0,0),(0,0),(0,0))
    box=[]
    for (lo,hi),(vlo,vhi),(dl,du) in zip(base,speeds,receiver_domain(a,s)):
        box.append((max((dl,0),(hi-r0,vhi-r1)),min((du,0),(lo+r0,vlo+r1))))
    widths = tuple((hi[0]-lo[0],hi[1]-lo[1]) for lo,hi in box)
    return {'lengths':lengths,'starts':starts,'N':N,'mass_speed':m*coefficient,
            'radius_jet':(r0,r1),'extrema_constant':base,'extrema_speed':speeds,
            'box':tuple(box),'width_jets':widths,
            'strictly_negative_constant_width':any(w[0]<0 for w in widths),
            'lex_positive_volume':all(w>(0,0) for w in widths)}
