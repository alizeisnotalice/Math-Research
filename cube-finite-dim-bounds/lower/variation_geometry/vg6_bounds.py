"""Explicit finite-epsilon error constants; all inputs must be certified rational bounds."""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from math import prod


def q(x):
    if isinstance(x,float):
        raise TypeError('Float inputs are not exact certificates')
    return F(x)


@dataclass(frozen=True)
class Endpoint:
    base: F
    velocity: F
    remainder: F = F(0)

    def __post_init__(self):
        for key in ('base','velocity','remainder'):
            object.__setattr__(self,key,q(getattr(self,key)))
        if self.remainder < 0:
            raise ValueError('remainder bound must be nonnegative')


def guarded_extreme(candidates,eta,kind='max'):
    """Certify a max/min endpoint expansion on the whole [0,eta].

    Every raw candidate must already have its own valid Endpoint bound on this
    interval. All candidates are required, including clipping constants. The
    returned remainder is valid even when baseline-tied candidates switch.
    Distinct-baseline candidates must pass a strict all-interval guard.
    """
    candidates=list(candidates);eta=q(eta)
    if not candidates or eta<=0 or kind not in ('max','min'):
        raise ValueError('nonempty candidates, eta>0, kind=max/min required')
    op=max if kind=='max' else min
    base=op(e.base for e in candidates)
    active=[e for e in candidates if e.base==base]
    velocity=op(e.velocity for e in active)
    carrier_B=min(e.remainder for e in active if e.velocity==velocity)
    output_B=max(e.remainder for e in active)
    guards=[]
    for e in candidates:
        if e.base==base:
            continue
        gap=(base-e.base) if kind=='max' else (e.base-base)
        slope=(e.velocity-velocity) if kind=='max' else (velocity-e.velocity)
        movement=max(F(0),slope)*eta+(e.remainder+carrier_B)*eta**2
        if gap<=movement:
            raise ValueError('raw candidate may change the extreme branch on [0,eta]')
        guards.append(dict(gap=gap,movement_upper=movement,strict_slack=gap-movement))
    return Endpoint(base,velocity,output_B),guards


def subsets_at_least_two():
    for k in range(2,5):
        yield from combinations(range(4),k)


def geometry_constants(boxes,tmax):
    """Boxes are (lower_endpoints,upper_endpoints); epsilon <= tmax**4.

    Returns universal error C*epsilon**2 after any separately certified V0,D.
    It does not enumerate union cells or certify the supplied derivative.
    Every Endpoint bound must hold on this whole interval. A clipped endpoint
    needs guarded_extreme at a cap >= eta, or a separately proved true bound;
    the output-endpoint gap check here cannot see omitted raw branches.
    """
    tmax=q(tmax)
    if not 0<tmax<=1 or not boxes:
        raise ValueError('nonempty finite family and 0<tmax<=1 required')
    eta=tmax**4
    axes=[[side[j] for box in boxes for side in box] for j in range(4)]
    if any(len(side)!=4 for box in boxes for side in box):
        raise ValueError('four-coordinate endpoints required')
    ranges=[];variations=[];speeds=[];remainders=[];gaps=[];heights=[]
    for ends in axes:
        cuts=sorted(set((e.base,e.velocity) for e in ends))
        bases=sorted(set(e.base for e in ends))
        R=bases[-1]-bases[0]
        T=sum((abs(b[1]-a[1]) for a,b in zip(cuts,cuts[1:])),F(0))
        V=max(abs(e.velocity) for e in ends)
        B=max(e.remainder for e in ends)
        delta=min((b-a for a,b in zip(bases,bases[1:])),default=None)
        if delta is not None and 4*V*eta>delta:
            raise ValueError('epsilon cap does not preserve distinct-base order')
        ranges.append(R);variations.append(T);speeds.append(V);remainders.append(B)
        gaps.append(delta);heights.append(R+2*eta*V+2*eta**2*B)
    Caff=sum((eta**(len(S)-2)*prod(variations[j] if j in S else ranges[j]
                                   for j in range(4))
              for S in subsets_at_least_two()),F(0))
    Crem=sum((sum(e.remainder for e in axes[j])
              *prod(heights[k] for k in range(4) if k!=j)
              for j in range(4)),F(0))
    return dict(eta=eta,tmax=tmax,ranges=ranges,velocity_variation=variations,
                max_speeds=speeds,max_remainders=remainders,min_distinct_gaps=gaps,
                ambient_widths=heights,C_affine=Caff,C_remainder=Crem,C=Caff+Crem,
                endpoint_speed_sums=[sum(abs(e.velocity)+e.remainder*eta for e in ends)
                                     for ends in axes])


def period_constants(direction,eta):
    d=list(map(q,direction));eta=q(eta)
    if len(d)!=4:
        raise ValueError('four period directions required')
    factors=[1-abs(x)*eta for x in d]
    if min(factors)<=0:
        raise ValueError('period positivity not certified')
    CJ=sum((eta**(len(S)-2)*prod(abs(d[j]) for j in S)
            for S in subsets_at_least_two()),F(0))
    return dict(j1=sum(d,F(0)),C_J=CJ,j_min=prod(factors))


def geometry_from_caps(ranges,velocity_variation,max_speeds,max_remainders,
                       remainder_sums,min_distinct_gaps,endpoint_speed_sums,tmax):
    """Same constants from certified rational caps, e.g. for algebraic n4 input.

    R,T,V,B,Bsum,S are upper bounds; each non-None gap is a strictly positive
    lower bound. T may conservatively be 2*sum(abs(endpoint velocities)).
    Algebraic order ranks must not be substituted for metric coordinate bounds.
    """
    arrays=[list(map(q,x)) for x in (ranges,velocity_variation,max_speeds,
                                   max_remainders,remainder_sums,endpoint_speed_sums)]
    if any(len(x)!=4 or min(x)<0 for x in arrays) or len(min_distinct_gaps)!=4:
        raise ValueError('four nonnegative coordinate caps required')
    R,T,V,B,BS,S=arrays;tmax=q(tmax);eta=tmax**4
    if not 0<tmax<=1:
        raise ValueError('0<tmax<=1 required')
    gaps=[None if d is None else q(d) for d in min_distinct_gaps]
    if any(d is not None and (d<=0 or 4*V[j]*eta>d) for j,d in enumerate(gaps)):
        raise ValueError('distinct-base endpoint order not certified')
    H=[R[j]+2*eta*V[j]+2*eta**2*B[j] for j in range(4)]
    CA=sum((eta**(len(A)-2)*prod(T[j] if j in A else R[j] for j in range(4))
            for A in subsets_at_least_two()),F(0))
    CR=sum((BS[j]*prod(H[k] for k in range(4) if k!=j) for j in range(4)),F(0))
    return dict(eta=eta,tmax=tmax,ranges=R,velocity_variation=T,max_speeds=V,
                max_remainders=B,min_distinct_gaps=gaps,ambient_widths=H,
                C_affine=CA,C_remainder=CR,C=CA+CR,endpoint_speed_sums=S)


def pure_constants(geometry,width_coefficients,center_speed_bounds,rho,
                   volume_coefficient,j1_abs=0,C_J=0,j_min=1):
    """One actual singleton (or exactly colocated mass), certified cone theta.

    Actual normalized side lengths <= w_j*epsilon**(1/4), center shifts <= P_j*epsilon.
    Ideal unshielded volume is volume_coefficient*epsilon/j(epsilon).
    rho_j is a lower bound on distance to every nonincident old endpoint in axis j,
    or None if no such endpoint exists. Pairwise pure-cube separation is a caller duty.
    Returns K*epsilon**(7/4)+C_period*epsilon**2.
    """
    w=list(map(q,width_coefficients));P=list(map(q,center_speed_bounds))
    if len(w)!=4 or len(P)!=4 or len(rho)!=4 or min(w+P)<0:
        raise ValueError('four nonnegative width/speed bounds required')
    eta=geometry['eta'];tmax=geometry['tmax']
    for j in range(4):
        if rho[j] is None:
            continue
        distance=q(rho[j])
        movement=geometry['max_speeds'][j]+geometry['max_remainders'][j]*eta+P[j]
        if distance<=0 or w[j]*tmax/2+movement*eta>distance/2:
            raise ValueError('pure cube leaves certified local-cone neighborhood')
    K=sum(((geometry['endpoint_speed_sums'][j]+2*P[j])
           *prod(w[k] for k in range(4) if k!=j) for j in range(4)),F(0))
    j_min=q(j_min);c=q(volume_coefficient);jA=q(j1_abs);CJ=q(C_J)
    if j_min<=0 or min(c,jA,CJ)<0:
        raise ValueError('positive Jacobian and nonnegative coefficient bounds required')
    Cp=c*(jA+CJ*eta)/j_min
    return dict(K=K,C_period=Cp)


def certify_gain(V0_upper,D_abs,g_lower,C,K,tmax,
                 j1_abs=0,C_J=0,m1_abs=0,C_M=0,ratio_scale=1):
    """j*V - V0*m >= g*epsilon - Cnet*epsilon**2 - K*epsilon**(7/4).

    g_lower must already enclose D+V0*(j1-m1), including all pure terms.
    j=J/J0, m=M/M0; ratio_scale=J0/M0. No floating derivatives are accepted.
    """
    V0,D,g,C,K,tmax,jA,CJ,mA,CM,scale=map(q,
        (V0_upper,D_abs,g_lower,C,K,tmax,j1_abs,C_J,m1_abs,C_M,ratio_scale))
    if min(V0,D,C,K,jA,CJ,mA,CM)<0 or g<=0 or scale<=0 or not 0<tmax<=1:
        raise ValueError('positive certified gain and nonnegative bounds required')
    eta=tmax**4
    if min(1-jA*eta-CJ*eta**2,1-mA*eta-CM*eta**2)<=0:
        raise ValueError('normalized mass/Jacobian positivity not certified')
    H=V0+D*eta+C*eta**2+K*tmax**7
    Cnet=C+jA*(D+C*eta+K*tmax**3)+CJ*H+V0*CM
    s=tmax
    while Cnet*s**4+K*s**3>g/2:
        s/=2
    epsilon=s**4
    normalized_margin=g*epsilon-Cnet*epsilon**2-K*s**7
    ratio_margin=scale*normalized_margin/(1+mA*epsilon+CM*epsilon**2)
    assert normalized_margin>=g*epsilon/2>0
    return dict(epsilon=epsilon,fourth_root=s,C_net=Cnet,K_net=K,
                normalized_margin_lower=normalized_margin,
                ratio_improvement_lower=ratio_margin)
