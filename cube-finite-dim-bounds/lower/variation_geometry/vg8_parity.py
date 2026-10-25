"""Exact original-lattice parity split into the two-coset superperiod.

Only finite selected witnesses are described. This module does not enumerate
receiver unions, certify old R4 coverage, or search perturbations.
"""
from fractions import Fraction as F
from math import prod,isqrt
from vg6_bounds import Endpoint,guarded_extreme


def q(x):
    if isinstance(x,float):raise TypeError('exact rational data required')
    return F(x)


def integers(xs,n):
    xs=tuple(xs)
    if len(xs)!=n or any(isinstance(x,(bool,float)) or int(x)!=x for x in xs):
        raise ValueError(f'{n} integer coordinates required')
    return tuple(map(int,xs))


def parity_counts(L,U):
    L,U=integers((L,U),2)
    if U<L:return 0,0
    even=U//2-(L-1)//2
    return even,U-L+1-even


def original_to_super(index):
    k=integers(index,4);quotient,parity=divmod(k[0],2)
    return parity,(quotient,)+k[1:]


def super_to_original(parity,index):
    if parity not in (0,1):raise ValueError('parity must be 0 or 1')
    k=integers(index,4)
    return (2*k[0]+parity,)+k[1:]


def unpack_key(key):
    key=integers(key,16)
    lo=[key[:4],key[4:8]];hi=[key[8:12],key[12:16]];counts=[]
    for a,b in zip(lo,hi):
        if any(y<x for x,y in zip(a,b)):
            if a!=(0,)*4 or b!=(-1,)*4:raise ValueError('noncanonical empty coset')
            counts.append(0)
        else:counts.append(prod(y-x+1 for x,y in zip(a,b)))
    return lo,hi,tuple(counts)


def split_block(lo,hi):
    lo=integers(lo,4);hi=integers(hi,4)
    lows=[];highs=[]
    empty=any(b<a for a,b in zip(lo,hi))
    for p in (0,1):
        first_lo=-((-(lo[0]-p))//2)  # ceil((lo1-p)/2), including negative indices
        first_hi=(hi[0]-p)//2
        if empty or first_hi<first_lo:
            lows.append((0,)*4);highs.append((-1,)*4)
        else:
            lows.append((first_lo,)+lo[1:]);highs.append((first_hi,)+hi[1:])
    key=tuple(x for row in lows+highs for x in row)
    ne,no=unpack_key(key)[2]
    return dict(key=key,N_even=ne,N_odd=no,N=ne+no,delta=ne-no)


def old_residue_keys(lo,hi):
    """Both old-period translation residues are necessary in the superperiod."""
    lo=integers(lo,4);hi=integers(hi,4)
    shifted_lo=(lo[0]+1,)+lo[1:];shifted_hi=(hi[0]+1,)+hi[1:]
    return split_block(lo,hi),split_block(shifted_lo,shifted_hi)


def finite_mass(key,epsilon):
    e=q(epsilon)
    if abs(e)>=1:raise ValueError('both orbit masses must remain strictly positive')
    ne,no=unpack_key(key)[2]
    return ne*(1+e)+no*(1-e)


def radius_descriptor(key,mass_direction,eta):
    alpha=q(mass_direction);eta=q(eta)
    if eta<=0 or abs(alpha)*eta>=1:raise ValueError('global masses 1 +/- tau*alpha need positivity')
    ne,no=unpack_key(key)[2];N=ne+no;delta=ne-no
    out=dict(N_even=ne,N_odd=no,N=N,delta=delta,alpha=alpha,eta=eta,
             complete_supercell_mass=F(2),complete_mass_slope=F(0))
    if not N:
        return dict(out,regime='empty')
    c=alpha*delta/N;rho=abs(c)*eta;h=1-rho
    if h<=0:raise ValueError('witness mass positivity not certified')
    return dict(out,regime='positive_baseline',relative_mass_slope=c,
                radius_velocity_factor=c/4,rho=rho,mass_lower=N*h,
                radius_remainder_factor=3*c*c/(32*h*h),
                identity='r1=lambda*r0; |rem|<=beta*r0*tau^2, beta uses h^-2 >= h^-7/4')


def relative_shift_to_super(h):
    h=tuple(map(q,h))
    if len(h)!=4:raise ValueError('four shift components required')
    return (h[0]/2,)+h[1:]


def compile_rational_receiver(key,P,mass_direction,shift_direction,period_direction,eta):
    """Exact rational baseline r0, signed mass direction and affine normalized phase.

    P(tau)=P0*(1+tau*kappa), epsilon=tau*alpha, d(tau)=P(tau)*tau*h.
    Q=(2P1,P2,P3,P4), odd superphase=(1/2,0,0,0)+tau*(h1/2,h2,h3,h4).
    Later receiver clipping and union/period-image completeness are caller duties.
    """
    info=radius_descriptor(key,mass_direction,eta);N=info['N']
    root=isqrt(isqrt(N))
    if not N or root**4!=N:raise ValueError('use signed descriptor with existing algebraic r0 for this N')
    P=tuple(map(q,P));kappa=tuple(map(q,period_direction));eta=q(eta)
    if len(P)!=4 or len(kappa)!=4 or min(P)<=0:raise ValueError('positive four-dimensional periods required')
    Q=(2*P[0],)+P[1:];u=relative_shift_to_super(shift_direction)
    phases=[(F(0),)*4,(F(1,2),F(0),F(0),F(0))]
    velocities=[(F(0),)*4,u]
    lo,hi,counts=unpack_key(key);active=[i for i,n in enumerate(counts) if n]
    r0=F(root,2);lam=info['radius_velocity_factor'];beta=info['radius_remainder_factor']
    lower=[];upper=[];guards=[]
    for j in range(4):
        hp=1-abs(kappa[j])*eta
        if hp<=0:raise ValueError('period positivity not certified')
        xmax,gmax=guarded_extreme([Endpoint(hi[i][j]+phases[i][j],velocities[i][j]) for i in active],eta,'max')
        xmin,gmin=guarded_extreme([Endpoint(lo[i][j]+phases[i][j],velocities[i][j]) for i in active],eta,'min')
        f0=r0/Q[j];f1=f0*(lam-kappa[j])
        Bf=f0*(beta+abs(lam*kappa[j])+kappa[j]**2)/hp
        lower.append(Endpoint(xmax.base-f0,xmax.velocity-f1,Bf))
        upper.append(Endpoint(xmin.base+f0,xmin.velocity+f1,Bf))
        guards.append(dict(max_atom=gmax,min_atom=gmin))
    return dict(descriptor=info,Q=Q,r0=r0,r1=lam*r0,lower=lower,upper=upper,
                super_shift_velocity=u,raw_atom_guards=guards)
