"""Exact four-coset checkerboard keys; no union search or R4 recomputation."""
from fractions import Fraction as F
from math import prod,isqrt
from vg6_bounds import Endpoint,guarded_extreme

COSETS=((0,0),(0,1),(1,0),(1,1))
SIGNS=(1,-1,-1,1)


def q(x):
    if isinstance(x,float):raise TypeError('exact rational input required')
    return F(x)


def integers(xs,n):
    xs=tuple(xs)
    if len(xs)!=n or any(isinstance(x,(bool,float)) or int(x)!=x for x in xs):
        raise ValueError(f'{n} integer coordinates required')
    return tuple(map(int,xs))


def original_to_super(index):
    k=integers(index,4);i,p=divmod(k[0],2);j,r=divmod(k[1],2)
    return 2*p+r,(i,j,k[2],k[3])


def super_to_original(coset,index):
    coset=integers((coset,),1)[0];k=integers(index,4)
    if not 0<=coset<4:raise ValueError('coset index must be 0..3')
    p,r=COSETS[coset]
    return 2*k[0]+p,2*k[1]+r,k[2],k[3]


def unpack_key(key):
    key=integers(key,32);lo=[key[4*i:4*i+4] for i in range(4)]
    hi=[key[16+4*i:20+4*i] for i in range(4)];counts=[]
    for a,b in zip(lo,hi):
        if any(y<x for x,y in zip(a,b)):
            if a!=(0,)*4 or b!=(-1,)*4:raise ValueError('noncanonical empty coset')
            counts.append(0)
        else:counts.append(prod(y-x+1 for x,y in zip(a,b)))
    return lo,hi,tuple(counts)


def split_block(lo,hi):
    lo=integers(lo,4);hi=integers(hi,4);lows=[];highs=[]
    empty=any(b<a for a,b in zip(lo,hi))
    for p,r in COSETS:
        a=tuple(-((-(lo[j]-v))//2) for j,v in enumerate((p,r)))+lo[2:]
        b=tuple((hi[j]-v)//2 for j,v in enumerate((p,r)))+hi[2:]
        if empty or any(y<x for x,y in zip(a,b)):
            lows.append((0,)*4);highs.append((-1,)*4)
        else:lows.append(a);highs.append(b)
    key=tuple(x for row in lows+highs for x in row);counts=unpack_key(key)[2]
    return dict(key=key,counts=counts,N=sum(counts),delta=sum(s*n for s,n in zip(SIGNS,counts)))


def old_residue_keys(lo,hi):
    lo=integers(lo,4);hi=integers(hi,4)
    return tuple(split_block((lo[0]+a,lo[1]+b)+lo[2:],(hi[0]+a,hi[1]+b)+hi[2:])
                 for a,b in COSETS)


def finite_mass(key,epsilon):
    epsilon=q(epsilon)
    if abs(epsilon)>=1:raise ValueError('all four masses must stay positive')
    counts=unpack_key(key)[2]
    return sum((n*(1+s*epsilon) for n,s in zip(counts,SIGNS)),F(0))


def radius_descriptor(key,alpha,eta):
    alpha=q(alpha);eta=q(eta)
    if eta<=0 or abs(alpha)*eta>=1:raise ValueError('global mass positivity not certified')
    counts=unpack_key(key)[2];N=sum(counts);delta=sum(s*n for s,n in zip(SIGNS,counts))
    info=dict(counts=counts,N=N,delta=delta,alpha=alpha,eta=eta,
              complete_supercell_mass=F(4),complete_mass_slope=F(0))
    if not N:return dict(info,regime='empty')
    c=alpha*delta/N;rho=eta*abs(c);h=1-rho
    if h<=0:raise ValueError('witness mass lower bound not positive')
    return dict(info,regime='positive_baseline',relative_mass_slope=c,
                radius_velocity_factor=c/4,mass_lower=N*h,rho=rho,
                radius_remainder_factor=3*c*c/(32*h*h))


def super_shift_direction(h):
    h=tuple(map(q,h))
    if len(h)!=4:raise ValueError('four shift components required')
    return h[0]/2,h[1]/2,h[2],h[3]


def coset_phases(h,anchor=False):
    u=super_shift_direction(h)
    return tuple(tuple((F(p,2) if j==0 else F(r,2) if j==1 else F(0))
                       +(s-(1 if anchor else 0))*u[j] for j in range(4))
                 for (p,r),s in zip(COSETS,SIGNS))


def rephase_key(key,integer_phase_additions):
    """If phase_i is replaced by phase_i+t_i, replace atom indices by k_i-t_i."""
    if len(integer_phase_additions)!=4:raise ValueError('four integer phase additions required')
    shifts=[integers(row,4) for row in integer_phase_additions]
    lo,hi,counts=unpack_key(key);lows=[];highs=[]
    for i in range(4):
        if counts[i]:
            lows.append(tuple(lo[i][j]-shifts[i][j] for j in range(4)))
            highs.append(tuple(hi[i][j]-shifts[i][j] for j in range(4)))
        else:lows.append((0,)*4);highs.append((-1,)*4)
    return tuple(x for row in lows+highs for x in row)


def compile_rational_receiver(key,P,alpha,shift_direction,period_direction,eta,anchor=False):
    """epsilon=tau*alpha, P=P0*(1+tau*kappa), d=P*tau*h; no cell clipping.

    Default phases move by s*u. anchor=True subtracts the SAME u from every
    coset, hence fixes positive cosets and moves negative cosets by -2*u.
    The compiler requires rational r0; generic signed descriptors do not.
    """
    info=radius_descriptor(key,alpha,eta);N=info['N'];root=isqrt(isqrt(N))
    if not N or root**4!=N:raise ValueError('supply certified algebraic r0 for this N')
    P=tuple(map(q,P));kappa=tuple(map(q,period_direction));eta=q(eta)
    if len(P)!=4 or len(kappa)!=4 or min(P)<=0:raise ValueError('positive four-dimensional periods required')
    Q=(2*P[0],2*P[1],P[2],P[3]);u=super_shift_direction(shift_direction)
    bases=coset_phases([0]*4);velocities=[tuple((s-(1 if anchor else 0))*x for x in u) for s in SIGNS]
    lo,hi,counts=unpack_key(key);active=[i for i,n in enumerate(counts) if n]
    r0=F(root,2);lam=info['radius_velocity_factor'];beta=info['radius_remainder_factor']
    lower=[];upper=[];guards=[]
    for j in range(4):
        hp=1-abs(kappa[j])*eta
        if hp<=0:raise ValueError('period positivity not certified')
        xmax,gmax=guarded_extreme([Endpoint(hi[i][j]+bases[i][j],velocities[i][j]) for i in active],eta,'max')
        xmin,gmin=guarded_extreme([Endpoint(lo[i][j]+bases[i][j],velocities[i][j]) for i in active],eta,'min')
        f0=r0/Q[j];f1=f0*(lam-kappa[j]);Bf=f0*(beta+abs(lam*kappa[j])+kappa[j]**2)/hp
        lower.append(Endpoint(xmax.base-f0,xmax.velocity-f1,Bf))
        upper.append(Endpoint(xmin.base+f0,xmin.velocity+f1,Bf))
        guards.append(dict(max_atom=gmax,min_atom=gmin))
    return dict(descriptor=info,Q=Q,r0=r0,r1=lam*r0,lower=lower,upper=upper,
                super_shift_velocity=u,anchor=anchor,raw_atom_guards=guards)
