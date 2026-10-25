"""Three-orbit 24-key descriptors and rational local witness compilation.

General N0 uses a field-compatible factor r1=lambda*r0, not a new radical.
The rational compiler is only for N0 a perfect fourth power. Real n4 can use
the descriptor with its existing exact algebraic r0/a field implementation.
"""
from fractions import Fraction as F
from math import prod,isqrt
from vg6_bounds import Endpoint,guarded_extreme


def q(x):
    if isinstance(x,float):raise TypeError('exact rational data required')
    return F(x)


def parse_key(key):
    if len(key)!=24 or any(isinstance(x,(bool,float)) or int(x)!=x for x in key):
        raise ValueError('expected 24 integer indices')
    key=tuple(map(int,key));lo=[key[4*i:4*i+4] for i in range(3)]
    hi=[key[12+4*i:16+4*i] for i in range(3)]
    counts=[]
    for a,b in zip(lo,hi):
        if any(y<x for x,y in zip(a,b)):
            if a!=(0,0,0,0) or b!=(-1,-1,-1,-1):
                raise ValueError('empty orbit must use canonical sentinel')
            counts.append(0)
        else:counts.append(prod(y-x+1 for x,y in zip(a,b)))
    return lo,hi,tuple(counts)


def descriptor(key,weights):
    lo,hi,counts=parse_key(key);w=tuple(map(q,weights))
    if len(w)!=2 or min(w)<0:raise ValueError('two nonnegative direction masses required')
    N0,k1,k2=counts;S=w[0]*k1+w[1]*k2
    out=dict(key=tuple(key),counts=counts,weights=w,captured_mass_slope=S,
             complete_cell_mass_slope=sum(w,F(0)),regime='N0>0' if N0 else 'pure')
    if N0:
        out.update(radius_velocity_factor=S/F(4*N0),
                   radius_remainder_factor=3*S*S/F(32*N0*N0),
                   identity='r1=lambda*r0=S/(8*N0^(3/4)); Br<=beta*r0 for S>=0')
    return out


def compile_rational(key,weights,phases,periods,direction,eta,phase_velocities=None):
    """Finite common-radius receiver lower=max atoms-f, upper=min atoms+f.

    Phases are unwrapped affine coordinates. Canonical phase wrapping requires
    compensating integer key shifts; wrapping phases alone changes atoms.
    All raw atom extrema pass audit013 guards on the supplied whole cap.
    No cell clipping is done here; later max/min clipping needs its own guard.
    """
    info=descriptor(key,weights);N0=info['counts'][0]
    root=isqrt(isqrt(N0))
    if not N0 or root**4!=N0:
        raise ValueError('use descriptor/exact algebraic radius for this N0 regime')
    eta=q(eta);a=tuple(map(q,periods));d=tuple(map(q,direction))
    phi=[tuple(map(q,row)) for row in phases]
    tau=[(F(0),)*4 for _ in range(3)] if phase_velocities is None else [tuple(map(q,row)) for row in phase_velocities]
    if len(a)!=4 or len(d)!=4 or len(phi)!=3 or len(tau)!=3 or any(len(z)!=4 for z in phi+tau):
        raise ValueError('expected three phase rows and four coordinate values')
    if min(a)<=0 or eta<=0:raise ValueError('positive periods and cap required')
    lo,hi,counts=parse_key(key);active=[i for i,k in enumerate(counts) if k]
    r0=F(root,2);lam=info['radius_velocity_factor'];beta=info['radius_remainder_factor']
    lowers=[];uppers=[];guards=[];radial=[]
    for j in range(4):
        h=1-abs(d[j])*eta
        if h<=0:raise ValueError('period positivity not certified')
        xmax,gmax=guarded_extreme([Endpoint(hi[i][j]+phi[i][j],tau[i][j]) for i in active],eta,'max')
        xmin,gmin=guarded_extreme([Endpoint(lo[i][j]+phi[i][j],tau[i][j]) for i in active],eta,'min')
        f0=r0/a[j];f1=f0*(lam-d[j])
        Bf=f0*(beta+abs(lam*d[j])+d[j]**2)/h
        lowers.append(Endpoint(xmax.base-f0,xmax.velocity-f1,Bf))
        uppers.append(Endpoint(xmin.base+f0,xmin.velocity+f1,Bf))
        guards.append(dict(max_atom=gmax,min_atom=gmin))
        radial.append(dict(f0=f0,f1=f1,B_f=Bf))
    return dict(descriptor=info,r0=r0,r1=lam*r0,lower=lowers,upper=uppers,
                radial=radial,raw_atom_guards=guards,eta=eta)
