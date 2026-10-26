"""Exact weight geometry and monotone-Darboux path-mean certificates.

No heat kernel, jet or ratio is computed here. All R values are supplied by
a qualified scalar catalog; this separates geometry, scalar work and reuse.
"""
from fractions import Fraction as F


def Q(x):
    if isinstance(x,float):raise TypeError('Exact rational input required')
    return F(x)


def mul(p,q):
    out=[F(0)]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q):out[i+j]+=a*b
    return out


def integral(poly,left,right):
    return sum((c*(right**(k+1)-left**(k+1))/F(k+1) for k,c in enumerate(poly)),F(0))


def weight_profile(offsets):
    """Piecewise polynomial W(rho)=E(max|U|-rho)+, plus J0=Emax|U|²/2.

    Independent uniform intervals may have zero width (fixed coordinates).
    Breakpoints are physical distances relative to the supplied anchor.
    """
    spans=[(Q(a),Q(b)) for a,b in offsets]
    if len(spans)!=4 or any(a>b for a,b in spans):raise ValueError('Four ordered intervals required')
    knots=sorted({F(0)}|{abs(x) for ab in spans for x in ab})
    pieces=[]
    for left,right in zip(knots,knots[1:]):
        mid=(left+right)/2;cdf=[F(1)]
        for lo,hi in spans:
            if lo==hi:
                if mid<abs(lo):cdf=[F(0)];break
                continue
            top=(hi,F(0)) if hi<mid else (F(0),F(1))
            bot=(lo,F(0)) if lo>-mid else (F(0),F(-1))
            c0,c1=top[0]-bot[0],top[1]-bot[1]
            if c0+c1*mid<=0:cdf=[F(0)];break
            cdf=mul(cdf,[c0/(hi-lo),c1/(hi-lo)])
        survival=[-c for c in cdf];survival[0]+=1
        pieces.append({'lo':left,'hi':right,'cdf':cdf,'survival':survival})
    tail=F(0)
    for piece in reversed(pieces):
        left,right=piece['lo'],piece['hi'];s=piece['survival']
        W=[tail+sum((c*right**(k+1)/F(k+1) for k,c in enumerate(s)),F(0))]
        W.extend(-c/F(k+1) for k,c in enumerate(s))
        piece['W']=W
        mass=integral(W,left,right);assert mass>=0
        piece['weight_integral']=mass
        tail+=integral(s,left,right)
    J0=sum((p['weight_integral'] for p in pieces),F(0))
    # Independent layer-cake expression for the same E r² / 2.
    M2half=sum((integral([F(0),*p['survival']],p['lo'],p['hi']) for p in pieces),F(0))
    assert J0==M2half
    mu=[(lo+hi)/2 for lo,hi in spans]
    return {'offsets':spans,'h':knots[-1],'pieces':pieces,'J0':J0,
            'M2':2*J0,'mu_infinity':max(map(abs,mu))}


def darboux_plan(offsets,time_lower,panels_per_piece=8):
    a=Q(time_lower)
    if a<=0 or not isinstance(panels_per_piece,int) or panels_per_piece<1:
        raise ValueError('Positive time and integer panels required')
    profile=weight_profile(offsets);cells=[]
    for piece in profile['pieces']:
        lo,hi=piece['lo'],piece['hi']
        for i in range(panels_per_piece):
            left=lo+(hi-lo)*F(i,panels_per_piece)
            right=lo+(hi-lo)*F(i+1,panels_per_piece)
            weight=integral(piece['W'],left,right);assert weight>=0
            if weight:
                cells.append({'left':left,'right':right,'x_left':left/a,
                              'x_right':right/a,'weight':weight})
    assert sum((c['weight'] for c in cells),F(0))==profile['J0']
    nodes=sorted({F(0),profile['h']/a}|{c[k] for c in cells for k in ('x_left','x_right')})
    return {'profile':profile,'a':a,'cells':cells,'ratio_nodes':nodes}


def evaluate_plan(plan,catalog,old_mean_factor=None):
    """Catalog maps exact str(x) to {lower,upper}, enclosing091 R_new(x).

    Every positive weight is exact. Old ratio values are consumed unchanged.
    A zero-radius plan uses the exact identity R(0)=1 and needs no catalog.
    Returned factors multiply the anchor band mass upper endpoint.
    """
    profile=plan['profile'];a=plan['a'];J0=profile['J0'];h=profile['h']
    if h==0:
        return {'J_lower':F(0),'J_upper':F(0),'J0':F(0),'old_J_upper':F(0),
                'factor_upper':F(1),'old085_factor_upper':F(1),'curvature_ratio_upper':F(0),
                'darboux_gap':F(0),'strict_curvature_improvement':False}
    def row(x):
        if x==0:return F(1),F(1)
        obj=catalog[str(x)];lo,hi=Q(obj['lower']),Q(obj['upper'])
        if not 1<=lo<=hi:raise ValueError('Invalid supplied ratio interval')
        return lo,hi
    lower=sum((c['weight']*row(c['x_left'])[0] for c in plan['cells']),F(0))
    upper=sum((c['weight']*row(c['x_right'])[1] for c in plan['cells']),F(0))
    R=row(h/a)[1];oldJ=J0*R
    lower=max(lower,J0);upper=min(upper,oldJ);assert lower<=upper
    linear=1+profile['mu_infinity']/a
    old085=min(R,linear+F(25,24)*oldJ/a**2)
    factor=min(R,linear+F(25,24)*upper/a**2,old085)
    if old_mean_factor is not None:
        old=Q(old_mean_factor)
        if old<0:raise ValueError('Invalid old mean factor')
        factor=min(factor,old)
    return {'J_lower':lower,'J_upper':upper,'J0':J0,'old_J_upper':oldJ,
            'factor_upper':factor,'old085_factor_upper':old085,
            'curvature_ratio_upper':upper/oldJ if oldJ else F(0),
            'darboux_gap':upper-lower,'strict_curvature_improvement':upper<oldJ}
