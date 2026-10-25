"""One locked 8-orbit rejection; independent tail primitives and quartic PV.

No imports from the primary checker, its Lb formula, or any Green cone code.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import permutations, product
import argparse, hashlib, json, time
import sympy as sy
from flint import arb, ctx

HERE = Path(__file__).resolve().parent

def ball(x):
    x = F(x)
    return arb(x.numerator) / x.denominator

def ends(x):
    return [str(x.lower().fmpq()), str(x.upper().fmpq())]

def inside(x, interval):
    lo, hi = map(F, ends(x)); a, b = map(F, interval)
    return a < lo <= hi < b

def dot(x, y):
    return sum(a*b for a,b in zip(x,y))

def orbit(seed):
    # Full product first; deduplicate only equal physical centers within this orbit.
    return sorted({tuple(s*t for s,t in zip(signs,p))
                   for p in permutations(seed) for signs in product((-1,1),repeat=3)})

V = [(F(1),F(j),F(k)) for j,k in product((-1,1),repeat=2)]

def unit_direction(y, v):
    """Integrate b'(s)/s over its complete line support, with PV at zero."""
    e, d = dot(y,y), dot(y,v)
    k = 1-e
    discr = d*d+3*k
    if discr <= 0:
        return arb(0), 'miss_or_tangent'
    root = ball(discr).sqrt()
    lo, hi = (-ball(d)-root)/3, (-ball(d)+root)/3
    poly = (ball(8*d*d-12*k)*(hi-lo)
            +ball(18*d)*(hi*hi-lo*lo)+12*(hi**3-lo**3))
    if k*d:
        assert not lo.contains(0) and not hi.contains(0)
        poly -= ball(4*k*d)*(abs(hi).log()-abs(lo).log())
    # If k=0, an endpoint is zero but its logarithmic coefficient is exactly zero.
    # If d=0, the same coefficient vanishes, without numerical limiting branches.
    return poly/arb.pi(), ('interior' if k>0 else 'boundary' if k==0 else 'exterior_intersection')

def tail_closed(e, d, ratio, sigma):
    E,D = ball(e),ball(d)
    discr = ball(3*e-d*d)
    assert discr > 0
    delta = discr.sqrt()
    q = (ball(d*d)+3*(ratio-E)).sqrt()
    upper = (q+sigma*D)/(ratio-E)
    assert upper > 0
    def primitive(u):
        den = E*u*u+2*sigma*D*u+3
        assert den > 0
        return ((1-ratio/E)*u + ratio*sigma*D/(E*E)*den.log()
                +ratio*ball(3*e-2*d*d)/(E*E*delta)
                 *((E*u+sigma*D)/delta).atan())
    return primitive(upper)-primitive(arb(0)), upper

def run(root):
    started=time.monotonic()
    snapshot=json.loads((HERE/'028_source_snapshot.json').read_text())
    for name,h in snapshot.items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==h,(name,'source drift')
    # Symbolic checks concern our independently chosen primitives, not the main code.
    u,E,D,R=sy.symbols('u E D R',real=True)
    delta=sy.sqrt(3*E-D*D)
    for sigma in (-1,1):
        den=E*u*u+2*sigma*D*u+3
        primitive=((1-R/E)*u+R*sigma*D/E**2*sy.log(den)
                   +R*(3*E-2*D*D)/(E**2*delta)*sy.atan((E*u+sigma*D)/delta))
        assert sy.cancel(sy.diff(primitive,u)-(1-R*u*u/den))==0
    s,k,d=sy.symbols('s k d',real=True)
    polynomial=(k-2*d*s-3*s*s)**2
    quotient,remainder=sy.div(sy.diff(polynomial,s),s,s)
    assert sy.expand(quotient-(8*d*d-12*k+36*d*s+36*s*s))==0
    assert remainder==-4*k*d
    data=json.loads((root/'results/round010c_source_input.json').read_text())
    assert len(data['point'])==3 and len(data['weights'])==len(data['specs'])==8
    point=list(map(F,data['point']));e=dot(point,point)
    assert e>0 and point==[F('211389546841383/250000000000000'),F('4485601305961609/10000000000000000'),F('4188889041543007/25000000000000000')]
    ctx.prec=256
    c=(4*arb(2).log()-1)/(24*arb.pi()**3)
    lower_constant=arb(5)/(16*arb.pi()**3)
    ratio=lower_constant/c
    assert ball(e)<ratio
    tails=[];tail_balls=[]
    for direction in V:
        dv=dot(point,direction)
        assert 3*e-dv*dv>0
        for sigma in (-1,1):
            val,upper=tail_closed(e,dv,ratio,sigma)
            assert val>0
            tail_balls.append(val)
            tails.append({'direction':list(map(str,direction)),'sigma':sigma,
                          'reciprocal_upper':ends(upper),'integral_interval':ends(val)})
    tail_sum=sum(tail_balls,arb(0))
    mu=c*tail_sum/arb.pi()
    contributions=[];contribution_balls=[];allcenters=0;effectivecenters=0;counts={}
    for weight,spec in zip(data['weights'],data['specs']):
        weight=F(weight);radius=F(spec['radius']);seed=tuple(map(F,spec['seed']))
        assert weight>=0 and radius>0
        centers=orbit(seed);allcenters+=len(centers)
        if weight:effectivecenters+=len(centers)
        subtotal=arb(0);localcounts={}
        for center in centers:
            y=[(x-a)/radius for x,a in zip(point,center)]
            for direction in V:
                term,kind=unit_direction(y,direction)
                localcounts[kind]=localcounts.get(kind,0)+1
                counts[kind]=counts.get(kind,0)+1
                subtotal+=term/ball(radius)
        weighted=ball(weight)*subtotal
        contribution_balls.append(weighted)
        contributions.append({'seed':list(map(str,seed)),'radius':str(radius),
                              'weight':str(weight),'center_count':len(centers),
                              'line_classes':localcounts,'weighted_Lb_interval':ends(weighted),
                              'weighted_Lb_ball':weighted.str(55)})
    assert [z['center_count'] for z in contributions]==[1,6,1,8,24,12,6,1]
    assert allcenters==59 and effectivecenters==51
    LB=sum(contribution_balls,arb(0))
    source_upper=mu-LB
    assert F(ends(source_upper)[1])<F(-202142,10**8)
    parent=json.loads((root/'results/round010c_source_certificate.json').read_text())
    assert parent['input_sha256']==snapshot['results/round010c_source_input.json']
    assert parent['source_sha256']==snapshot['src/round010c_source_certificate.py']
    assert inside(mu,parent['mu_upper_expression_interval'])
    assert inside(LB,parent['LB_interval'])
    assert inside(source_upper,parent['corrected_density_upper_interval'])
    return {'status':'PASS independent strict negative source for this fixed rational eight-orbit correction',
            'source_sha256':snapshot,'point':list(map(str,point)),'norm_squared':str(e),
            'precision_bits':ctx.prec,'fixed_c':'(4log(2)-1)/(24pi^3)',
            'Jensen_lower_inherited_from_027':'G(y)>=5/(16pi^3 |y|^2); no rerun',
            'point_in_Jensen_interior':True,'eight_tails':tails,
            'mu_upper_expression_interval':ends(mu),'mu_upper_expression_ball':mu.str(55),
            'basis_contributions':contributions,'expanded_center_count':allcenters,
            'nonzero_weight_center_count':effectivecenters,'all_line_classes':counts,
            'LB_interval':ends(LB),'LB_ball':LB.str(55),
            'corrected_density_upper_interval':ends(source_upper),
            'corrected_density_upper_ball':source_upper.str(55),
            'strict_source_upper_bound':'-202142/100000000',
            'all_three_intervals_strictly_inside_primary_enclosures':True,
            'implementation':'8 exact antiderivative evaluations; full orbit expansion and independent quartic derivative/PV endpoint primitive; no primary code imports',
            'shared_libraries':{'SymPy':sy.__version__,'Arb':'python-flint, 256bit','Fraction':'stdlib'},
            'scope':'only this fixed point/c/weights/balls; local continuous negative density on an open set; no g=edge or numerical cost assumption; no family or C_n bound',
            'seconds':time.monotonic()-started}

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--project',type=Path,default=HERE.parent)
    parser.add_argument('--save',action='store_true')
    args=parser.parse_args();result=run(args.project)
    result['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.save:(HERE/'028_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
