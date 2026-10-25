"""One fixed exact finite-epsilon certificate; no parameter search."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import hashlib,json
import sympy as sy
from vg6_bounds import Endpoint as E,geometry_constants,geometry_from_caps,period_constants,pure_constants,certify_gain

P=Path(__file__).resolve().parent


def encode(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [encode(v) for v in x]
    return x


def volume(boxes):
    ans=F(0)
    for k in range(1,len(boxes)+1):
        for subset in combinations(boxes,k):
            v=F(1)
            for j in range(4):
                v*=max(F(0),min(b[1][j] for b in subset)-max(b[0][j] for b in subset))
            ans+=(-1)**(k+1)*v
    return ans


def main():
    boxes=[([E(-2,0),E(0,0),E(0,0),E(0,0)],
            [E(0,-1),E(1,1,1),E(1,2),E(1,3)]),
           ([E(0,2,1),E(0,0),E(0,0),E(0,0)],
            [E(2,0),E(1,-1),E(1,1),E(1,-2)]),
           ([E(0,F(-1,2)),E(0,0),E(0,0),E(0,0)],
            [E(0,F(1,2)),E(1,0),E(1,0),E(1,0)])]
    geom=geometry_constants(boxes,F(1,2))
    caps=geometry_from_caps(geom['ranges'],geom['velocity_variation'],geom['max_speeds'],
        geom['max_remainders'],[F(1),F(1),F(0),F(0)],geom['min_distinct_gaps'],
        geom['endpoint_speed_sums'],geom['tmax'])
    assert caps==geom
    period=period_constants([0,0,0,0],geom['eta'])
    drift=[F(1,3),F(-1,4),F(1,5),F(-1,6)]
    pure=pure_constants(geom,[1,1,1,1],list(map(abs,drift)),[2,1,1,1],1)
    V0=F(4);D=F(6)+F(15,16);g=D-V0
    cert=certify_gain(V0,D,g,geom['C']+pure['C_period'],pure['K'],F(1,2),m1_abs=1)
    eps=cert['epsilon'];s=cert['fourth_root'];r=s/2
    actual=[]
    for lo,hi in boxes:
        actual.append((tuple(e.base+eps*e.velocity+eps**2*e.remainder for e in lo),
                       tuple(e.base+eps*e.velocity+eps**2*e.remainder for e in hi)))
    center=[F(-2)+drift[0]*eps]+[drift[j]*eps for j in range(1,4)]
    satellite=(tuple(c-r for c in center),tuple(c+r for c in center))
    main_volume=volume(actual)
    full_volume=volume(actual+[satellite])
    x=sy.symbols('x')
    polynomial=sy.expand((2-x)*(1+x+x*x)*(1+2*x)*(1+3*x)
                         +(2-2*x-x*x)*(1-x)*(1+x)*(1-2*x)+x)
    assert sy.diff(polynomial,x).subs(x,0)==6
    exact_poly=F(int(sy.numer(polynomial.subs(x,sy.Rational(eps.numerator,eps.denominator)))),
                 int(sy.denom(polynomial.subs(x,sy.Rational(eps.numerator,eps.denominator)))))
    assert main_volume==exact_poly
    overlap=F(1)
    for c in drift:overlap*=r+c*eps
    assert full_volume-main_volume==eps-overlap
    error=abs(full_volume-V0-D*eps)
    bound=geom['C']*eps**2+pure['K']*s**7
    assert error<=bound
    actual_improvement=full_volume/(1+eps)-V0
    assert actual_improvement>=cert['ratio_improvement_lower']>0
    out=dict(status='one exact finite-family certificate PASS; not a maximal-operator input',
             geometry=geom,period=period,pure=pure,V0=V0,D=D,g_lower=g,certificate=cert,
             main_volume_polynomial=str(polynomial),main_volume=main_volume,
             pure_added_volume=full_volume-main_volume,full_volume=full_volume,
             actual_error=error,error_bound=bound,actual_ratio_improvement=actual_improvement,
             source_sha256={n:hashlib.sha256((P/n).read_bytes()).hexdigest()
                            for n in ['vg6_bounds.py','vg6_toy.py']})
    (P/'vg6_results.json').write_text(json.dumps(encode(out),indent=2)+'\n')
    print(json.dumps(encode(out),indent=2))


if __name__=='__main__':main()
