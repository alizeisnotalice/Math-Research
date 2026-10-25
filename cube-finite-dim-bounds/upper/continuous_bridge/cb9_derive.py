"""Exact twenty-cone algebra only; no cost or face optimization."""
from pathlib import Path
import sympy as s,json
from itertools import product
OUT=Path(__file__).resolve().parent

def cones():
    cc=[]
    for i in range(3):
        j,k=[a for a in range(3) if a!=i]
        for sign in (-1,1):
            def ray(axis,ss):
                r=[0]*3;r[i]=sign;r[axis]=ss;return tuple(r)
            A,C=ray(j,1),ray(j,-1);B,D=ray(k,1),ray(k,-1)
            cc.extend([(A,B,C),(A,C,D)])
    for e in product((-1,1),repeat=3):cc.append(((e[0],e[1],0),(e[0],0,e[2]),(0,e[1],e[2])))
    return cc

def canonical(r):
    first=next(v for v in r if v)
    return tuple(v*(1 if first>0 else -1) for v in r)
if __name__=='__main__':
    cc=cones();assert len(cc)==20
    dets=[s.det(s.Matrix.hstack(*(s.Matrix(r) for r in c))) for c in cc];assert all(abs(d)==2 for d in dets)
    rays=sorted({r for c in cc for r in c});lines=sorted({canonical(r) for r in rays});assert len(rays)==12 and len(lines)==6
    x=s.symbols('x y z',real=True);u=s.symbols('u',real=True)
    L=[sum(a*b for a,b in zip(r,x)) for r in lines];N=0
    for cone in cc:
        ell=[sum(a*b for a,b in zip(r,x)) for r in cone];ids={lines.index(canonical(r)) for r in cone};assert len(ids)==3
        pair=ell[0]*ell[1]+ell[0]*ell[2]+ell[1]*ell[2]
        N+=(u-pair)*s.prod(u+L[j]**2 for j in range(6) if j not in ids)
    poly=s.Poly(s.expand(N),u);coeff={str(deg[0]):str(s.factor(c)) for deg,c in poly.terms()}
    record={'cones':[[list(r) for r in c] for c in cc],'determinants':list(map(int,dets)),'rays':list(map(list,rays)),'unoriented_lines':list(map(list,lines)),'spectral_unit_ball_volume':str(sum(abs(d) for d in dets)/(6*4**3)),'N_coefficients_in_u':coeff,'formula':'sum_cone Re product_j(a-i ell_j)^(-1) = a*N(a^2)/prod_6(a^2+L_j^2)'}
    (OUT/'cb9_cones.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
