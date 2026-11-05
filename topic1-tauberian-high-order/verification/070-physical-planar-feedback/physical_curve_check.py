"""Exact low-order check for the PHYSICAL inverse-capacity curve recursion.

Adapted, with scope changed, from the Fraction curve bookkeeping in
r069-computation/formal/curve_formal.py (read-only).
Here kappa=1, rho=g=(1+2u)^(-1), delta=u(1+2u)/8,
eta=a*delta**2 with a=4u**2.  Omitting the Lambda*u^6 correction to w is
valid for height coefficients through u^6 only.  This script is evidence
for finite algebra, not a proof of the all-layer geometric reduction.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import argparse, json, time

def clean(p): return {k:v for k,v in p.items() if v}
def add(*ps):
    out={}
    for p in ps:
        for k,v in p.items(): out[k]=out.get(k,F(0))+v
    return clean(out)
def sc(p,c): return clean({k:v*c for k,v in p.items()})
def mul(p,q,R):
    out={}
    for (i,j),x in p.items():
        for (k,l),y in q.items():
            if j+l<=R: out[(i+k,j+l)]=out.get((i+k,j+l),F(0))+x*y
    return clean(out)
def diff(p): return {(i-1,j):i*v for (i,j),v in p.items() if i}
def integ(p): return {(i+1,j):v/F(i+1) for (i,j),v in p.items()}
def at1(p):
    out={}
    for (i,j),v in p.items(): out[(0,j)]=out.get((0,j),F(0))+v
    return clean(out)
def at0(p): return {k:v for k,v in p.items() if k[0]==0}
def rev(p):
    out={}
    for (i,j),v in p.items():
        for k in range(i+1):
            out[(k,j)]=out.get((k,j),F(0))+v*comb(i,k)*(-1)**k
    return clean(out)
def const(x): return {(0,0):F(x)} if x else {}
def seq(p,R): return [str(p.get((0,j),F(0))) for j in range(R+1)]

def run(R,depth,nonlinear):
    assert R<=6, "The rho=g low-order specialization is only certified through height order 6."
    delta={(0,1):F(1,8),(0,2):F(1,4)}
    aa={(0,2):F(4)}
    eta=mul(aa,mul(delta,delta,R),R) if nonlinear else {}
    def psi(z):
        # Catalan expansion: delta*sum (-1)^m C_m eta^m z^(m+1).
        term=z; out={}
        for m in range(R+1):
            if not term: break
            cm=F(comb(2*m,m),m+1)*(-1)**m
            out=add(out,sc(mul(delta,term,R),cm))
            term=mul(eta,mul(term,z,R),R)
        return out
    r={(1,0):F(1)}; h=psi(r)
    pieces=[{'r':r,'h':h,'P':integ(h)}]
    rows=[]
    for j in range(depth+1):
        total=at1(pieces[-1]['P'])
        rows.append({'layer':j,'pieces':len(pieces),
                     'A_integral_h_dq_coefficients':seq(total,R),
                     'twice_A_one_corner_area_coefficients':seq(sc(total,2),R),
                     'terms':sum(len(p[x]) for p in pieces for x in ('r','h','P'))})
        if j==depth: break
        lows=[]; highs=[]
        for p in pieces:
            r,h,P=(p[x] for x in ('r','h','P'))
            q=h
            mass=add(sc(mul(q,q,R),F(1,2)),
                     mul(q,add(const(1),sc(r,-1)),R),P)
            lows.append({'r':q,'h':psi(add(q,sc(mass,-1)))})
            q=add(const(1),sc(r,-1),h)
            mass=add(sc(total,2),sc(P,-1),sc(mul(h,h,R),F(1,2)))
            highs.append({'r':rev(q),'h':rev(psi(add(q,sc(mass,-1))))})
        new=lows+list(reversed(highs))
        prefix={}
        for k,p in enumerate(new):
            if k:
                assert at1(new[k-1]['r'])==at0(p['r'])
                assert at1(new[k-1]['h'])==at0(p['h'])
            p['P']=add(prefix,integ(mul(p['h'],diff(p['r']),R)))
            prefix=at1(p['P'])
        assert at0(new[0]['r'])=={}
        assert at1(new[-1]['r'])==const(1)
        pieces=new
    return rows

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--order',type=int,default=6)
    ap.add_argument('--depth',type=int,default=4)
    args=ap.parse_args(); start=time.monotonic()
    physical=run(args.order,args.depth,True)
    plane=run(args.order,args.depth,False)
    differences=[]
    for p,l in zip(physical,plane):
        diffc=[F(x)-F(y) for x,y in zip(p['A_integral_h_dq_coefficients'],
                                     l['A_integral_h_dq_coefficients'])]
        # For this actual parameter path, Aphys-Aplane=-u^5/384+O(u^6).
        assert all(diffc[k]==0 for k in range(min(5,len(diffc))))
        if args.order>=5: assert diffc[5]==-F(1,384)
        differences.append({'layer':p['layer'],'A_physical_minus_plane':[str(x) for x in diffc]})
    result={'scope':'kappa=1, rho=g; actual physical inverse, finite Taylor check only',
            'order':args.order,'depth':args.depth,'elapsed_seconds':time.monotonic()-start,
            'physical':physical,'plane':plane,'differences':differences,
            'all_assertions_passed':True}
    path=Path(__file__).with_name('physical_curve_results.json')
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__': main()
