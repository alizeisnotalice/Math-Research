#!/usr/bin/env python3
"""Portable verification of three exact Green point anchors; no search.
Reuses the completed independent all-cone/Taylor audit algorithm.
Requires python-flint and SymPy; the main theorem checker needs neither.
"""
from fractions import Fraction as Q
from collections import Counter
from itertools import permutations,product,combinations
from math import comb,isqrt,prod
from pathlib import Path
from flint import fmpq_poly as Poly,fmpq
import sympy as sp
import json,time
HERE=Path(__file__).resolve().parent

def qp(p):return [Q(str(c)) for c in p]

def lin(b):return Poly([fmpq(str(b)),1])

def poly_product(polys):
    out=Poly([1])
    for p in polys:out*=p
    return out

def determinant(rows):
    a=[[Q(c) for c in row] for row in rows];out=Q(1)
    for k in range(len(a)):
        p=next(j for j in range(k,len(a)) if a[j][k])
        if p!=k:a[k],a[p]=a[p],a[k];out=-out
        out*=a[k][k]
        for j in range(k+1,len(a)):
            ratio=a[j][k]/a[k][k]
            for l in range(k,len(a)):a[j][l]-=ratio*a[k][l]
    return out

def dot(v,w):return sum(a*b for a,b in zip(v,w))

def build_cones():
    body=list(product((1,),(-1,1),(-1,1),(-1,1)))
    rays=[(1,0,0,0),(1,1,0,0),(1,1,1,0),(1,1,1,1),(2,1,1,0),(3,1,1,1)]
    ids=[(0,1,4,5),(4,1,2,3),(4,1,3,5)]
    # Explicit inverse rows, checked against each basis column matrix.
    inverses=[[(1,-1,-1,-1),(0,1,-1,0),(0,0,1,-1),(0,0,0,1)],
              [(1,-1,0,0),(0,1,-1,0),(-1,1,1,-1),(0,0,0,1)],
              [(0,0,1,-1),(0,1,-1,0),(Q(-1,2),Q(1,2),Q(1,2),Q(1,2)),(Q(1,2),Q(-1,2),Q(-1,2),Q(1,2))]]
    base=[]
    for inds,inverse,expected_det in zip(ids,inverses,[1,1,2]):
        cone=tuple(rays[j] for j in inds)
        assert abs(determinant(cone))==expected_det
        assert [[dot(row,col) for col in cone] for row in inverse]==[[int(i==j) for j in range(4)] for i in range(4)]
        base.append((cone,expected_det))
    cones=[]
    # Act by sending old coordinate j to position permutation[j]. This is
    # an independently chosen action convention; the entire group is used.
    for perm in permutations(range(4)):
        for signs in product((-1,1),repeat=4):
            for cone,det in base:
                transformed=[]
                for ray in cone:
                    v=[0]*4
                    for j in range(4):v[perm[j]]=signs[j]*ray[j]
                    transformed.append(tuple(v))
                cones.append((tuple(transformed),det))
    assert len(cones)==1152 and len({tuple(sorted(c)) for c,d in cones})==1152
    allrays={r for c,d in cones for r in c};assert len(allrays)==240
    rate={r:sum(abs(dot(v,r)) for v in body) for r in allrays}
    assert [rate[r] for r in rays]==[8,8,12,12,16,24]
    for cone,det in cones:
        for v in body:
            signs=[dot(v,r) for r in cone]
            assert min(signs)>=0 or max(signs)<=0
    raw_weight=sum((Q(det,16*prod(rate[r] for r in cone)) for cone,det in cones),Q(0))
    assert raw_weight==Q(11,3072)
    return cones,rate,raw_weight

def density(point,cones,rate):
    raw=[];common_powers=Counter();zero_cones=0;repeated_cones=0
    for cone,det in cones:
        bs=[Q(dot(point,r),rate[r]) for r in cone]
        poles=Counter(b*b for b in bs)
        for pole,m in poles.items():common_powers[pole]=max(common_powers[pole],m)
        zero_cones+=int(Q(0) in poles)
        repeated_cones+=int(len(set(bs))<4)
        e1=sum(bs,Q(0));e2=sum((x*y for x,y in combinations(bs,2)),Q(0))
        e3=sum((prod(triple) for triple in combinations(bs,3)),Q(0));e4=prod(bs)
        real=Poly([fmpq(str(e4)),fmpq(str(-e2)),1])
        imaginary_over_t=Poly([fmpq(str(-e3)),fmpq(str(e1))])
        raw.append((poles,Q(det,16*prod(rate[r] for r in cone)),real,imaginary_over_t))
    H=poly_product([lin(pole)**m for pole,m in common_powers.items()])
    numerator=Poly([0]);imaginary=Poly([0])
    for poles,weight,real,im in raw:
        rawden=poly_product([lin(pole)**m for pole,m in poles.items()])
        multiplier,remainder=divmod(H,rawden);assert not remainder
        numerator+=fmpq(str(weight))*real*multiplier
        imaginary+=fmpq(str(weight))*im*multiplier
    assert not imaginary
    # First cancellation occurs here, after all 1152 terms were summed.
    common=numerator.gcd(H);N=numerator//common;H=H//common
    assert H[0]!=0 and N[0]!=0
    # Normalize denominator leading coefficient to 1 for stable recording.
    lead=H.leading_coefficient();N=N/lead;H=H/lead
    multiplicities={};remaining=H
    for pole in sorted(common_powers):
        m=0
        while remaining.degree()>0:
            quotient,remainder=divmod(remaining,lin(pole))
            if remainder:break
            m+=1;remaining=quotient
        if m:multiplicities[pole]=m
    assert remaining.degree()==0 and Q(0) not in multiplicities
    assert H.degree()-N.degree()==2
    leading=Q(str(N.leading_coefficient()/H.leading_coefficient()))
    assert leading==Q(11,3072)
    return N,H,multiplicities,{'raw_common_denominator_degree':sum(common_powers.values()),
      'raw_zero_pole_multiplicity':common_powers.get(Q(0),0),
      'cones_with_zero_projection':zero_cones,'cones_with_repeated_signed_projection':repeated_cones,
      'complete_imaginary_numerator_zero':True,'reduced_numerator_degree_u':N.degree(),
      'reduced_denominator_degree_u':H.degree(),'t0_pi_minus4':str(Q(str(N[0]/H[0]))),
      'leading_t_minus4_coefficient_pi_minus4':str(leading),
      'N_coefficients_ascending_u':list(map(str,qp(N))),'D_coefficients_ascending_u':list(map(str,qp(H))),
      'pole_squared_multiplicities':{str(p):m for p,m in multiplicities.items()}}

def shift_coefficients(poly,root,count):
    # Taylor coefficients of poly(root+v), via the binomial theorem.
    coeffs=qp(poly)
    return [sum((coeffs[j]*comb(j,k)*root**(j-k) for j in range(k,len(coeffs))),Q(0)) for k in range(count)]

def integrate_rational(N,H,poles):
    assert N.degree()<H.degree()
    reconstructed=Poly([0]);coefficient=Q(0);terms=[]
    for squared,m in sorted(poles.items()):
        assert squared>0
        sn,sd=isqrt(squared.numerator),isqrt(squared.denominator)
        assert sn*sn==squared.numerator and sd*sd==squared.denominator
        a=Q(sn,sd)
        other,remainder=divmod(H,lin(squared)**m);assert not remainder
        numerator=shift_coefficients(N,-squared,m)
        denominator=shift_coefficients(other,-squared,m);assert denominator[0]
        expansion=[]
        for j in range(m):
            expansion.append((numerator[j]-sum((denominator[k]*expansion[j-k] for k in range(1,j+1)),Q(0)))/denominator[0])
        for j,c in enumerate(expansion):
            k=m-j
            if not c:continue
            termden=lin(squared)**k
            multiplier,remainder=divmod(H,termden);assert not remainder
            reconstructed+=fmpq(str(c))*multiplier
            coefficient+=c*Q(comb(2*k-2,k-1),2*4**(k-1))/a**(2*k-1)
            terms.append({'a':str(a),'k':k,'coefficient':str(c)})
    assert reconstructed==N
    return coefficient,terms

def parse_author_density(text):
    t=sp.Symbol('t');expr=sp.sympify(text,locals={'t':t})
    N,H=map(lambda x:sp.Poly(x,t),sp.fraction(expr))
    assert all(e[0]%2==0 for e,c in N.terms()+H.terms())
    def convert(p):
        out=[Q(0)]*(p.degree()//2+1)
        for (power,),c in p.terms():out[power//2]=Q(c)
        return Poly([fmpq(str(c)) for c in out])
    return convert(N),convert(H)

def main():
    start=time.perf_counter();saved=json.loads((HERE/'green_anchors.json').read_text())
    cones,rate,peak=build_cones();forms={}
    for name,x in [('A',(1,1,1,0)),('B',(1,1,0,0))]:
        N,H,poles,metadata=density(x,cones,rate)
        savedN,savedH=parse_author_density(saved['p4_pi_minus4'][name]);assert N*savedH==H*savedN
        forms[name]=(N,H,poles)
    answers={}
    targets={'AA':Q(10495703317,2149908480000),'AB':Q(15547486331,1920360960000),'BB':Q(143,10584)}
    for label in targets:
        (N,H,poles),(n,h,other)=[forms[c] for c in label]
        merged=Counter(poles);merged.update(other)
        value,pf=integrate_rational(N*n,H*h,merged);assert value==targets[label]
        answers[label]=str(value)
    assert peak==Q(11,3072) and peak*peak==Q(121,9437184)
    print(json.dumps({'status':'PASS','scope':'exact point values and global peak only; no g_min, I or new C8','cone_count':len(cones),'G8_pi_minus7':answers,'heat_peak_C_pi_minus8':str(peak*peak),'wall_seconds':time.perf_counter()-start},indent=2))
if __name__=='__main__':main()
