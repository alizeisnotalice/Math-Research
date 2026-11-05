#!/usr/bin/env python3
"""Directed-decimal certificate for the r002 finite staircase witnesses.

Source JSON decimals are exact inputs. Each arithmetic operation is rounded
outwards at 90 digits. log and exp use Taylor series with explicit tails.
"""
from decimal import Decimal, Context, ROUND_FLOOR, ROUND_CEILING, localcontext
from fractions import Fraction
from math import comb, factorial
from pathlib import Path
import json
import hashlib

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'rounds/r002/multilayer.json'
OUT=ROOT/'rounds/r004/certificates.json'
PREC=90; TERMS=60
DOWN=Context(prec=PREC,rounding=ROUND_FLOOR); UP=Context(prec=PREC,rounding=ROUND_CEILING)
ZERO=Decimal(0); ONE=Decimal(1); TWO=Decimal(2)
def op(fun,rounding):
    with localcontext(DOWN if rounding==ROUND_FLOOR else UP): return +fun()
def lo(fun): return op(fun,ROUND_FLOOR)
def hi(fun): return op(fun,ROUND_CEILING)
def iv(a,b): assert a<=b; return a,b
def add(x,y): return iv(lo(lambda:x[0]+y[0]),hi(lambda:x[1]+y[1]))
def sub(x,y): return iv(lo(lambda:x[0]-y[1]),hi(lambda:x[1]-y[0]))
def mul(x,y):
    vals=[lo(lambda a=a,b=b:a*b) for a in x for b in y]+[hi(lambda a=a,b=b:a*b) for a in x for b in y]
    return iv(min(vals),max(vals))
def scale(x,a): return iv(lo(lambda:x[0]*a),hi(lambda:x[1]*a)) if a>=0 else iv(lo(lambda:x[1]*a),hi(lambda:x[0]*a))
def powi(x,n): assert x[0]>=0; return iv(lo(lambda:x[0]**n),hi(lambda:x[1]**n))
def divi(x,y): assert x[0]>=0 and y[0]>0; return iv(lo(lambda:x[0]/y[1]),hi(lambda:x[1]/y[0]))
def neg(x): return x.copy_negate()
def decimal_fraction(x): return Fraction(x)
def fraction_i(q):
    return iv(lo(lambda:Decimal(q.numerator)/Decimal(q.denominator)),
              hi(lambda:Decimal(q.numerator)/Decimal(q.denominator)))

def atanh(z):
    """2 atanh(z); tail <= 2z^(2N+3)/((2N+3)(1-z^2))."""
    total=ZERO; power=z
    for k in range(TERMS+1):
        total=lo(lambda total=total,power=power,k=k:total+power/Decimal(2*k+1)); power=lo(lambda power=power:power*z*z)
    lower=lo(lambda:TWO*total)
    total=ZERO; power=z
    for k in range(TERMS+1):
        total=hi(lambda total=total,power=power,k=k:total+power/Decimal(2*k+1)); power=hi(lambda power=power:power*z*z)
    zsq=hi(lambda:z*z)
    denominator=lo(lambda:Decimal(2*TERMS+3)*(ONE-zsq))
    upper=hi(lambda:TWO*(total+power/denominator))
    return iv(lower,upper)
def atanh_i(z):
    """Monotonicity turns point enclosures into an enclosure on [z0,z1]."""
    a=atanh(z[0]); b=atanh(z[1]); return iv(a[0],b[1])
LN2=atanh_i(iv(lo(lambda:ONE/Decimal(3)),hi(lambda:ONE/Decimal(3))))
def ln(x):
    assert x>0
    xf=decimal_fraction(x); p=Fraction(1); e=0
    while xf>=2*p: p*=2; e+=1
    while xf<p: p/=2; e-=1
    m=divi(iv(x,x),fraction_i(p))
    # The exact mantissa is in [1,2), so z is nonnegative.  Use interval
    # subtraction/addition/division for the upper endpoint; direct formulas
    # would round the correlated denominator in the wrong direction.
    numerator=sub(m,iv(ONE,ONE))
    numerator=iv(max(ZERO,numerator[0]),numerator[1])
    z=divi(numerator,add(m,iv(ONE,ONE)))
    return add(atanh_i(z),scale(LN2,Decimal(e)))
def exp_pos(u):
    """u in [0,1/2]: Taylor lower and positive geometric upper tail."""
    total=ONE; term=ONE
    for k in range(1,TERMS+1):
        term=lo(lambda term=term,k=k:term*u/Decimal(k)); total=lo(lambda total=total,term=term:total+term)
    lower=total; total=ONE; term=ONE
    for k in range(1,TERMS+1):
        term=hi(lambda term=term,k=k:term*u/Decimal(k)); total=hi(lambda total=total,term=term:total+term)
    first=hi(lambda:term*u/Decimal(TERMS+1))
    quotient=hi(lambda:u/Decimal(TERMS+2))
    denominator=lo(lambda:ONE-quotient)
    return iv(lower,hi(lambda:total+first/denominator))
def exp(x):
    if x==ZERO:return iv(ONE,ONE)
    if x<ZERO:
        a=exp(neg(x)); return iv(lo(lambda:ONE/a[1]),hi(lambda:ONE/a[0]))
    xf=decimal_fraction(x); p=Fraction(1); s=0
    while xf>Fraction(1,2)*p: p*=2; s+=1
    u=divi(iv(x,x),fraction_i(p))
    a0=exp_pos(u[0]); a1=exp_pos(u[1]); a=iv(a0[0],a1[1])
    for _ in range(s): a=mul(a,a)
    return a

def tail(n,r,t,lnr):
    """Exact finite-mixture identity, with directed elementary enclosures."""
    den=powi(add(iv(ONE,ONE),iv(r,r)),n); emt=exp(neg(t)); ans=iv(ZERO,ZERO)
    for k in range(n+1):
        x=add(iv(t,t),scale(lnr,Decimal(k))); c=Decimal(comb(n,k))
        if x[1]<=ZERO: ans=add(ans,iv(lo(lambda:c*r**k),hi(lambda:c*r**k)))
        else:
            if x[0]<=ZERO: raise ArithmeticError('x crosses zero')
            poly=iv(ZERO,ZERO)
            for j in range(n):
                reciprocal=divi(iv(ONE,ONE),iv(Decimal(factorial(j)),Decimal(factorial(j))))
                poly=add(poly,mul(powi(x,j),reciprocal))
            ans=add(ans,scale(mul(emt,poly),c))
    return divi(ans,den)
def phi(n,h):
    if h==ZERO:return iv(ZERO,ZERO)
    if h<=ONE:return iv(h,h)
    return scale(powi(add(iv(ONE,ONE),ln(h)),n-1),h)
def fmt(x): return str(x)
def certify(row):
    n=int(row['n']); r=row['r']; ts=row['thresholds']; hs=row['heights']
    assert n>=2 and ZERO<r<=ONE and len(ts)==len(hs)>0
    assert ts[0]>=ZERO and all(a<b for a,b in zip(ts,ts[1:]))
    assert hs[0]>=ZERO and all(a<=b for a,b in zip(hs,hs[1:]))
    lnr=ln(r)
    tails=[tail(n,r,t,lnr) for t in ts]+[iv(ZERO,ZERO)]
    ps=[sub(tails[j],tails[j+1]) for j in range(len(ts))]
    ex=iv(ZERO,ZERO); cost=iv(ZERO,ZERO)
    for j,(t,h,p) in enumerate(zip(ts,hs,ps)):
        ex=add(ex,scale(p,h)); w=exp(neg(t)) if j+1==len(ts) else sub(exp(neg(t)),exp(neg(ts[j+1])))
        cost=add(cost,mul(w,phi(n,h)))
    assert cost[0]>ZERO and all(p[0]>=ZERO for p in ps)
    logscore=scale(ln(cost[1]),Decimal(-1))
    return {'n':n,'expectation_interval':[fmt(ex[0]),fmt(ex[1])],
      'strict_excess_lower':fmt(lo(lambda:ex[0]-ONE)), 'cost_interval':[fmt(cost[0]),fmt(cost[1])],
      'log_score_lower':fmt(logscore[0]), 'certified':ex[0]>ONE}
def main():
    with SOURCE.open() as f:data=json.load(f,parse_float=Decimal)
    rows=[certify(x) for x in data['rows']]
    result={'method':'exact decimal inputs; every +,-,*,/,power outward-rounded at 90 digits; log=2atanh((m-1)/(m+1)) with 60 terms and tail; exp has 60 Taylor terms, a positive geometric tail, and outward repeated squaring',
      'source':'rounds/r002/multilayer.json','source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'precision_digits':PREC,'terms':TERMS,'rows':rows,'all_certified':all(x['certified'] for x in rows)}
    OUT.write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result,indent=2))
    assert result['all_certified'], 'A witness constraint was not certified'
if __name__=='__main__':main()
