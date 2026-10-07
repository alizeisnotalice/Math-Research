#!/usr/bin/env python3
"""Actual B rational symbols and principal signs; no finite-R/actual-FIRST test."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'ordered_common_excessivity_exact_guard_20261007_results.json'
ORDER=4
ONE=[F(1)]+[F(0)]*ORDER
ZERO=[F(0)]*(ORDER+1)

def add(a,b): return [x+y for x,y in zip(a,b)]
def scale(a,c): return [c*x for x in a]
def mul(a,b):
    return [sum((a[j]*b[k-j] for j in range(k+1)),F(0)) for k in range(ORDER+1)]
def inverse(a):
    assert a[0]
    out=[1/a[0]]
    for k in range(1,ORDER+1):
        out.append(-sum((a[j]*out[k-j] for j in range(1,k+1)),F(0))/a[0])
    return out
def shifted_ratio(a,b):
    assert a[0]==b[0]==0 and b[1]>0
    # Only principal constant and degree-two coefficient are needed.
    q0=a[1]/b[1]
    return q0,(a[2]-b[2]*q0)/b[1]
def product(parts):
    out=ONE
    for p in parts: out=mul(out,p)
    return out
def hash_fraction(x):
    digest=hashlib.sha256()
    for value in (x.numerator,x.denominator):
        absolute=abs(value)
        data=absolute.to_bytes(max(1,(absolute.bit_length()+7)//8),'big')
        digest.update(b'-' if value<0 else b'+')
        digest.update(len(data).to_bytes(8,'big'))
        digest.update(data)
    return digest.hexdigest()

# Direct inverse of log(1+t)/t, not a hard-coded B coefficient formula.
log_over_t=[F((-1)**k,k+1) for k in range(ORDER+1)]
B=add(inverse(log_over_t),scale(ONE,-1))
assert B[1:3]==[F(1,2),F(-1,12)]
checks_total=1
rows=[]

def b_interval(x):
    if x==0: return F(0),F(0)  # No log(1+0)/0 or B(0)/0.
    assert 0<x<1
    lo=sum(((-1)**(j+1)*x**j/j for j in range(1,13)),F(0))
    hi=lo+x**13/13
    return x/hi-1,x/lo-1

for n in (8,32,128):
    checks=0
    eps=F(1,32*n)
    vectors={'axis':[1]+[0]*(n-1), 'balanced':[1]*n,
             'mixed':[1+(j%3) for j in range(n)]}
    signatures=[]
    representatives=[]
    for c in (F(1,2),F(3,4),F(1)):
      for r in (F(0),F(1,4),F(1,2),F(3,4),F(1)):
       axis_bracket=F(n-1)
       diagonal_bracket=F(n+2,n)-3
       physical_b=-c*c*r*(1-r/2)*axis_bracket/2
       physical_a=r*r*c*diagonal_bracket/4
       assert (physical_b<0 and physical_a<0) if r else (physical_b==physical_a==0)
       checks+=2
       for name,vec in vectors.items():
        bs=[[B[k]*F(v*v)**k for k in range(ORDER+1)] for v in vec]
        As=[]
        for b in bs:
            cb=scale(b,c)
            As.append(mul(cb,inverse(add(ONE,cb))))
        sum_b=sum_a=ZERO
        for b,a in zip(bs,As):
            sum_b=add(sum_b,b); sum_a=add(sum_a,a)
        T=product([add(ONE,scale(a,-r)) for a in As])
        numerator=add(ONE,scale(T,-1))
        eB=[ONE]+[ZERO for _ in range(ORDER)]
        eA=[ONE]+[ZERO for _ in range(ORDER)]
        for b,a in zip(bs,As):
            for k in range(ORDER,0,-1):
                eB[k]=add(eB[k],mul(eB[k-1],b))
                eA[k]=add(eA[k],mul(eA[k-1],a))
        G=product([inverse(add(ONE,scale(b,c))) for b in bs])
        sumB=sumA=ZERO
        for k in range(1,ORDER+1):
            sumB=add(sumB,scale(eB[k],c**k*(1-(1-r)**k)))
            sumA=add(sumA,scale(eA[k],(-1)**(k+1)*r**k))
        rhsB=mul(G,sumB)
        for left,right in zip(numerator,rhsB): assert left==right; checks+=1
        for left,right in zip(numerator,sumA): assert left==right; checks+=1
        s2=sum(F(v*v) for v in vec)
        s4=sum(F(v**4) for v in vec)
        cross=(s2*s2-s4)/2
        expected_b=-c*c*r*(1-r/2)*s4/(2*s2)-r*r*c*c*s2/4
        expected_a=-r*r*c*cross/(2*s2)
        ratio_b=shifted_ratio(numerator,sum_b)
        ratio_a=shifted_ratio(numerator,sum_a)
        assert ratio_b==(r*c,expected_b)
        assert ratio_a==(r,expected_a)
        checks+=4
        # Certified finite Fourier points for actual B, using rational log enclosures.
        lob=hib=loa=hia=F(0)
        product_lo=product_hi=F(1)
        for v in vec:
            blo,bhi=b_interval(eps*eps*v*v)
            assert 0<=blo<=bhi
            alo=c*blo/(1+c*blo); ahi=c*bhi/(1+c*bhi)
            assert 0<=alo<=ahi<1
            lob+=blo; hib+=bhi; loa+=alo; hia+=ahi
            product_lo*=1-r*ahi; product_hi*=1-r*alo
            checks+=2
        numlo=1-product_hi; numhi=1-product_lo
        rb=(numlo/hib,numhi/lob)
        ra=(numlo/hia,numhi/loa)
        coefficient_b=((rb[0]-r*c)/(eps*eps),(rb[1]-r*c)/(eps*eps))
        coefficient_a=((ra[0]-r)/(eps*eps),(ra[1]-r)/(eps*eps))
        tolerance=10*eps*eps*s2*s2
        for bounds,expected in ((coefficient_b,expected_b),(coefficient_a,expected_a)):
            assert expected-tolerance<bounds[0]<=bounds[1]<expected+tolerance
            checks+=2
            signatures.extend(hash_fraction(x) for x in bounds)
        if c==1 and r==F(1,2) and name=='balanced':
            representatives.append({'c':str(c),'r':str(r),'direction':name,
                                    'DB_degree2':str(expected_b),'DA_degree2':str(expected_a),
                                    'finite_probe_tolerance':str(tolerance)})
    rows.append({'n':n,'epsilon':str(eps),'checks':checks,
                 'actual_symbol_probe_count':45,'formal_case_count':45,
                 'principal_coefficient_representatives':representatives,
                 'enclosure_hashes_sha256':hashlib.sha256(''.join(signatures).encode()).hexdigest(),
                 'status':'PASS_EXACT_ORIGINAL_SYMBOL_PRINCIPAL_COMPONENTS_ONLY'})
    checks_total+=checks

result={'scope':'Original B/Gc low-frequency coefficients, identities and angular principal signs; not finite-R Green inequality, uniform-n remainder, L1 obstacle, actual FIRST/history or weak endpoint',
        'arithmetic':'Fraction formal-series inverse and rational alternating-log enclosures',
        'random_seed':'none; deterministic','checks':checks_total,'rounds':rows,
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
OUT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'checks':checks_total,'rounds':[{'n':row['n'],'checks':row['checks'],'status':row['status']} for row in rows],'output':str(OUT)},ensure_ascii=False))
