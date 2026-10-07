#!/usr/bin/env python3
"""Exact rational consistency certificates, not a numerical proof in all dimensions."""
from fractions import Fraction as F
import json
ATOMS=[(F(0),F(3,20)),(F(11,40),F(3,40)),(F(81,400),F(3,400)),(F(10),F(307,400))]
X,Z=F(1,5),-F(3,50)
def radius(j):return F(4,2**j)
def trace(t):return {j:sum((w for a,w in ATOMS if abs(t-a)<radius(j)),F(0))/(2*radius(j)) for j in range(1,11)}
def maximal(t):
    ds=sorted(set(abs(t-a) for a,w in ATOMS))
    return max(sum((w for a,w in ATOMS if abs(t-a)<=d),F(0))/(2*d) for d in ds)
def main():
    assert sum(w for a,w in ATOMS)==1
    tx,tz=trace(X),trace(Z)
    assert maximal(X)==F(3,2) and maximal(Z)==F(5,4)
    assert next(j for j,g in tx.items() if g>F(1,8))==3
    assert next(j for j,g in tz.items() if g>F(1,8))==3
    assert next(j for j,g in tx.items() if g>F(1,2))==10
    assert next(j for j,g in tz.items() if g>F(1,2))==5
    W=min(F(1,2),tx[4],tz[5])-max(F(1,4),max(tx[j] for j in range(1,4)),max(tz[j] for j in range(1,5)))
    lost=sum((w for a,w in ATOMS if abs(X-a)<radius(4) and abs(Z-a)>radius(4)),F(0))
    assert W==F(33,200) and lost==F(33,400) and W*2*radius(4)/lost==1
    smooth=[]
    for batch,eps in enumerate((F(1,10**5),F(1,10**6),F(1,10**7))):
        h=eps/10
        # Interval certificate: every atom may move independently within h.
        for t in (X,Z):
            margin=min(abs(abs(t-a)-radius(j)) for a,w in ATOMS for j in range(1,11))
            assert margin>eps+h
            dmin=min(abs(t-a) for a,w in ATOMS)
            q=(eps+h)/dmin
            # For every observer and every smooth source displacement,
            # all cumulative mass/distance ratios obey these multiplicative bounds.
            assert maximal(t)/(1+q)>1 and maximal(t)/(1-q)<2
            for v in (t-eps,t+eps):assert trace(v)==trace(t)
        assert X-Z-2*eps>radius(4)
        assert abs(X)+eps+h<radius(4) and abs(Z)+eps+h<radius(5)
        fee=W/F(1,4)*F(3,20)/(2*radius(4)*2*radius(5))*4*eps**2
        assert fee==F(396,125)*eps**2>0
        smooth.append(dict(batch=batch,epsilon=str(eps),smoothing_support=str(h),theta='1',positive_fee=str(fee)))
    assert F(7,8)**8<F(1,2) # implies 2^(-2/n)>7/8 for n>=16
    annuli=[]
    for batch,ns in enumerate(((512,1024,2048),(4096,8192,16384),(32768,65536,131072))):
        for n in ns:
            A=F(n,2);lower=(1-F(256,n))*F(n*(n-1),8)
            ratio=lower/A**2
            assert 0<ratio<F(1,2)
            annuli.append(dict(batch=batch,n=n,A=str(A),outer_lower=str(lower),ratio_lower=str(ratio)))
    return dict(status='passed',smooth_cases=smooth,annuli_cases=annuli,theta='1',W=str(W),lost=str(lost),dimension_free_bound_proved=False,scope='Fraction interval inequalities; annuli are arbitrary masks, smooth example is a local entropy obstruction only.')
if __name__=='__main__':print(json.dumps(main(),indent=2))
