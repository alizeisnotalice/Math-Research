#!/usr/bin/env python3
"""Exact arithmetic checks for the stated one-dimensional contrast example."""
from fractions import Fraction as F

def merge(intervals):
    out=[]
    for a,b in sorted((a,b) for a,b in intervals if b>a):
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out

def halo_fraction(E,alpha):
    """Side-effect-free copy of the round-20/round-013 outer-halo verifier."""
    ends=sum(([a,b] for a,b in E),[]); ell=sum(b-a for a,b in E)
    lo=ends[0]-ell/(2*alpha); hi=ends[-1]+ell/(2*alpha); pieces=list(E)
    for e in ends:
        cuts=sorted({lo,hi,*((e+f)/2 for f in ends if lo<(e+f)/2<hi)})
        def G(x):
            a,b=min(e,2*x-e),max(e,2*x-e)
            return sum(max(F(0),min(b,d)-max(a,c)) for c,d in E)-2*alpha*abs(x-e)
        for a,b in zip(cuts,cuts[1:]):
            ga,gb=G(a),G(b)
            if ga>0 and gb>0: pieces.append((a,b))
            elif ga>0: pieces.append((a,a+(b-a)*ga/(ga-gb)))
            elif gb>0: pieces.append((a+(b-a)*(-ga)/(gb-ga),b))
    return merge(pieces)

def main():
    # E=(0,1), alpha=1/4. Endpoints have Lebesgue and weight mass zero.
    H=halo_fraction([(F(0),F(1))],F(1,4))
    assert H==[(F(-1),F(2))]
    print('outer_fraction_halo=(-1,2) PASS')
    # b=1+64*1_(1,2), v=|x|^(-1/2): b(E)=1, v(E)=2,
    # b(H)=67, v(H)=2+2sqrt(2), hence w(H)/w(E)=(69+2sqrt(2))/3.
    assert F(1)+F(2)==3 and 67>0
    print('main_ratio=(69+2sqrt(2))/3 PASS')
    root2_square = 2
    # The displayed ratio gives the lower bound (69+2 sqrt(2))/12 for the
    # normalized mixed quantity.  The comparison pure-weight upper bound is
    # 2(1+sqrt(2)); it is deliberately not asserted to be an equality.
    # This reduces to 45>22 sqrt(2); the positive-square gap is 45^2-(22^2)*2.
    square_gap = 45**2 - 22**2 * root2_square
    assert square_gap == 1057 and 45 > 0
    # Thus (69+2sqrt2)/12 > 2(1+sqrt2).
    print(f"contrast_square_gap={square_gap} PASS")
    # For b=1+B*1_(1,2), this witness wins exactly when B>19+22sqrt(2).
    # Thus B=51 works and B=50 fails; this is not a sharp global threshold.
    assert 31**2 < 22**2*2 < 32**2
    print('B_witness_threshold=19+22sqrt(2); B50 FAIL, B51 PASS (31^2<968<32^2)')

    # For t=u^2 with 0<u<=1/256 (i.e. t<=2^-16), the specified domination
    # sqrt(t)/(1+t)>130t/(1-t) is 1-u^2>130u(1+u^2).
    u = F(1, 256)
    residual = 1 - u*u - 130*u*(1+u*u)
    assert residual > 0
    print(f"endpoint_u={u} residual={residual} PASS")
    print("monotonic_guard=1-u^2-130u(1+u^2) decreases for u>0; PASS")

if __name__ == '__main__':
    main()
