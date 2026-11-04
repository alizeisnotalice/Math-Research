"""Exact, radius-continuous optimization for full 2D square holes at stage S1."""

from fractions import Fraction as F
from itertools import product


def best(u, p, q, x):
    w = u / (1 + 2*u)
    a = 4*u*u
    cutoff = F(6, 5)
    # For these inputs ||T k||_infty >= (1-max(p,q))||k||_infty >= 2||k||/3.
    # Thus the range of indices contains every hole intersecting r<=6/5.
    centres = [(i+p*j,j+q*i) for i,j in product(range(-3,4), repeat=2)]
    centres = [z for z in centres if all(abs(zi-xi) <= cutoff+w for zi,xi in zip(z,x))]
    knots = {F(0),cutoff}
    for z in centres:
        for zi, xi in zip(z,x):
            for endpoint in (zi-w,zi+w):
                d = abs(endpoint-xi)
                if d <= cutoff:
                    knots.add(d)
    knots = sorted(knots)

    def area(r):
        result = F(0)
        for z in centres:
            lengths = [max(F(0), min(xi+r,zi+w)-max(xi-r,zi-w)) for zi,xi in zip(z,x)]
            result += lengths[0]*lengths[1]
        return result

    options = set(knots[1:])
    for left,right in zip(knots,knots[1:]):
        mid = (left+right)/2
        vol_left, vol_mid, vol_right = area(left),area(mid),area(right)
        sl = (vol_mid-vol_left)/(mid-left)
        sr = (vol_right-vol_mid)/(right-mid)
        aa = (sr-sl)/(right-left)
        bb = sl-aa*(mid+left)
        cc = vol_left-aa*left*left-bb*left
        if bb:
            stationary = -2*cc/bb
            if left < stationary < right:
                options.add(stationary)
    r = min(options, key=lambda s: (area(s)/(a*s*s),s))
    minimum = area(r)/(a*r*r)
    # On r>=6/5 use full rows.  The deviation from (i+p*j,j) is (0,q*i).
    # |i|<=3m/2 yields at most 6u^3 m of extra range. Section S2 then gives
    # 2m(2m+1)g^2/r^2 >= (3/2)g^2/(1+6u^3+w)^2 for every m>=1.
    assert cutoff >= 1+6*u**3+2*w
    g = 1/(1+2*u)
    assert minimum < F(3,2)*g*g/(1+6*u**3+w)**2
    return r, minimum


def run():
    for u in (F(1,32),F(1,64),F(1,128)):
        g = 1/(1+2*u)
        w = u*g
        q = 4*u**3
        for tag,p,x1 in (("macro",F(1,3),F(0)),("intermediate",u*u,-w+u*u/8)):
            for h,expected in ((q/64,False),(q/8,tag == "intermediate"),(2*q,True)):
                r,ratio = best(u,p,q,(x1,w-h))
                assert (ratio<1) == expected
                if tag == "macro":
                    assert r == g+h-q
                print(f"u={u} type={tag:<12} depth/q={h/q}  r={r} min(deficit/threshold)={ratio} captured={expected}")


if __name__ == "__main__":
    run()
