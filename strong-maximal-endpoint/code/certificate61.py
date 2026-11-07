#!/usr/bin/env python3
"""Directed-Decimal finite-multiplier certificate for r061.

This verifier is deliberately independent of r060's quadrature evaluator.
It treats every displayed input decimal as an exact rational, constructs a
finite nonnegative multiplier, and proves a lower bound for its Fenchel dual
by rectangular *lower* sums for q(a)=a-Phi^*_[0,2](a).  The source geometry
is the 33-node gamma=1/5 parabola selected in r060.

The only elementary transcendental evaluations are Decimal exp/ln calls at
90 digits, expanded by one representable directed endpoint.  All remaining
operations are directed Decimal interval operations.  This is a certificate
for the stated finite 33 constraints only; it says nothing about an all-v
family or an endpoint constant.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from bisect import bisect_right
from dataclasses import dataclass
from decimal import Decimal, Context, Inexact, ROUND_CEILING, ROUND_FLOOR, localcontext
from pathlib import Path

PREC = 90
DOWN = Context(prec=PREC, rounding=ROUND_FLOOR)
UP = Context(prec=PREC, rounding=ROUND_CEILING)
NEAR = Context(prec=PREC)
EXACT = Context(prec=500)
Z, O = Decimal(0), Decimal(1)
N, B, CAP = 64, 32, Decimal(2)

# Exact archived decimal parameters.  The coefficients are the r060 selected
# gamma=1/5 corrected weights, rounded *down* to these finite decimals; all
# are nonnegative.  A common multiplicative rescaling would also be valid,
# but is unnecessary for the lower dual bound below.
C = Decimal('1.2360667450799148')
TAU = Decimal('66.43897513452868')
GAMMA = Decimal('0.2')
WGROUP = [Decimal(x) for x in (
    '0.9027983090612243', '0.021024307061216656',
    '0.01594485961321416', '0.002507817248475967', '0',
    '0.07361479558522702', '0.13847413030461514',
    '0.18159938252811883', '0.18799750402983434',
    '0.13912249795573456', '0.016745853784823114', '0', '0', '0',
    '0', '0', '0.16150601653160584')]
LO, HI = Decimal('-15.2'), Decimal('160')
FACT = [Decimal(math.factorial(j)) for j in range(B + 1)]
BIN = [Decimal(math.comb(B, k)) for k in range(B + 1)]


def down(f):
    with localcontext(DOWN): return +f()
def up(f):
    with localcontext(UP): return +f()
def exact(f):
    with localcontext(EXACT) as ctx:
        ctx.clear_flags(); ans = +f()
        if ctx.flags[Inexact]:
            raise ArithmeticError('exact rational preprocessing exhausted 500 digits')
        return ans
def outward(x: Decimal, upper: bool) -> Decimal:
    with localcontext(UP if upper else DOWN):
        return x.next_plus() if upper else x.next_minus()


@dataclass(frozen=True)
class IV:
    lo: Decimal
    hi: Decimal
    def __post_init__(self):
        if self.lo > self.hi: raise ValueError((self.lo, self.hi))

def pt(x): return IV(x, x)
def add(x, y): return IV(down(lambda: x.lo+y.lo), up(lambda: x.hi+y.hi))
def sub(x, y): return IV(down(lambda: x.lo-y.hi), up(lambda: x.hi-y.lo))
def mul(x, y):
    ls = [down(lambda a=a,b=b: a*b) for a in (x.lo,x.hi) for b in (y.lo,y.hi)]
    hs = [up(lambda a=a,b=b: a*b) for a in (x.lo,x.hi) for b in (y.lo,y.hi)]
    return IV(min(ls), max(hs))
def div(x, y):
    if y.lo <= Z <= y.hi: raise ZeroDivisionError(y)
    ls = [down(lambda a=a,b=b: a/b) for a in (x.lo,x.hi) for b in (y.lo,y.hi)]
    hs = [up(lambda a=a,b=b: a/b) for a in (x.lo,x.hi) for b in (y.lo,y.hi)]
    return IV(min(ls), max(hs))
def scal(x, a): return mul(x, pt(a))
def powi(x, n):
    if n < 0 or x.lo < Z: raise ValueError((x,n))
    ans, base = pt(O), x
    while n:
        if n & 1: ans = mul(ans, base)
        n >>= 1
        if n: base = mul(base, base)
    return ans
def expiv(x):
    with localcontext(NEAR): lo, hi = x.lo.exp(), x.hi.exp()
    return IV(outward(lo, False), outward(hi, True))
def lniv(x):
    if x.lo <= Z: raise ValueError(x)
    with localcontext(NEAR): lo, hi = x.lo.ln(), x.hi.ln()
    return IV(outward(lo, False), outward(hi, True))
def posdiff(x, y):
    z = sub(x, y)
    return IV(max(Z,z.lo), max(Z,z.hi))


class QUnder:
    """A finite, certified concave chord minorant of q=a-Phi^*_[0,2](a)."""
    def __init__(self, hgrid=512):
        xs, ys = [Z, O, Decimal(N)], [Z, O, O]
        last = Decimal(N)
        for i in range(1, hgrid):
            h = exact(lambda i=i: O + Decimal(i)/Decimal(hgrid))
            L = lniv(pt(h)); oneL = add(pt(O), L)
            phi = mul(pt(h), powi(oneL, N-1))
            apr = mul(powi(oneL, N-2), add(pt(Decimal(N)), L))
            # x is left of the tangent abscissa.  The adverse endpoint in
            # (1-h)a yields a lower q value valid at that node.
            x = apr.lo
            y = add(pt(phi.lo), mul(pt(exact(lambda: O-h)), pt(apr.hi))).lo
            if x <= last: raise ArithmeticError(('q nodes not ordered',i,last,x))
            xs.append(x); ys.append(y); last = x
        L2 = lniv(pt(CAP)); oneL2 = add(pt(O), L2)
        phi2 = mul(pt(CAP), powi(oneL2, N-1))
        acap = mul(powi(oneL2, N-2), add(pt(Decimal(N)), L2))
        if acap.hi <= last: raise ArithmeticError('cap q node order')
        xs.append(acap.hi); ys.append(down(lambda: phi2.lo-acap.hi))
        self.xs, self.ys, self.phi2lo = xs, ys, phi2.lo

    def value(self, a: Decimal) -> Decimal:
        # This returns a lower value.  q is concave; on each tabulated
        # interval its chord lies below q.  The final affine branch is exact
        # up to the directed lower phi(2) endpoint.
        if a <= Z: return Z
        if a <= O: return a
        if a <= Decimal(N): return O
        if a >= self.xs[-1]: return down(lambda: self.phi2lo-a)
        i = bisect_right(self.xs, a)-1
        x0,x1,y0,y1 = self.xs[i],self.xs[i+1],self.ys[i],self.ys[i+1]
        slope = div(pt(exact(lambda: y1-y0)), pt(exact(lambda: x1-x0)))
        return add(pt(y0), mul(slope, pt(exact(lambda: a-x0)))).lo


def geometry():
    vs = list(range(17))
    ds = [exact(lambda v=v: GAMMA*Decimal(v*v)/Decimal(N)) for v in vs]
    ux = [Z] + [exact(lambda v=v: ds[v]+Decimal(v)) for v in range(1,17)] + [exact(lambda v=v: ds[v]-Decimal(v)) for v in range(1,17)]
    uy = [Z] + [exact(lambda v=v: ds[v]-Decimal(v)) for v in range(1,17)] + [exact(lambda v=v: ds[v]+Decimal(v)) for v in range(1,17)]
    w = WGROUP + WGROUP[1:]
    if not (len(ux)==len(uy)==len(w)==33): raise AssertionError()
    if min(ux+uy) != LO: raise AssertionError((min(ux+uy),LO))
    if sorted(zip(ux,uy,w)) != sorted(zip(uy,ux,w)): raise AssertionError('reflection failure')
    return ux,uy,w


def axis_values(xs, u):
    """P=e^(x-u)density and CDF for u+Gamma(32)+c Bin(32,p), enclosed."""
    ec = expiv(pt(-C)); den = powi(add(pt(O),ec), B)
    epows = [pt(O)]
    for _ in range(B): epows.append(mul(epows[-1],ec))
    Ps, Fs = [], []
    for x in xs:
        prow, frow = [], []
        for uj in u:
            sump, sumf = pt(Z), pt(Z)
            for k in range(B+1):
                z = exact(lambda k=k,uj=uj: x-uj-C*Decimal(k))
                if z <= Z: continue
                zI = pt(z); power, series = pt(O), pt(Z)
                pterm = None
                for j in range(B):
                    term = div(power, pt(FACT[j])); series = add(series, term)
                    if j == B-1: pterm = term
                    power = mul(power, zI)
                sump = add(sump, scal(pterm, BIN[k]))
                tail = expiv(pt(exact(lambda uj=uj: uj-x)))
                # exp(-ck) - exp(-(x-u)) sum_{j<32} z^j/j!
                sumf = add(sumf, scal(posdiff(epows[k], mul(tail,series)), BIN[k]))
            prow.append(div(sump, den)); frow.append(div(sumf, den))
        Ps.append(prow); Fs.append(frow)
    return Ps, Fs


def cell_sum(Px, Py, dFx, dFy, ex, coeff, weights, q, m):
    """Returns an interval enclosing e^tau int_[LO,HI]^2 q(S)dmu from below."""
    total = pt(Z)
    for i in range(m):
        if i % 16 == 0: print(json.dumps({'phase':'rectangles','row':i,'rows':m}), flush=True)
        for j in range(i, m):
            # S is coordinatewise increasing because every P and coefficient
            # is nonnegative.  Thus its rectangle range is [al,au].
            with localcontext(DOWN):
                alv=Z; isv=Z
                for r in range(33):
                    alv += coeff[r].lo*Px[i][r].lo*Py[j][r].lo
                    isv += weights[r]*dFx[i][r].lo*dFy[j][r].lo
                alv=+alv; isv=+isv
            with localcontext(UP):
                auv=Z; ish=Z
                for r in range(33):
                    auv += coeff[r].hi*Px[i+1][r].hi*Py[j+1][r].hi
                    ish += weights[r]*dFx[i][r].hi*dFy[j][r].hi
                auv=+auv; ish=+ish
            L,U = max(Z,alv), auv
            qL,qU = q.value(L),q.value(U)
            if U == L:
                slope, intercept = pt(Z), pt(qL)
            else:
                slope = div(pt(exact(lambda:qU-qL)), pt(exact(lambda:U-L)))
                intercept = sub(pt(qL), mul(slope,pt(L)))
            # e^tau times exact base-measure mass of the cell.
            mass = mul(expiv(pt(TAU)), mul(ex[i],ex[j]))
            cell = add(mul(slope, IV(isv,ish)), mul(intercept,mass))
            total = add(total, cell if i == j else scal(cell, Decimal(2)))
    return total


def tail_loss(Fx, Fy, w):
    # q(S)>=-S since Phi^*_[0,2](S)<=2S.  The complement of the square is
    # bounded by the union of its two coordinate tails for each kernel.
    ans = pt(Z)
    for r in range(33):
        tx = IV(max(Z,down(lambda r=r: O-Fx[r].hi)), max(Z,up(lambda r=r: O-Fx[r].lo)))
        ty = IV(max(Z,down(lambda r=r: O-Fy[r].hi)), max(Z,up(lambda r=r: O-Fy[r].lo)))
        ans = add(ans, scal(add(tx,ty),w[r]))
    return ans


def focused_grid(m):
    """Exact nonuniform mesh: resolve the only non-affine q region densely.

    The Gamma--Bernoulli laws have their material mass around 32+c*32/(1+e^c)
    plus an anchor shift.  The active soft q(S) region is consequently in the
    central coordinate band.  The outer bands are retained (not discarded)
    but use coarse cells; their q contribution is overwhelmingly in an exact
    affine branch.  This is a predetermined mesh, not an optimizer output.
    """
    if m % 8: raise ValueError('focused mesh requires m divisible by 8')
    # These binary fractions make every grid point a terminating exact
    # decimal: 16/ (m/8), 60/(3m/4), and 90/(m/8) all terminate for m=2^k.
    counts=(m//8, 3*m//4, m//8)
    ends=(LO, Decimal('10'), Decimal('70'), HI)
    xs=[]
    for block,(a,b,k) in enumerate(zip(ends[:-1],ends[1:],counts)):
        step=exact(lambda a=a,b=b,k=k:(b-a)/Decimal(k))
        part=[exact(lambda i=i,a=a,step=step:a+Decimal(i)*step) for i in range(k+1)]
        xs += part if block==0 else part[1:]
    if len(xs)!=m+1 or xs[0]!=LO or xs[-1]!=HI: raise AssertionError('mesh construction')
    return xs


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--m',type=int,default=256,help='uniform grid per axis (power of two recommended)')
    ap.add_argument('--out',type=Path,default=Path('rounds/r061/certificate61.json'))
    args=ap.parse_args()
    if args.m < 8: raise SystemExit('m must be >=8')
    if args.m & (args.m-1): raise SystemExit('m must be a power of two')
    if args.m % 8: raise SystemExit('m must be divisible by 8')
    source=Path(__file__); source_sha=hashlib.sha256(source.read_bytes()).hexdigest()
    ux,uy,W = geometry()
    xs=focused_grid(args.m)
    Px,Fx=axis_values(xs,ux); Py,Fy=axis_values(xs,uy)
    dFx=[ [posdiff(Fx[i+1][r],Fx[i][r]) for r in range(33)] for i in range(args.m) ]
    dFy=[ [posdiff(Fy[i+1][r],Fy[i][r]) for r in range(33)] for i in range(args.m) ]
    ex=[posdiff(expiv(pt(-xs[i])), expiv(pt(-xs[i+1]))) for i in range(args.m)]
    etau=expiv(pt(TAU))
    coeff=[mul(pt(W[r]), expiv(pt(exact(lambda r=r: ux[r]+uy[r]-TAU)))) for r in range(33)]
    q=QUnder()
    window=cell_sum(Px,Py,dFx,dFy,ex,coeff,W,q,args.m)
    loss=tail_loss(Fx[-1],Fy[-1],W)
    lower=sub(window,loss)
    # Exact finite cone-union volume, included only as a separately certified
    # geometry value for the stated 33 anchors.
    points=sorted(zip(ux,uy)); floor=None; vol=pt(Z); i=0
    while i<len(points):
        x=points[i][0]; j=i; newfloor=floor
        while j<len(points) and points[j][0]==x:
            newfloor=points[j][1] if newfloor is None else min(newfloor,points[j][1]); j+=1
        nx=points[j][0] if j<len(points) else None
        floor=newfloor
        xmass=expiv(pt(-x)) if nx is None else posdiff(expiv(pt(-x)),expiv(pt(-nx)))
        vol=add(vol,mul(xmass,expiv(pt(-floor))))
        i=j
    target=Decimal('0.1000')
    result={
      'status':'certified' if lower.lo >= target else 'valid directed lower-bound run; target not reached',
      'scope':'Finite 33-node gamma=1/5 parabola only: a rigorous lower bound for the cap-2 Fenchel source cost with all 33 target-one constraints. No continuous/all-v family claim and no endpoint-constant claim.',
      'arithmetic':'90-digit directed Decimal intervals; finite input decimals are exact rationals; Decimal exp/ln endpoints are expanded outward by one representable number; all algebraic interval operations are directed.',
      'parameters_exact_decimal':{'N':N,'B':B,'cap':'2','target':'1','c':str(C),'tau':str(TAU),'gamma':str(GAMMA),'window':[str(LO),str(HI)],'weights_grouped':list(map(str,WGROUP))},
      'grid_per_axis':args.m,'mesh':'predetermined three-band exact grid [-15.2,10],[10,70],[70,160] with cell counts m/8,3m/4,m/8','q_chord_hgrid':512,
      'window_dual_scaled_interval':[str(window.lo),str(window.hi)],
      'tail_loss_upper':str(loss.hi),
      'dual_scaled_lower':str(lower.lo),
      'dual_scaled_interval':[str(lower.lo),str(lower.hi)],
      'target_dual_scaled':'0.1000','target_met':lower.lo>=target,
      'source_cost_lower':'exp(-tau)*('+str(lower.lo)+')',
      'cone_union_volume_interval':[str(vol.lo),str(vol.hi)],
      'cone_union_volume_lower':str(vol.lo),
      'source_sha256':source_sha,
      'method':'Fenchel weak duality; q(S)=S-Phi*; concave chord minorants on monotone rectangles; exact Gamma--Bernoulli CDF increments for the linear S term; union-bound tail using q>=-S.'
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    if not result['target_met']: raise SystemExit('valid run did not reach 0.1000')

if __name__=='__main__': main()
