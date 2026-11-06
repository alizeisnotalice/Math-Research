#!/usr/bin/env python3
"""Directed-decimal lower certificate for the r052 fixed 65-receiver D1 test.

Scope: the fixed cap-2, target-one 65-receiver geometry from r053 only.
It proves a lower bound for its Fenchel dual diagnostic D1; it is neither a
strong-maximal halo estimate nor a bound for C_N.

All JSON finite decimals are interpreted as exact decimal rationals.  Every
non-rational elementary operation is enclosed outwards.  Decimal exp/ln are
used at 90 digits and each returned endpoint is expanded by one representable
decimal in the required direction.
"""
from __future__ import annotations

import argparse
import bisect
import hashlib
import json
import math
import multiprocessing as mp
from dataclasses import dataclass
from decimal import Decimal, Context, Inexact, ROUND_CEILING, ROUND_FLOOR, localcontext
from pathlib import Path
from typing import Iterable

ROOT = Path('/Users/zhengzhihao/Desktop/T')
REPO = ROOT/'Math-Research-strong-maximal/strong-maximal-endpoint'
FAIR = REPO/'rounds/r052/fair52.json'
SEED = REPO/'rounds/r030/receiver_seed.json'
OUTDIR = ROOT/'work/r054'

PREC = 90
DOWN = Context(prec=PREC, rounding=ROUND_FLOOR)
UP = Context(prec=PREC, rounding=ROUND_CEILING)
NEAR = Context(prec=PREC)  # Decimal ln/exp are correctly rounded; expand one ulp.
EXACT = Context(prec=500)
Z = Decimal(0); O = Decimal(1)
N = 64; B = 32; CAP = Decimal(2); HGRID = 256

def down(f):
    with localcontext(DOWN): return +f()
def up(f):
    with localcontext(UP): return +f()
def exact(f):
    # Finite decimal inputs and their rational combinations are required to
    # stay exact.  Do not merely assume that 500 digits suffice: fail if this
    # preprocessing context reports any inexact rounding.
    with localcontext(EXACT) as ctx:
        ctx.clear_flags()
        answer=+f()
        if ctx.flags[Inexact]:
            raise ArithmeticError('rational preprocessing exceeded EXACT precision')
        return answer
def outward(x: Decimal, upper: bool) -> Decimal:
    with localcontext(UP if upper else DOWN):
        return x.next_plus() if upper else x.next_minus()

@dataclass(frozen=True)
class IV:
    lo: Decimal
    hi: Decimal
    def __post_init__(self):
        if self.lo > self.hi: raise ValueError((self.lo, self.hi))

def pt(x: Decimal) -> IV: return IV(x, x)
def add(x: IV, y: IV) -> IV: return IV(down(lambda:x.lo+y.lo), up(lambda:x.hi+y.hi))
def sub(x: IV, y: IV) -> IV: return IV(down(lambda:x.lo-y.hi), up(lambda:x.hi-y.lo))
def neg(x: IV) -> IV: return IV(-x.hi, -x.lo)
def mul(x: IV, y: IV) -> IV:
    vals_lo=[down(lambda a=a,b=b:a*b) for a in (x.lo,x.hi) for b in (y.lo,y.hi)]
    vals_hi=[up(lambda a=a,b=b:a*b) for a in (x.lo,x.hi) for b in (y.lo,y.hi)]
    return IV(min(vals_lo), max(vals_hi))
def scal(x: IV, a: Decimal) -> IV:
    return mul(x, pt(a))
def div(x: IV, y: IV) -> IV:
    if y.lo <= 0 <= y.hi: raise ZeroDivisionError(y)
    vals_lo=[down(lambda a=a,b=b:a/b) for a in (x.lo,x.hi) for b in (y.lo,y.hi)]
    vals_hi=[up(lambda a=a,b=b:a/b) for a in (x.lo,x.hi) for b in (y.lo,y.hi)]
    return IV(min(vals_lo), max(vals_hi))
def powi(x: IV, n: int) -> IV:
    if n < 0: raise ValueError(n)
    if n == 0: return pt(O)
    # Integer powers are built only from directed interval products.  In
    # particular, this does not trust Decimal.__pow__ to have a usable
    # directed-rounding contract.
    if x.lo < 0:
        raise ValueError('certificate powers require nonnegative base')
    answer, base, exponent = pt(O), x, n
    while exponent:
        if exponent & 1: answer=mul(answer,base)
        exponent >>= 1
        if exponent: base=mul(base,base)
    return answer
def expiv(x: IV) -> IV:
    with localcontext(NEAR): a=x.lo.exp(); b=x.hi.exp()
    return IV(outward(a,False), outward(b,True))
def lniv(x: IV) -> IV:
    if x.lo <= 0: raise ValueError(x)
    with localcontext(NEAR): a=x.lo.ln(); b=x.hi.ln()
    return IV(outward(a,False), outward(b,True))
def posdiff(x: IV, y: IV) -> IV:
    z=sub(x,y)
    return IV(max(Z,z.lo), max(Z,z.hi))

FACT=[Decimal(math.factorial(j)) for j in range(B+1)]
BIN=[Decimal(math.comb(B,k)) for k in range(B+1)]

def djson(path: Path):
    return json.loads(path.read_text(), parse_float=Decimal, parse_int=Decimal)

def finite(x) -> Decimal:
    # json parses integral data as Decimal as well; Decimal is an exact finite rational.
    return x if isinstance(x,Decimal) else Decimal(str(x))

class QUnder:
    """Piecewise chord lower minorant of q(a)=a-Phi^*_[0,2](a)."""
    def __init__(self):
        xs=[Z, O, Decimal(N)]
        ys=[Z, O, O]
        last=Decimal(N)
        for i in range(1,HGRID):
            h=exact(lambda:Decimal(1)+Decimal(i)/Decimal(HGRID))
            hh=pt(h); L=lniv(hh); oneL=add(pt(O),L)
            phi=mul(hh,powi(oneL,N-1))
            apr=mul(powi(oneL,N-2),add(pt(Decimal(N)),L))
            x=apr.lo  # a rational decimal left of the true tangent abscissa
            # Here x=apr.lo >= 64 and q is nonincreasing on [64,infinity).
            # Thus q(x)>=q(a_true); the adverse right endpoint in the
            # negative factor (1-h) already gives a lower value at a_true.
            q_at_x_lower=add(pt(phi.lo),mul(pt(exact(lambda:O-h)),pt(apr.hi))).lo
            if x <= last: raise ArithmeticError(('nonincreasing q node',i,last,x))
            xs.append(x); ys.append(q_at_x_lower); last=x
        h=Decimal(2); L=lniv(pt(h)); oneL=add(pt(O),L)
        phi2=mul(pt(h),powi(oneL,N-1))
        acap=mul(powi(oneL,N-2),add(pt(Decimal(N)),L))
        xcap=acap.hi  # safely to the right of the true cap breakpoint
        if xcap <= last: raise ArithmeticError(('cap ordering',last,xcap))
        xs.append(xcap); ys.append(down(lambda:phi2.lo-xcap))
        self.xs,self.ys,self.phi2lo=xs,ys,phi2.lo
    def value(self,a: Decimal) -> Decimal:
        if a <= Z: return Z
        if a <= O: return a  # q(a)=a on [0,1]
        if a <= Decimal(N): return O  # q(a)=1 on [1,64]
        if a >= self.xs[-1]: return down(lambda:self.phi2lo-a)
        j=bisect.bisect_right(self.xs,a)-1
        x0,x1=self.xs[j],self.xs[j+1]; y0,y1=self.ys[j],self.ys[j+1]
        # A line here can have either sign of slope.  General IV arithmetic,
        # not a chain of ROUND_FLOOR scalar operations, preserves direction.
        da=pt(exact(lambda:a-x0)); dx=pt(exact(lambda:x1-x0))
        slope=div(pt(exact(lambda:y1-y0)),dx)
        return add(pt(y0),mul(slope,da)).lo

def geometry():
    fair=djson(FAIR); seed=djson(SEED)
    v=[finite(x) for x in seed['v_nonnegative']]
    d0=[finite(x) for x in seed['receiver_d']]
    kappa=Decimal('-0.08')
    d=[exact(lambda dd=dd,vv=vv:dd+kappa*(vv/Decimal(16))*(vv/Decimal(16))) for dd,vv in zip(d0,v)]
    # Exact ordering of r053/chord_probe53.py: all (+,-) anchors, then all mirrored (-,+) anchors.
    ux=[Z]+[exact(lambda j=j:d[j]+v[j]) for j in range(1,33)]+[exact(lambda j=j:d[j]-v[j]) for j in range(1,33)]
    uy=[Z]+[exact(lambda j=j:d[j]-v[j]) for j in range(1,33)]+[exact(lambda j=j:d[j]+v[j]) for j in range(1,33)]
    raww=[finite(x) for x in fair['shared_weights']]
    scale=fair['winner_N256_soft64']['scale']
    w0=[exact(lambda x=x:x*scale) for x in raww]
    c0=[finite(x) for x in fair['matched_c_vector']]
    w=w0+w0[1:]; c=c0+c0[1:]
    if not (len(ux)==len(uy)==len(w)==len(c)==65): raise AssertionError()
    if any(x < Z for x in c): raise AssertionError('nonpositive-support argument needs c>=0')
    if any(x < Z for x in w): raise AssertionError('weights must be nonnegative')
    tau=finite(seed['tau_ref'])
    # The list is explicitly invariant under (x,y) swap after a permutation.
    sig=sorted((ux[j],uy[j],w[j],c[j]) for j in range(65))
    sigswap=sorted((uy[j],ux[j],w[j],c[j]) for j in range(65))
    if sig != sigswap: raise AssertionError('receiver symmetry failed')
    return ux,uy,w,c,tau

def axis_values(xs: list[Decimal], u: list[Decimal], c: list[Decimal]):
    """P=e^(x-u) density and F=CDF of u+Gamma(32)+c Bin(32,p_c), enclosed."""
    Ps=[]; Fs=[]
    for x in xs:
        rowP=[]; rowF=[]
        for uj,cj in zip(u,c):
            ec=expiv(pt(-cj)); zden=add(pt(O),ec); den=powi(zden,B)
            sump=pt(Z); sumf=pt(Z); ecpow=pt(O)
            for kk in range(B+1):
                z=exact(lambda:x-uj-cj*Decimal(kk))
                if z > Z:
                    zI=pt(z); power=pt(O); poly=pt(Z); pterm=None
                    # Every integral power is generated by a directed
                    # multiplication, never by Decimal's general power.
                    for jj in range(B):
                        term=div(power,pt(FACT[jj]))
                        poly=add(poly,term)
                        if jj == B-1: pterm=term
                        power=mul(power,zI)
                    sump=add(sump,scal(pterm,BIN[kk]))
                    # This is F_Gamma(x-u-c k) times e^{-c k}:
                    # e^{-c k}-e^{-(x-u)} sum_{j<32}(x-u-c k)^j/j!.
                    # The u shift is essential; P above is deliberately the
                    # e^{x-u}-rescaled density, as in the r053 prototype.
                    tail=expiv(pt(exact(lambda:uj-x)))
                    bracket=posdiff(ecpow,mul(tail,poly))
                    sumf=add(sumf,scal(bracket,BIN[kk]))
                ecpow=mul(ecpow,ec)
            rowP.append(div(sump,den)); rowF.append(div(sumf,den))
        Ps.append(rowP); Fs.append(rowF)
    return Ps,Fs

def positive_delta(vals: list[list[IV]]):
    ans=[]
    for i in range(len(vals)-1):
        row=[]
        for j in range(65): row.append(posdiff(vals[i+1][j],vals[i][j]))
        ans.append(row)
    return ans

STATE={}

def lower_row(task):
    i0,i1=task
    Px,Py,dFx,dFy,ex,Q,coeff=(STATE['Px'],STATE['Py'],STATE['dFx'],STATE['dFy'],
      STATE['ex'],STATE['Q'],STATE['coeff'])
    tauiv=STATE['etau']; weights=STATE['w']; M=STATE['M']; total=pt(Z)
    for i in range(i0,i1):
        # exploit the verified x/y symmetry: only j>=i and double off diagonal cells.
        for j in range(i,M):
            # All factors are nonnegative.  Keep the two directed positive
            # sums scalar: this is equivalent to the interval expression but
            # avoids allocating thousands of short IV objects per cell.
            with localcontext(DOWN):
                alv=Z; iav=Z
                for r in range(65):
                    alv += coeff[r].lo*Px[i][r].lo*Py[j][r].lo
                    iav += weights[r]*dFx[i][r].lo*dFy[j][r].lo
                alv=+alv; iav=+iav
            with localcontext(UP):
                auv=Z; iah=Z
                for r in range(65):
                    auv += coeff[r].hi*Px[i+1][r].hi*Py[j+1][r].hi
                    iah += weights[r]*dFx[i][r].hi*dFy[j][r].hi
                auv=+auv; iah=+iah
            al=pt(alv); au=pt(auv); ia=IV(iav,iah)
            # This is already I_a=e^tau int_cell a(x,y)e^{-x-y}dxdy.
            # The e^{-tau} in a cancels the displayed e^tau exactly.
            ixi=mul(tauiv,mul(ex[i],ex[j]))
            L=max(Z,al.lo); U=au.hi
            mult=Decimal(1 if i==j else 2)
            if U < L:
                raise ArithmeticError(('reversed source enclosure',i,j,L,U))
            if U == L:
                cell=scal(ixi,Q.value(L))
            else:
                qL=Q.value(L); qU=Q.value(U)
                # qL and qU are explicit safe decimal values, while L,U are
                # explicit enclosing endpoints.  Enclose the quotient itself;
                # rounding both numerator and denominator downward is invalid
                # for a positive chord slope.
                numerator=exact(lambda:qU-qL); denominator=exact(lambda:U-L)
                slope=div(pt(numerator),pt(denominator))
                intercept=sub(pt(qL),mul(slope,pt(L)))
                cell=add(mul(slope,ia),mul(intercept,ixi))
            total=add(total,scal(cell,mult))
    return total

def tail_bound(Fx160: list[IV], Fy160: list[IV], weights: list[Decimal]) -> IV:
    ans=pt(Z)
    for r in range(65):
        tx=IV(max(Z,down(lambda r=r:O-Fx160[r].hi)), max(Z,up(lambda r=r:O-Fx160[r].lo)))
        ty=IV(max(Z,down(lambda r=r:O-Fy160[r].hi)), max(Z,up(lambda r=r:O-Fy160[r].lo)))
        ans=add(ans,scal(add(tx,ty),weights[r]))
    return ans

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--m',type=int,default=512); ap.add_argument('--workers',type=int,default=4)
    ap.add_argument('--out',type=Path,default=OUTDIR/'d1_certificate.json')
    ap.add_argument('--snapshot',type=Path,required=False)
    args=ap.parse_args()
    if args.m<8: raise SystemExit('m too small')
    source_bytes=Path(__file__).read_bytes()
    source_sha=hashlib.sha256(source_bytes).hexdigest()
    if args.snapshot is not None:
        if args.snapshot.exists() and args.snapshot.read_bytes()!=source_bytes:
            raise SystemExit('refusing to overwrite a distinct source snapshot')
        if not args.snapshot.exists(): args.snapshot.write_bytes(source_bytes)
    ux,uy,w,c,tau=geometry()
    lo=min(ux+uy); hi=Decimal(160); M=args.m
    # All constructed inputs are finite decimal rationals; the denominator M is a power of two in the intended runs.
    if M & (M-1): raise SystemExit('use power-of-two M for exact finite-decimal grid')
    step=exact(lambda:(hi-lo)/Decimal(M)); xs=[exact(lambda i=i:lo+Decimal(i)*step) for i in range(M+1)]
    if xs[-1] != hi: raise ArithmeticError(('grid endpoint',xs[-1],hi))
    Px,Fx=axis_values(xs,ux,c); Py,Fy=axis_values(xs,uy,c)
    dFx=positive_delta(Fx); dFy=positive_delta(Fy)
    ex=[]
    for i in range(M): ex.append(posdiff(expiv(pt(-xs[i])),expiv(pt(-xs[i+1]))))
    etau=expiv(pt(tau))
    coeff=[]
    for rr in range(65):
        coeff.append(mul(pt(w[rr]),expiv(pt(exact(lambda rr=rr:ux[rr]+uy[rr]-tau)))))
    Q=QUnder()
    STATE.update(Px=Px,Py=Py,dFx=dFx,dFy=dFy,ex=ex,Q=Q,coeff=coeff,etau=etau,w=w,M=M)
    workers=max(1,args.workers); cuts=[]; base=M//workers; start=0
    for z in range(workers):
        end=M if z==workers-1 else start+base
        cuts.append((start,end)); start=end
    if workers==1: pieces=[lower_row(cuts[0])]
    else:
        with mp.get_context('fork').Pool(workers) as pool: pieces=pool.map(lower_row,cuts)
    window=pt(Z)
    for x in pieces: window=add(window,x)
    tail=tail_bound(Fx[-1],Fy[-1],w)
    lower=sub(window,tail)
    result={
      'scope':'Directed-decimal lower bound for the fixed r052 kappa=-0.08 65-receiver cap-2 target-one D1 diagnostic only; not a halo or C_N statement.',
      'inputs_exact':'All JSON finite decimals are treated as exact decimal rationals; geometry uses r052 kappa=-0.08 exactly.',
      'm':M,'workers':workers,'precision_digits':PREC,'window':[str(lo),str(hi)],
      'minorant_window_functional_interval':[str(window.lo),str(window.hi)],
      'tail_loss_upper':str(tail.hi),'D1_lower':str(lower.lo),
      'target_check_D1_gt_1_over_20': lower.lo>Decimal(1)/Decimal(20),
      'fair_sha256':hashlib.sha256(FAIR.read_bytes()).hexdigest(),
      'seed_sha256':hashlib.sha256(SEED.read_bytes()).hexdigest(),
      'source_sha256':source_sha,
      'source_snapshot':str(args.snapshot) if args.snapshot is not None else None,
      'method':'outward Decimal interval arithmetic; exp/ln one-ulp expanded; exact finite Gamma-Bernoulli CDF/density mixtures; concave chord minorants of q; union-bound tail loss'
    }
    args.out.parent.mkdir(parents=True,exist_ok=True); args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if not result['target_check_D1_gt_1_over_20']:
        raise SystemExit('certificate did not reach 1/20')
if __name__=='__main__': main()
