#!/usr/bin/env python3
"""Portable directed-Decimal verifier for the explicit scalar54 staircase."""
from __future__ import annotations
import argparse,hashlib,json,math
from dataclasses import dataclass
from decimal import Decimal,Context,Inexact,ROUND_CEILING,ROUND_FLOOR,localcontext
from pathlib import Path
P=110; DOWN=Context(prec=P,rounding=ROUND_FLOOR); UP=Context(prec=P,rounding=ROUND_CEILING); NEAR=Context(prec=P); EXACT=Context(prec=500)
Z=Decimal(0); O=Decimal(1); N=64; M=512
def ex(f):
 with localcontext(EXACT) as c:
  c.clear_flags();r=+f()
  if c.flags[Inexact]:raise ArithmeticError('non-exact rational preprocessing')
  return r
def dn(f):
 with localcontext(DOWN):return +f()
def up(f):
 with localcontext(UP):return +f()
def outward(x,hi):
 with localcontext(UP if hi else DOWN):return x.next_plus() if hi else x.next_minus()
@dataclass(frozen=True)
class IV:
 lo:Decimal;hi:Decimal
 def __post_init__(self):
  if not(self.lo.is_finite() and self.hi.is_finite() and self.lo<=self.hi):raise ValueError((self.lo,self.hi))
def pt(x):
 if not x.is_finite():raise ValueError(x)
 return IV(x,x)
def add(a,b):return IV(dn(lambda:a.lo+b.lo),up(lambda:a.hi+b.hi))
def sub(a,b):return IV(dn(lambda:a.lo-b.hi),up(lambda:a.hi-b.lo))
def mul(a,b):return IV(min(dn(lambda x=x,y=y:x*y)for x in(a.lo,a.hi)for y in(b.lo,b.hi)),max(up(lambda x=x,y=y:x*y)for x in(a.lo,a.hi)for y in(b.lo,b.hi)))
def div(a,b):
 if b.lo<=Z<=b.hi:raise ZeroDivisionError(b)
 return IV(min(dn(lambda x=x,y=y:x/y)for x in(a.lo,a.hi)for y in(b.lo,b.hi)),max(up(lambda x=x,y=y:x/y)for x in(a.lo,a.hi)for y in(b.lo,b.hi)))
def powi(a,n):
 if n<0 or a.lo<Z:raise ValueError((a,n))
 r=pt(O);q=a
 while n:
  if n&1:r=mul(r,q)
  n//=2
  if n:q=mul(q,q)
 return r
def expiv(a):
 with localcontext(NEAR):x=a.lo.exp();y=a.hi.exp()
 return IV(outward(x,False),outward(y,True))
def lniv(a):
 if a.lo<=Z:raise ValueError(a)
 with localcontext(NEAR):x=a.lo.ln();y=a.hi.ln()
 return IV(outward(x,False),outward(y,True))
def dec(x):
 s=str(x).split('`')[0];v=Decimal(s)
 if not v.is_finite():raise ValueError(s)
 return v
def pair(a):return {'lo':str(a.lo),'hi':str(a.hi)}
def gpoly(x):
 t=pt(O);s=pt(O);xx=pt(x)
 for j in range(1,N):t=div(mul(t,xx),pt(Decimal(j)));s=add(s,t)
 return s
def Phi(h):
 if h<=O:return pt(h)
 return mul(pt(h),powi(add(pt(O),lniv(pt(h))),N-1))
def main():
 d=Path(__file__).parent;ap=argparse.ArgumentParser();ap.add_argument('--infile',default=str(d.parent/'rounds/r054/scalar54.json'));ap.add_argument('--out',default='/tmp/r054_scalar54.json');a=ap.parse_args()
 raw=json.loads(Path(a.infile).read_text());c=dec(raw['c_exact_decimal']);tau=dec(raw['tau_metadata_exact_decimal']);bs=[dec(x)for x in raw['boundaries_exact_decimal']]
 if not(c>Z and len(bs)==M+1 and bs[0]>=Z and all(bs[j]<bs[j+1]for j in range(M))):raise ValueError('parameter guards')
 hs=[O]+[ex(lambda j=j:O+Decimal(j)/Decimal(M))for j in range(1,M+1)];step=ex(lambda:O/Decimal(M));ec=expiv(pt(c.copy_negate()));den=powi(add(pt(O),ec),N)
 ep=[pt(O)]
 for _ in range(N):ep.append(mul(ep[-1],ec))
 eb=[expiv(pt(b.copy_negate()))for b in bs];tails=[]
 for b,e in zip(bs,eb):
  s=pt(Z)
  for k in range(N+1):
   x=ex(lambda b=b,k=k:b-c*Decimal(k))
   term=mul(pt(Decimal(math.comb(N,k))),ep[k]) if x<=Z else mul(pt(Decimal(math.comb(N,k))),mul(e,gpoly(x)))
   s=add(s,div(term,den))
  tails.append(s)
 mean=tails[0]
 for t in tails[1:]:mean=add(mean,mul(pt(step),t))
 J=pt(Z)
 for j,h in enumerate(hs):J=add(J,mul(Phi(h),sub(eb[j],eb[j+1])if j<M else eb[j]))
 L=Decimal('100000000');alpha=ex(lambda:O-Decimal(N-1)/L)
 eab=[expiv(pt(ex(lambda b=b:alpha*b).copy_negate()))for b in bs];Jalpha=pt(Z)
 for j,h in enumerate(hs):Jalpha=add(Jalpha,mul(Phi(h),sub(eab[j],eab[j+1])if j<M else eab[j]))
 Jalpha=div(Jalpha,pt(alpha))
 rate=dec('71.79');threshold=expiv(pt(rate.copy_negate()));margin=ex(lambda:mean.lo-O)
 ans={'status':'directed Decimal certificate for explicit scalar n=64 staircase','scope':'An explicit finite-L input is certified below; no optimality claim.','arithmetic':'110-digit directed intervals; elementary exp/ln one-ulp outward; all input decimals exact rationals and every preprocessing arithmetic is Inexact-trapped.','n':N,'cap':'2','levels':M,'height_rule':'0 below b0; 1+j/512 on [bj,b(j+1)), j=0..511; 2 on [b512,infinity)','c_exact_decimal':str(c),'tau_metadata_exact_decimal':str(tau),'boundaries_exact_decimal':list(map(str,bs)),'boundary_sha256':hashlib.sha256(('\n'.join(map(str,bs))).encode()).hexdigest(),'mean_interval':pair(mean),'budget_J_interval':pair(J),'finite_L':str(L),'alpha':str(alpha),'finite_L_budget_Jalpha_interval':pair(Jalpha),'exp_minus_71_79_interval':pair(threshold),'mean_lower_minus_1':str(margin),'mean_lower_exceeds_1_000001':mean.lo>dec('1.000001'),'budget_upper_below_exp_minus_71_79':J.hi<threshold.lo,'finite_L_budget_upper_below_exp_minus_71_79':Jalpha.hi<threshold.lo,'claimed_logratio_strictly_exceeds':'71.79','finite_input_transfer':'For f_L(x)=H(sum_i -log|x_i|-L)1_[-1,1]^64 with L=1e8, fixed-r rectangles give average >=E H(T)>1 on volume K_L. Its Orlicz integral is <=K_L Jalpha because (1+t/L)^63<=exp(63t/L). Thus its actual ratio is >=1/Jalpha.','code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_parameter_sha256':hashlib.sha256(json.dumps({'c':str(c),'tau':str(tau),'b':list(map(str,bs))},sort_keys=True).encode()).hexdigest(),'checks':{'tail':'finite Gamma-binomial upper-tail sum','mean':'P(T>=b0)+sum_{j=1}^512 P(T>=bj)/512','budget':'exact step integral against e^-t','finite_L_budget':'Jalpha=(1/alpha) times the same explicit step sum with exp(-alpha*b)'}}
 if not(ans['mean_lower_exceeds_1_000001'] and ans['budget_upper_below_exp_minus_71_79'] and ans['finite_L_budget_upper_below_exp_minus_71_79']):raise ArithmeticError(ans)
 Path(a.out).write_text(json.dumps(ans,indent=2)+'\n');print(json.dumps({'mean':ans['mean_interval'],'J':ans['budget_J_interval'],'e^-71.79':ans['exp_minus_71_79_interval'],'pass':True},indent=2))
if __name__=='__main__':main()
