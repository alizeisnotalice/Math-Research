#!/usr/bin/env python3
"""Registered small original-symbol diagnostic. Not interval arithmetic."""
from pathlib import Path
from fractions import Fraction as F
import hashlib, json, math
import numpy as np
P=Path(__file__).resolve().parent
REG=P/"ordered_frozen_spectral_difference_registration_20261007.json"
reg=json.loads(REG.read_text())
assert reg["status"]=="registered_before_execution"
exact=[]
numeric=[]
def ex(label,ok):
    assert ok,label
    exact.append(label)
def nu(label,ok):
    assert ok,label
    numeric.append(label)
ex("small integral27/4",F(math.factorial(5),2**6)+4*F(math.factorial(4),2**5)+5*F(math.factorial(3),2**4)==F(27,4))
ex("late33",F(1)+4*F(2)+4*F(6)==33)
ex("sum<40",F(27,4)+33<40)
ex("small S late sum<40",F(27,4)+17<40)
for n in reg["dimensions"]:
    for power in (2,3,4):
        h=F(n**power)
        bound=(h-1)**2/(n*h)
        ex(f"spike{n}^{power} energy>kappaW",bound>1)
        ex(f"spike{n}^{power} exact expansion",bound==(h-2+1/h)/n)
    h=F(10**6)
    ex(f"spike{n} fixedheight",(h-1)**2/(n*h)>100)
    for k in range(1,11):
        ex(f"factorial{k}moment",F(math.factorial(k),2**(k+1))>0)

def B(x):
    x=np.asarray(x,dtype=float)
    ans=np.zeros_like(x)
    q=x>0
    small=q & (x<1e-4)
    # Stable local expansion B=x/2-x^2/12+x^3/24-19x^4/720.
    z=x[small]
    ans[small]=z/2-z*z/12+z**3/24-19*z**4/720
    big=q & ~small
    ans[big]=x[big]/np.log1p(x[big])-1
    return ans
def frequencies(n,pattern):
    if pattern=="zero":return np.zeros(n)
    if pattern=="all frequency 1/1000":return np.full(n,.001)
    if pattern=="all frequency 1":return np.ones(n)
    if pattern=="all frequency64":return np.full(n,64.)
    return np.resize(np.array([0.,.125,1.,8.,64.]),n)
def calc(a,r):
    S=float(a.sum())
    r=np.asarray(r,dtype=float)
    if S==0:return np.zeros_like(r),np.zeros_like(r)
    x=r[:,None]*a[None,:]
    d=-np.log1p(-x)-x
    small=x<1e-4
    d[small]=sum(x[small]**k/k for k in range(2,8))
    D=d.sum(axis=1)
    RD=(x*x/(1-x)).sum(axis=1)
    u=r*S
    e=np.exp(-D)
    om=-np.expm1(-D)
    m=-math.sqrt(S)*np.exp(-u)*om
    rm=math.sqrt(S)*np.exp(-u)*(u*om-RD*e)
    return m,rm
def integrate(a,order):
    nodes,weights=np.polynomial.legendre.leggauss(order)
    total=0.
    for lo,hi in zip(np.linspace(-40,0,17)[:-1],np.linspace(-40,0,17)[1:]):
        log_r=(hi+lo)/2+(hi-lo)*nodes/2
        m,rm=calc(a,np.exp(log_r))
        total+=(hi-lo)/2*float(np.dot(weights,m*m+rm*rm))
    return total
rows=[]
for n in reg["dimensions"]:
    for cstr in reg["c"]:
        c=float(F(cstr))
        for pat in reg["patterns"]:
            b=B(frequencies(n,pat)**2)
            a=c*b/(1+c*b)
            S=float(a.sum()); Q=float(np.dot(a,a))
            ints=[integrate(a,o) for o in (16,32,64)]
            bound=0 if S==0 else 40*min(S,1/S)
            r0=math.exp(-40)
            tail=S*Q*Q*r0**4*(1+(r0*S+2)**2)/4
            nu(f"{n}/{cstr}/{pat}:spectrum",bool(np.all((a>=0)&(a<=1))))
            nu(f"{n}/{cstr}/{pat}:Q",Q<=min(S,S*S)+1e-12)
            nu(f"{n}/{cstr}/{pat}:integral",ints[-1]+tail<=bound+1e-12)
            nu(f"{n}/{cstr}/{pat}:nonnegative",ints[-1]>=-1e-14)
            rr=np.array([1/1024,1/16,.25,.5,.75,1-1/1024])
            m,rm=calc(a,rr)
            K=np.prod(1-rr[:,None]*a[None,:],axis=1)
            KP=-K*np.sum(a[None,:]/(1-rr[:,None]*a[None,:]),axis=1)
            H=np.exp(-rr*S)
            em=math.sqrt(S)*(K-H)
            erm=math.sqrt(S)*rr*(KP+S*H)
            residual=max(float(np.max(np.abs(m-em))),float(np.max(np.abs(rm-erm))))
            nu(f"{n}/{cstr}/{pat}:joint derivative",residual<1e-10)
            rows.append({"n":n,"c":cstr,"pattern":pat,"S":S,"Q":Q,
                         "integral_GL16_32_64":ints,"bound":bound,"small_r_tail_bound":tail,
                         "nested_differences":[abs(ints[1]-ints[0]),abs(ints[2]-ints[1])],
                         "derivative_residual":residual})
res={"status":"PASS_EXACT_CONSTANTS_AND_ORIGINAL_SYMBOL_DIAGNOSTIC",
     "exact_count":len(exact),"numeric_count":len(numeric),"total_checks":len(exact)+len(numeric),
     "exact_checks":exact,"numeric_checks":numeric,"rows":rows,
     "registration_sha256":hashlib.sha256(REG.read_bytes()).hexdigest(),
     "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     "limitations":["Double original-symbol quadrature, not interval",
                    "Low-symbol Taylor evaluation is a diagnostic, not interval",
                    "Spike guard certifies a necessary lower bound, no obstacle solver",
                    "No source-energy upper, no A_ord or cube endpoint conclusion"]}
OUT=P/"ordered_frozen_spectral_difference_results_20261007.json"
OUT.write_text(json.dumps(res,indent=2)+"\n")
print(json.dumps({k:res[k] for k in ("status","exact_count","numeric_count","total_checks",
                                    "registration_sha256","script_sha256")},indent=2))
print("max nested difference",max(max(r["nested_differences"]) for r in rows))
print("max derivative residual",max(r["derivative_residual"] for r in rows))
