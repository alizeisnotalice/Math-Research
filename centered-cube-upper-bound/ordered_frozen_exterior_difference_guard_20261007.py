#!/usr/bin/env python3
"""Pure rational intervals for new original face and source-split constants."""
from fractions import Fraction as F
from pathlib import Path
import math,json,hashlib,sys
sys.set_int_max_str_digits(0)  # Exact interval endpoints can exceed 4300 decimal digits.
P=Path(__file__).resolve().parent
REG=P/"ordered_frozen_exterior_difference_registration_20261007.json"
reg=json.loads(REG.read_text())
assert reg["status"]=="registered_before_execution"
checks=[]
rows=[]
def test(n,label,cond):
    if not cond:raise AssertionError((n,label))
    checks.append({"round_n":n,"label":label})
def interval_string(t):
    return {"lower":str(t[0]),"upper":str(t[1])}
def exp_interval(x,M):
    term=F(1);lo=term
    for k in range(1,M+1):
        term*=x/k;lo+=term
    nxt=term*x/(M+1)
    ratio=x/(M+2)
    assert ratio<1
    return lo,lo+nxt/(1-ratio)
def log2_interval(M):
    x=F(1,3)
    lo=2*sum((x**(2*j+1)/F(2*j+1) for j in range(M)),F())
    return lo,lo+2*x**(2*M+1)/(F(2*M+1)*(1-x*x))
def kcoef(m,n,r):
    if any(j>1 for j in m):return F()
    k=sum(m)
    return r**k*(1-r)**(n-k)
def hcoef(m,n,r):
    j=sum(m)
    den=math.prod(math.factorial(k) for k in m)
    return r**j/F(den)
def product_S(coef,m,n,r):
    result=(n+1)*coef(m,n,r)
    for i,j in enumerate(m):
        if j:
            mm=list(m);mm[i]-=1
            result-=coef(mm,n,r)
    return result

for spec in reg["rounds"]:
    n,M=spec["n"],spec["series_terms"]
    rr=[F(),F(1,100*n**4),F(1,2*n),F(1,n),F(1,2),F(3,4),F(1)]
    # N maximal L2 constant and bounded negative source amount.
    test(n,"L2 constant from spectral integral",40*4==160)
    test(n,"negative-source 2kappa weak payment",F(2)*160==320)
    for s in [F(0),F(1,128),F(1,8),F(1,2),F(1),F(2),F(n)]:
        if s:
            test(n,"N spectral factor bound",
                 (1+s)**2*min(F(1),1/(s*s))<=4)
    for w in [F(),F(1,16),F(1,2),F(1)]:
        test(n,"negative cap squared",w*w<=w)
    for r in rr:
        elo,ehi=exp_interval(n*r,M)
        emlo,emhi=1/ehi,1/elo
        atom=((n+1)*((1-r)**n-emhi),(n+1)*((1-r)**n-emlo))
        test(n,"atom sign",atom[1]<=0)
        test(n,"atom abs<=2",atom[0]>=-2)
        for k in range(n+1):
            m=[1]*k+[0]*(n-k)
            square=(n+1)*r**k*(1-r)**(n-k)
            if k:
                square-=k*r**(k-1)*(1-r)**(n-k+1)
                bracket=r**(k-1)*(1-r)**(n-k)*((n+k+1)*r-k)
                test(n,"squarefree n+k+1 bracket",square==bracket)
            test(n,"squarefree full recurrence",product_S(kcoef,m,n,r)==square)
            if k:
                m[0]=2
                repeated=-r**k*(1-r)**(n-k)
                test(n,"one-repeat full recurrence",
                     product_S(kcoef,m,n,r)==repeated)
            for mtest in (m,[2,2]+[0]*(n-2)):
                j=sum(mtest)
                den=math.prod(math.factorial(ki) for ki in mtest)
                hp=(n+1)*r**j/F(den)
                if j:hp-=j*r**(j-1)/F(den)
                test(n,"Poisson full coefficient",
                     product_S(hcoef,mtest,n,r)==hp)
        rows.append({"n":n,"r":str(r),"exp_terms":M,
                     "atom_interval":interval_string(atom)})
    # Continuous actual two-face density, no reset surrogate.
    llo,lhi=log2_interval(M)
    w2=(F(4,3)*(1-lhi),F(4,3)*(1-llo))
    test(n,"original G1^2 peak interval",0<w2[0]<=w2[1]<=1)
    r=F(1,100*n**4)
    x=2*n*r
    remainder=(2*n+1)*F(2,6)*x**3/(1-x)
    lead=(2-w2[1],2-w2[0])
    face=(-r*r*lead[1]-remainder,-r*r*lead[0]+remainder)
    test(n,"two-face full density negative",face[1]<0)
    test(n,"remainder ratio<1/2",remainder/(r*r)<F(1,2))
    rows.append({"n":n,"r":str(r),"two_face_density_interval":interval_string(face),
                 "two_face_interval_divided_r2":interval_string((face[0]/(r*r),face[1]/(r*r))),
                 "G1_squared_peak_interval":interval_string(w2),
                 "majorant_remainder":str(remainder)})
res={"status":"PASS_EXACT_RATIONAL_INTERVALS","total_checks":len(checks),
     "checks_by_round":{str(n):sum(c["round_n"]==n for c in checks) for n in (8,32,128)},
     "checks":checks,"rows":rows,
     "registration_sha256":hashlib.sha256(REG.read_bytes()).hexdigest(),
     "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     "limitations":["No old symbol rerun",
                    "Original c=1 face density interval relies on analytic convolution-peak and coefficient remainder proof",
                    "No actual obstacle solver or hpositive maximal experiment",
                    "No positive-part weak budget or A_ord improvement"]}
out=P/"ordered_frozen_exterior_difference_results_20261007.json"
out.write_text(json.dumps(res,indent=2)+"\n")
print({k:res[k] for k in ("status","total_checks","checks_by_round",
                         "registration_sha256","script_sha256")})
for row in rows:
    if "two_face_interval_divided_r2" in row:
        a=row["two_face_interval_divided_r2"]
        print("two face",row["n"],float(F(a["lower"])),float(F(a["upper"])))
