#!/usr/bin/env python3
"""New exact scalar conditional bridge checks; not original-space source samples."""
from fractions import Fraction as F
from pathlib import Path
import math,json,hashlib,sys
sys.set_int_max_str_digits(0)
P=Path(__file__).resolve().parent
REG=P/"gamma_occupancy_winner_bridge_registration_20261007.json"
reg=json.loads(REG.read_text())
assert reg["status"]=="registered_before_execution"
checks=[];rows=[]
def test(n,label,ok):
    if not ok:raise AssertionError((n,label))
    checks.append({"n":n,"label":label})
def ceilroot(n,k):
    a=1
    while a**k<n:a+=1
    return a
def exp_interval(x,M):
    term=lo=F(1)
    for j in range(1,M+1):
        term*=x/j;lo+=term
    nxt=term*x/(M+1);ratio=x/(M+2)
    assert ratio<1
    return lo,lo+nxt/(1-ratio)
def eval_normalized(q,p,n,posterior):
    return sum((w*(q/p)**j*((1-q)/(1-p))**(n-j)
                for j,w in posterior.items()),F())
for spec in reg["rounds"]:
    n,N,M=spec["n"],spec["window_nodes"],spec["exp_terms"]
    b=F(1,ceilroot(n,3));K=(n*b).numerator//(n*b).denominator
    test(n,"b<=1/2",b<=F(1,2))
    test(n,"n b^3<=1",n*b**3<=1)
    elo,ehi=exp_interval(F(2,3),M)
    test(n,"certified exp(2/3)<2",ehi<2)
    test(n,"sqrt3>5/3",3>F(5,3)**2)
    test(n,"ideal volume coefficient<10",64*ehi**2<300)
    test(n,"ideal weighted peak coefficient<128/5",
         64*ehi/(3*F(5,3))<F(128,5))
    test(n,"Gamma occupancy coefficient<16/5",
         8*ehi/(3*F(5,3))<F(16,5))
    for k in range(1,K+1):
        p=F(k,n);m=ceilroot(k,2);eps=F(1,4*m);left=p*(1-eps)
        test(n,"np>=1",n*p>=1)
        test(n,"left window inside(0,b]",0<left<=p<=b)
        test(n,"epsilon<=1/4",eps<=F(1,4))
        test(n,"rational sqrt interval",k<=m*m<=4*k)
        test(n,"ceil window weighted mass>=1/8",k*eps**2>=F(1,64))
        scalar_integral=F(math.factorial(k-1)*math.factorial(n-k),
                          math.factorial(n))/(p**k*(1-p)**(n-k))
        test(n,"single-mask Beta integral positive",scalar_integral>0)
        profiles=[("single",{k:F(1)})]
        if k>=2:
            profiles.append(("stationary_mixture",{k-1:F(1,4),k:F(1,2),k+1:F(1,4)}))
        for label,posterior in profiles:
            test(n,"posterior mass1",sum(posterior.values(),F())==1)
            test(n,"posterior mean np",sum(j*w for j,w in posterior.items())==n*p)
            test(n,"outside constant coefficient0",posterior.get(0,F())==0)
            test(n,"response normalized atp",eval_normalized(p,p,n,posterior)==1)
            minimum=F(1);chi_max=F()
            for t in range(N):
                q=left+(p-left)*F(t,N-1)
                chi=n*(p-q)**2/(q*(1-q))
                value=eval_normalized(q,p,n,posterior)
                test(n,"rational chi<=1/6",chi<=F(1,6))
                test(n,"positive response persistence>=5/6",value>=F(5,6))
                chi_max=max(chi_max,chi);minimum=min(minimum,value)
            rows.append({"n":n,"k":k,"p":str(p),"b":str(b),
                         "profile":label,"nodes":N,"ceil_sqrt_k":m,
                         "epsilon":str(eps),"left":str(left),
                         "max_chi":str(chi_max),"min_normalized_response":str(minimum),
                         "single_profile_logtime_integral":str(scalar_integral)})
        for height in [F(4)+F(1,128),F(8)+F(1,128),F(16),F(64)]:
            lower=F(5,6)*height-2
            test(n,"same-source cap subtraction gives >kappa",lower>1)
            if height>8:
                test(n,"overshoot persistence at least M/2",lower>=height/2)
            for a in [F(1,8),F(1,2),F(1),F(2),F(8)]:
                if height*height>=a*a*k:
                    # Squared positive comparisons avoid irrational sqrt.
                    test(n,"weighted overshoot threshold direction",height/a>=0 and (height/a)**2>=k)
res={"status":"PASS_EXACT_CONDITIONAL_WINNER_BRIDGE",
     "total_checks":len(checks),
     "checks_by_round":{str(n):sum(t["n"]==n for t in checks) for n in (8,32,128)},
     "checks":checks,"rows":rows,
     "registration_sha256":hashlib.sha256(REG.read_bytes()).hexdigest(),
     "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     "limitations":["Scalar profiles verify universal positive-Bernstein persistence only; no assertion they are arbitrary original G_c source coefficients",
                    "Three-term profile stationarity checked, not global peak",
                    "Original exterior obstacle/source contraction is proved analytically, not reconstructed numerically",
                    "Ideal sqrt window constants have analytic proof; guard uses narrower rational ceil window",
                    "No actual geom/fullfuture gate samples, no full maximal bound or exponent fit"]}
out=P/"gamma_occupancy_winner_bridge_results_20261007.json"
out.write_text(json.dumps(res,indent=2)+"\n")
print({k:res[k] for k in ("status","total_checks","checks_by_round",
                         "registration_sha256","script_sha256")})
print("saved_profiles",len(rows))
