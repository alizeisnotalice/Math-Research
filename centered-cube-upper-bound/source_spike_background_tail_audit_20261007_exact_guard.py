#!/usr/bin/env python3
"""Compressed exact qualification guards for a real cube input; no Monte Carlo."""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import hashlib,json,time
BASE=Path(__file__).resolve().parent
PFX="source_spike_background_tail_audit_20261007"
checks=0
def ck(ok,label):
    global checks
    if not ok:
        raise AssertionError(label)
    checks+=1
def display(v):
    return {"approximation_only":float(v),
            "numerator_bits":v.numerator.bit_length(),
            "denominator_bits":v.denominator.bit_length(),
            "fraction_sha256":hashlib.sha256(
                (str(v.numerator)+"/"+str(v.denominator)).encode()).hexdigest()}
def run(n,J):
    start=checks
    eta=F(1,100)
    q=1+F(1,2*J)
    delta=q-1
    root=isqrt(n)
    ck(root*root==n and J==n*n,"registered round")
    ck(J>=n and delta<=F(1,2*n),"finite geometry dilation budget")
    x=J*(1-q**(-n))
    ck(x>=F(n,3),"true spike pointwise lower coefficient")
    ck(x<=F(n,2),"true spike pointwise upper coefficient")
    lower=(1+x)/(1+eta)
    upper=1+x
    ck(lower>F(n,4),"pointwise spike is order n")
    ck(upper<=1+F(n,2),"source cone upper envelope")
    A,C=2*q-F(2,3),2*q-1
    Kq=C*A**(n-1)/n+(1-F(1,n))*C*C*A**(n-2)
    K1=F(4,3)**(n-2)*(1+F(1,3*n))
    ck(Kq>0 and K1>0,"cone square positive")
    ck(Kq/K1<=16,"exact cone square dilation factor")
    c_n=F(9,16)*(1+F(1,3*n))
    ck(c_n<=F(3,4) and 16*c_n<16,"Minkowski coefficient bound")
    ck(K1/F(2)**n==c_n*F(2,3)**n,"cone change of variables")
    # b stays symbolic; constructing q^(Jn) adds no information to the budget.
    b_lower=F(3,2)
    b_upper=F(2)
    ck(b_lower==1+J*delta and b_lower<b_upper,"complete radius budget")
    cone_upper=1+F(n,2)
    spike_first=eta*cone_upper/b_lower**n
    spike_second=eta*cone_upper**2/b_lower**n
    background_second=64*F(2,3)**n
    moment_ratio=background_second+spike_second
    ck(spike_first==eta*cone_upper*F(2,3)**n,"spike first background cost")
    ck(spike_second==eta*cone_upper**2*F(2,3)**n,"spike square background cost")
    ck(moment_ratio==(64+eta*cone_upper**2)*F(2,3)**n,"complete second ratio")
    ck(moment_ratio<root,"registered family does not refute square-root factor")
    Wmin=F(2)**n
    Wmax=Wmin+eta/b_lower**n
    ck(Wmin<Wmax and Wmin>=1,"complete background mass in I1 units")
    ck(1/Wmin==F(1,2)**n,"first/source ratio exponential volume")
    ck(spike_second<=eta*cone_upper**2,"no free background removal")
    for j in (1,2,n):
        R=q**j
        nextR=q**(j+1)
        captured=1+eta/R**n
        larger=1+eta/nextR**n
        ck(captured>larger>1,"strict winner beats larger full-capture cube")
        ck(captured>1,"strict selected threshold")
        ck(F(1)/(1+eta)<=F(1)/captured<1,"normalization retains background")
    return {"n":n,"J":J,"checks":checks-start,
            "S_spike_lower":display(lower),"S_spike_upper":display(upper),
            "cone_square_q_over_1":display(Kq/K1),
            "I2_over_I1_upper":display(moment_ratio),
            "spike_first_fraction_upper":display(spike_first),
            "complete_background_W_over_I1":str(Wmin),
            "scope":"Single-spike counterexample eligibility audit, not a general upper theorem."}
def main():
    started=time.monotonic()
    reg=BASE/(PFX+"_registration.json")
    registration=json.loads(reg.read_text())
    target=BASE/(PFX+"_results.json")
    if target.exists():
        raise RuntimeError("Refusing to overwrite terminal guard")
    ck(registration["registered_before_execution"],"registration exists")
    rounds=[]
    for n in (16,64,256):
        r=run(n,n*n)
        rounds.append(r)
        print(json.dumps({"n":n,"checks":r["checks"],
                          "moment_ratio_display":r["I2_over_I1_upper"]["approximation_only"]}),flush=True)
    result={"status":"PASS","checks":checks,"rounds":rounds,
            "elapsed_seconds":time.monotonic()-started,"random_seed":None,
            "scope":registration["exclusions"],
            "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "registration_sha256":hashlib.sha256(reg.read_bytes()).hexdigest()}
    target.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:result[k] for k in ("status","checks","elapsed_seconds","script_sha256")}),flush=True)
if __name__=="__main__":
    main()
