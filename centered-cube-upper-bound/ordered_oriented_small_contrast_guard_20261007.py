#!/usr/bin/env python3
"""Pure Fraction checks of analytic scalar envelopes, not of actual sources."""
from pathlib import Path
from fractions import Fraction as F
import hashlib, json
HERE=Path(__file__).resolve().parent
REG=HERE/"ordered_oriented_small_contrast_registration_20261007.json"
reg=json.loads(REG.read_text())
assert reg["status"]=="registered_before_execution"
counts={}
min_margins={}
def check(round,label,b):
    if not b: raise AssertionError((round,label))
    counts[round]=counts.get(round,0)+1
node=[F(0),F(1,16),F(1,4),F(1,2),F(1),F(2),F(8),F(128)]
for N in reg["rounds"]:
    for delta in map(F,reg["delta"]):
        for j in range(N+1):
            r=1+delta*F(j,N)
            Q=(r-1)**2/r
            check(N,"near reciprocal contrast",r-1/r<=2*delta)
            check(N,"near cubic remainder envelope",(r-1/r)/12<=delta/6)
            check(N,"scaled entropy remainder",(r-1/r)/6<=delta/3)
            for rfar in (1+delta*(1+F(j,N)),1+delta*2**(j%13)):
                Qfar=(rfar-1)**2/rfar
                check(N,"far logarithm rational envelope",
                      rfar-1<=((1+delta)/delta)*Qfar)
                for k in range(17):
                    dt=F(k,8)
                    # D= lambda Q for a pair of reverse edges of unit measure.
                    check(N,"double-edge far coefficient",
                          dt*(rfar-1)<=2*(1+delta)/delta*Qfar)
        for lam in map(F,reg["lambda"]):
            for u in node:
                for v in node:
                    # Node values already in units of lambda.
                    mu=max(F(1),u)+F(1,N)
                    mv=max(F(1),v)+F(1,N)
                    for factor in (F(1),F(2)):
                        M=lam*mu*factor
                        Mp=lam*mv*factor
                        for gate in (0,1):
                            for gatep in (0,1):
                                H=gate*lam/M
                                Hp=gatep*lam/Mp
                                U,V=lam*u,lam*v
                                t=H*(1+U/lam)
                                tp=Hp*(1+V/lam)
                                check(N,"bounded cap",0<=t<=2 and 0<=tp<=2)
                                check(N,"H response cap",
                                      H<=2*lam/(lam+U) and Hp<=2*lam/(lam+V))
                                check(N,"exact product difference",
                                      tp-t==(H+Hp)*(V-U)/(2*lam)
                                      +(1+(U+V)/(2*lam))*(Hp-H))
                                if H==Hp:
                                    check(N,"same-H signed free",(tp-t)*(V-U)>=0)
                                if gate==gatep==1 and V>U:
                                    check(N,"both-alive maximum order",
                                          (t>tp)==(Mp/M>(lam+V)/(lam+U)))
        # Sharp elementary endpoint log(2)<3/4 via exp lower polynomial.
        x=F(3,4)
        check(N,"exp(3/4)>2",1+x+x*x/2+x*x*x/6>2)
res={"status":"PASS_EXACT_FRACTION_ENVELOPES","counts_by_round":counts,
     "total_checks":sum(counts.values()),
     "registration_sha256":hashlib.sha256(REG.read_bytes()).hexdigest(),
     "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     "notes":["Log envelopes themselves have analytic proofs in Markdown",
              "Scalar edge choices do not assert actual evolution realizability",
              "Only a paid high-contrast branch and an unpaid signed near interface",
              "No actual source simulation or cube fee"]}
out=HERE/"ordered_oriented_small_contrast_results_20261007.json"
out.write_text(json.dumps(res,indent=2)+"\n")
print(json.dumps(res,indent=2))
