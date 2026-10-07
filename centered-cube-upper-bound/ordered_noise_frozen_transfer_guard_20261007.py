#!/usr/bin/env python3
"""Registered local transfer guard on the original B/m/G_c symbol.

Fraction checks are exact. Decimal identities are high precision diagnostics,
not intervals; uniform operator/source statements are proved analytically.
"""
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import math

HERE=Path(__file__).resolve().parent
REG=HERE/"ordered_noise_frozen_transfer_registration_20261007.json"
OUT=HERE/"ordered_noise_frozen_transfer_guard_results_20261007.json"
reg=json.loads(REG.read_text())
predicates=0

def check(ok,name):
    global predicates
    predicates+=1
    if not ok:
        raise AssertionError(name)

def dec(f):
    f=Fraction(f)
    return Decimal(f.numerator)/Decimal(f.denominator)

rounds=[]
with localcontext() as ctx:
    ctx.prec=reg["decimal_precision"]
    tol=Decimal("1e-70")
    for n in reg["rounds"]:
        m=math.isqrt(n)
        if m*m<n:
            m+=1
        r0=Fraction(1,2*m)
        check(m*m>=n and (m-1)*(m-1)<n,"ceil sqrt")
        check(n*r0*r0<=Fraction(1,4),"window size budget")
        check(Fraction(n)*r0*r0/(2*(1-r0)**2)<=Fraction(1,2),
              "log dominance < 1/2")
        mask_records=[]
        symbols=[]
        budgets=[]
        for cstr in reg["c0"]:
            c0=Fraction(cstr)
            check(Fraction(1,2)<=c0<=1,"actual early c range")
            for r in (Fraction(0),r0/2,r0):
                a=r/(1-r)
                cd=dec(c0)
                rd=dec(r)
                ad=dec(a)
                constant=(Decimal(n)*((1-rd).ln()+ad)).exp()
                logupper=Fraction(n)*r*r/(2*(1-r)**2)
                check(constant<=dec(logupper).exp()+tol,"constant direction")
                check(constant<2,"uniform local cost")
                if r:
                    for j in range(n+1):
                        bern=(1-rd)**(n-j)*rd**j
                        poisson=(-Decimal(n)*ad).exp()*ad**j
                        check(abs(bern-constant*poisson)<=tol*max(bern,Decimal("1e-200")),
                              "all mask no-repeat coefficient")
                    no_repeat=(-Decimal(n)*ad).exp()*(1+ad)**n
                    check(abs(no_repeat*constant-1)<tol,"whole no-repeat mass")
                else:
                    check(constant==1,"r0 holding identity")
                mask_records.append({"c0":cstr,"r":str(r),"clock_over_c0":str(a),
                    "constant":str(constant),"log_constant_upper":str(logupper)})
                cs=c0*(1-r)
                delta=c0*r
                for lstr in reg["symbol_lambda"]:
                    v=dec(Fraction(lstr))
                    log=(1+v).ln()
                    b=v/log-1
                    w=log/v
                    g=1/(1+cd*b)
                    nold=cd+(1-cd)*w
                    ns=dec(cs)+(1-dec(cs))*w
                    kt=1-rd+rd*g
                    check(abs(ns/nold-kt)<tol,"original cocycle multiplier")
                    check(0<kt<=1 and 0<g<1,"original positivity multiplier")
                    psi0=b/(1+cd*b)
                    psis=b/(1+dec(cs)*b)
                    forcing=(psis-psi0)*kt
                    rhs=dec(delta)/(cd*cd)*(1-g)**2
                    check(abs(forcing-rhs)<tol,"free forcing symbol identity")
                    check(abs(g*nold-w)<tol,"actual replacement identity")
                    symbols.append({"c0":cstr,"r":str(r),"lambda":lstr,
                        "transition_residual":str(abs(ns/nold-kt)),
                        "forcing_residual":str(abs(forcing-rhs))})
                free=2*n*r*r
                killed_upper=2*n*r*r/(1-r)
                killed_actual=4*Decimal(n)*(-rd-(1-rd).ln())
                check(killed_actual<=dec(killed_upper)+tol,"killed integrated variation")
                check(free<=Fraction(1,2),"free local source mass")
                check(killed_upper<=1,"killed local source mass")
                budgets.append({"c0":cstr,"r":str(r),
                    "free_forcing_mass_over_W0":str(free),
                    "killed_forcing_mass_upper_over_W0":str(killed_upper),
                    "killed_exact_log_budget_decimal":str(killed_actual)})
        # One piecewise approximation on the actual early interval. No fresh
        # input N_t nu is assigned a new W in this spacetime occupation bound.
        c=Fraction(1)
        partition=[c]
        while c>Fraction(1,2):
            nextc=max(Fraction(1,2),c*(1-r0))
            relative=(c-nextc)/c
            check(0<relative<=r0,"early window coverage")
            check(Fraction(1,2)<=nextc<=c<=1,"full early c interval")
            max_generator_variation=4*n*relative/nextc
            check(max_generator_variation<=8*n*r0,"single live occupation coefficient")
            c=nextc
            partition.append(c)
        total_occupation_mass_upper=4*n*r0
        check(total_occupation_mass_upper==Fraction(2*n,m),"single early time integral")
        check((total_occupation_mass_upper/2)**2<=n,"O sqrt n forcing variation only")
        rounds.append({"n":n,"ceil_sqrt":m,"relative_window_max":str(r0),
            "mask_records":mask_records,"symbols":symbols,"budgets":budgets,
            "early_c_partition":[str(x) for x in partition],
            "early_forcing_spacetime_TV_over_W_upper":str(total_occupation_mass_upper),
            "not_an_output_error_bound":True})

result={"status":"PASS_EXACT_AND_DECIMAL","predicates":predicates,
    "registration_sha256":hashlib.sha256(REG.read_bytes()).hexdigest(),
    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "decimal_scope":"80 digit diagnostics of original logarithmic symbol; not interval certificates.",
    "actual_scope":"early initial time <=1/2; no receiver/source reassignment, full FIRST or CPGP certification.",
    "rounds":rounds}
OUT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
print(result["status"],predicates,"predicates")
