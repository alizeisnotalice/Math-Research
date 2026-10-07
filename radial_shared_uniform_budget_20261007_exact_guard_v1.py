#!/usr/bin/env python3
"""Exact original centered-hard shared-U count and L1 margin guards."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import hashlib,json,time
BASE=Path(__file__).resolve().parent
PFX="radial_shared_uniform_budget_20261007"
checks=0
def ck(ok,name):
    global checks
    if not ok:
        raise AssertionError(name)
    checks+=1
def digest_fraction(v):
    def sha_int(k):
        b=k.to_bytes(max(1,(k.bit_length()+7)//8),"big")
        return hashlib.sha256(b).hexdigest()
    return {"numerator_sha256_unsigned_bytes":sha_int(v.numerator),
            "denominator_sha256_unsigned_bytes":sha_int(v.denominator),
            "numerator_bits":v.numerator.bit_length(),
            "denominator_bits":v.denominator.bit_length(),
            "approximation_not_used_for_assertions":float(v)}
def summarize(v):
    if v.numerator.bit_length()+v.denominator.bit_length()<600:
        return str(v)
    return digest_fraction(v)
def round_guard(n,J):
    startchecks=checks
    q=1+F(1,2*J)
    L=[q**j for j in range(J+1)]
    V=[x**n for x in L]
    tau=1/(2*V[-1])
    p=1-q**(-n)
    eps=F(1,1000000*n*J*J)
    probs=[tau*V[j] for j in range(1,J+1)]
    ck(L[0]==1 and L[-1]<2,"physical window")
    ck(all(tau<V[j]**(-1) and tau*V[j]<=F(1,2) for j in range(J+1)),
       "strict threshold including A0")
    ck(p*sum(V[1:])==V[-1]-1,"annular volume telescoping")
    # S is genuinely Poisson-binomial, no independent H_j replacement.
    S=[F(1)]
    for prob in probs:
        ck(tau<=prob<=F(1,2),"actual winner probability")
        out=[F(0)]*(len(S)+1)
        for k,val in enumerate(S):
            out[k]+=(1-prob)*val
            out[k+1]+=prob*val
        S=out
    # Full count law includes B0 regardless of the shared activation event.
    rare=[p*val for val in S]
    rare[0]+=1-p
    law=[F(0)]*(J+2)
    for k,val in enumerate(rare):
        law[k]+=(1-tau)*val
        law[k+1]+=tau*val
    ck(sum(S)==1 and min(S)>=0 and sum(law)==1 and min(law)>=0,"full probability law")
    mean=sum(F(k)*val for k,val in enumerate(law))
    stop=1-law[0]
    deficit=sum(F(max(0,k-1))*val for k,val in enumerate(law))
    pair=sum(F(comb(k,2))*val for k,val in enumerate(law))
    ck(mean==F(1,2) and mean==tau*V[-1],"source normalized Lebesgue mean")
    ck(stop==tau+(1-tau)*p*(1-S[0]),"full stop includes initial winner")
    ck(mean==stop+deficit,"overlap identity")
    tails=[sum(law[m:]) for m in range(1,J+2)]
    ck(mean==sum(tails) and deficit==sum(tails[1:]),"count tails")
    ck(pair==sum(F(m-1)*tails[m-1] for m in range(2,J+2)),"pair tails")
    palm=[F(k)*val/mean for k,val in enumerate(law)]
    reciprocal=sum(palm[k]/k for k in range(1,J+2))
    ck(sum(palm)==1 and reciprocal==2*stop,"full Palm law")
    cutoff_records=[]
    for K in (1,2,max(2,J//4)):
        small=sum(palm[:K+1])
        bound=2*(1-p)*tau+2*p*K*sum(S[:K+1])
        ck(small<=bound,"Palm smallcount rareactivation bound")
        ck(small<=K*reciprocal,"Palm reciprocal consequence")
        cutoff_records.append({"K":K,"Palm_small_probability":summarize(small),
                               "analytic_upper_bound":summarize(bound)})
    # No spatial discretization: exact U norm volumes and arbitrary-y margins.
    ck(4*eps<=1-1/q and 1/q-4*eps>=0,"L1 margin eligibility")
    ulo,uhi=1/(2*q)+2*eps,F(1,2)-2*eps
    ck(ulo<uhi,"nonempty simultaneous L1 activation core")
    for j in range(1,J+1):
        ck(L[j]*ulo-eps>=L[j-1]/2+eps,"arbitrary source lower margin")
        ck(L[j]*uhi+eps<=L[j]/2-eps,"arbitrary source upper margin")
        ck(L[j]*(1/(2*q)-2*eps)+eps<=L[j-1]/2-eps,
           "inactive smaller cube is fully captured")
    ck(uhi+eps<=F(1,2)-eps,"L1 A0 retained")
    bad=(1/q+4*eps)**n-(1/q-4*eps)**n+1-(1-4*eps)**n
    ck(0<=bad<=12*n*eps,"exact bad U volume")
    ck(J*bad<=12*n*J*eps and 12*n*J*eps==F(12,1000000*J),
       "L1 complete joint-law mean error")
    activated=(1-4*eps)**n-(1/q+4*eps)**n
    ck(activated>=0 and activated<=p and p-activated<=8*n*eps,
       "common activated core probability")
    ck(tau<=stop<=tau+p and reciprocal<=2*(tau+p),"rareactivation stop range")
    # Shared-U versus independent activation law: this is the source of tails.
    ck(sum(F(k)*val for k,val in enumerate(S))==sum(probs),"conditional S mean")
    ck(sum(probs)>=J*tau,"genuine large-count conditional mean")
    records=[summarize(x) for x in law]
    canonical=json.dumps(records,sort_keys=True,separators=(",",":")).encode()
    return {"n":n,"J":J,"checks":checks-startchecks,"q":str(q),
            "tau":str(tau),"activation_probability":str(p),
            "mean":str(mean),"stop":summarize(stop),"deficit":summarize(deficit),
            "pair":summarize(pair),"Palm_reciprocal":summarize(reciprocal),
            "count_law_entries":records,
            "count_law_records_sha256":hashlib.sha256(canonical).hexdigest(),
            "cutoffs":cutoff_records,"epsilon":str(eps),
            "bad_U_probability":summarize(bad),
            "L1_TV_upper":str(12*n*eps),"L1_mean_error_upper":str(12*n*J*eps),
            "L1_scope":"Analytic uniform coupling certificate, not spatial sampled law."}
def main():
    started=time.monotonic()
    reg=BASE/(PFX+"_registration.json")
    registration=json.loads(reg.read_text())
    ck(registration["registered_before_execution"],"preregistration")
    target=BASE/(PFX+"_results.json")
    if target.exists():
        raise RuntimeError("Refusing to overwrite terminal results")
    rounds=[]
    for n,J in ((2,8),(4,16),(8,32)):
        r=round_guard(n,J)
        rounds.append(r)
        print(json.dumps({"n":n,"J":J,"checks":r["checks"],
                          "Palm_reciprocal":r["Palm_reciprocal"]}),flush=True)
    result={"status":"PASS","checks":checks,"rounds":rounds,"random_seed":None,
            "elapsed_seconds":time.monotonic()-started,"scope":registration["scope"],
            "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "registration_sha256":hashlib.sha256(reg.read_bytes()).hexdigest(),
            "actual_history_probe":"NOT_EXECUTED"}
    target.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:result[k] for k in ("status","checks","elapsed_seconds","script_sha256")}),flush=True)
if __name__=="__main__":
    main()
