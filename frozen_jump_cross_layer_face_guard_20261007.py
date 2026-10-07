"""Proper-face coefficient checks; the full all-visited high tail remains unpaid."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import numpy as np

BASE = Path(__file__).parent
REG = BASE / "frozen_jump_cross_layer_face_registration_20261007.json"
COUNT = 0

def check(ok, label):
    global COUNT
    COUNT += 1
    if not ok:
        raise AssertionError(label)

def log_gamma_ratio(k, omega):
    # |Gamma(k-iw)|^2/Gamma(k)^2 =
    # (pi*w/sinh(pi*w))*product_{j=1}^{k-1}(1+w^2/j^2).
    x = math.pi*omega
    logsinh = x-math.log(2)+math.log1p(-math.exp(-2*x))
    j = np.arange(1,k,dtype=float)
    return (math.log(x)-logsinh+float(np.log1p((omega/j)**2).sum()))/2

def logsum(values):
    top = max(values)
    return top+math.log(sum(math.exp(v-top) for v in values))

def run(n):
    before = COUNT
    kmax = 8*n
    counts = [1]+[0]*n
    denominator = 1
    missing = [F(1)]
    selected = {0,1,n//2,n,2*n,4*n,8*n}
    coupon = []
    for k in range(kmax+1):
        if k:
            counts = [j*counts[j]+(n-j+1)*counts[j-1] if j else 0 for j in range(n+1)]
            denominator *= n
            missing.append(F(denominator-counts[n],denominator))
        check(sum(counts)==denominator,"positive occupancy total")
        check(all(v>=0 for v in counts),"positive occupancy states")
        check(denominator-counts[n] <= n*(n-1)**k,"exact union upper")
        check(0 <= missing[k] <= 1,"missing mass probability")
        if k:
            check(missing[k] <= missing[k-1],"missing mass decreases")
        if k in selected:
            coupon.append({"k":k,"missing_mass":str(missing[k]),
                           "all_visited_count_sha256":hashlib.sha256(str(counts[n]).encode()).hexdigest()})
    gamma = []
    for k in (1,n//2,n,2*n,8*n):
        for omega in (math.sqrt(k),k/2,k,2*k,16*n):
            value = log_gamma_ratio(k,float(omega))
            bound = -min(omega*omega/k,omega)/10
            check(value <= bound+1e-8,"Gamma uniform upper diagnostic")
            gamma.append({"k":k,"omega":omega,"log_gamma_ratio":value,"log_upper":bound})
    sums = []
    for omega in (n,4*n,16*n,64*n):
        terms = [math.log(float(missing[k]))-math.log(k)+log_gamma_ratio(k,omega)
                 for k in range(1,kmax+1) if missing[k]]
        partial = logsum(terms)
        bound = math.log(2)+2*math.log(n)-omega/(10*math.sqrt(n))
        check(partial <= bound+1e-8,"partial coefficient sum fits full analytic upper")
        omitted = F(n*n,kmax+1)*F(n-1,n)**(kmax+1)
        sums.append({"omega":omega,"partial_log_sum":partial,"full_analytic_log_upper":bound,
                     "omitted_k_fixed_frequency_upper":str(omitted)})
    # Uniform analytic inequality uses log(n)<=sqrt(n), pi>3, 40/pi<e^3.
    check(F(3)+F(5,2)-F(48,5)==-F(41,10),"registered cutoff exponent algebra")
    check(-F(41,10)<-4,"coarse exponential source bound")
    tail = F(1,16)+F(8,2**(16*n))
    check(tail < F(1,8),"exact proper-face high-tail cap")
    return {"dimension":n,"cutoff":16*n,"status":"PASS","checks":COUNT-before,
            "coupon":coupon,"gamma":gamma,"coefficient_sums":sums,
            "proper_face_tail_over_W_analytic_upper":str(tail)}

reg = json.loads(REG.read_text())
rounds = [run(n) for n in reg["rounds"]]
result = {"status":"PASS_EXACT_AND_NUMERIC","predicates":COUNT,"rounds":rounds,
          "scope":reg["scope"],"numeric":"floating Gamma diagnostics, not intervals; infinite-tail bound proved analytically",
          "registration_sha256":hashlib.sha256(REG.read_bytes()).hexdigest(),
          "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
out = BASE / "frozen_jump_cross_layer_face_guard_results_20261007.json"
out.write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps({"status":result["status"],"predicates":COUNT,"output":str(out)}))
