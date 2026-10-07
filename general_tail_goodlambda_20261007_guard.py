"""Registered scalar guards; NOT an original spatial/gate or order certificate."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
from datetime import datetime, timezone

HERE = Path(__file__).resolve().parent
PREFIX = "general_tail_goodlambda_20261007"

def log2_bounds(terms):
    z = F(1, 3)
    lower = 2 * sum((z ** (2*k+1) / (2*k+1) for k in range(terms)), F(0))
    upper = lower + 2*z**(2*terms+1)/((2*terms+1)*(1-z*z))
    return lower, upper

def interval_mul(k, interval):
    lo, hi = interval
    return (k*lo, k*hi) if k >= 0 else (k*hi, k*lo)

rounds = []
previous = None
total = 0
for n, s, terms in ((16, 4, 16), (64, 8, 24), (256, 16, 32)):
    checks = []
    def check(name, condition):
        assert condition, (n, name)
        checks.append(name)
    li = log2_bounds(terms)
    lo, hi = li
    check("log2_positive_interval", F(2, 3) < lo < hi < F(7, 10))
    check("log_interval_width", hi-lo < F(1, 10**(terms//2)))
    if previous is not None:
        check("log_interval_nested", previous[0] < lo and hi < previous[1])
    previous = li
    exponent = s + 2
    T = 2**exponent
    Hlo, Hhi = interval_mul(2*exponent, li)
    check("H_positive", Hlo > 0)
    check("H_below_strong_envelope", Hhi < 1+n*lo)
    # ell(r) <= log(T)+r/T.  All positive tested r are powers of 2.
    points = (("zero", F(0), None), ("half", F(1,2), -1),
              ("one", F(1), 0), ("T_half", F(T,2), exponent-1),
              ("T", F(T), exponent), ("T_double", F(2*T), exponent+1),
              ("T_square", F(T*T), 2*exponent))
    for name, r, e in points:
        if r == T:
            check("tangent_exact_at_T", F(1)+F(T,T) == 2)
            # Both sides are exactly 1+log(T).
            continue
        ell_hi = r if r <= 1 else 1+interval_mul(e,li)[1]
        rhs_lo = exponent*lo+r/T
        check("tangent_"+name, ell_hi < rhs_lo)
    eta, delta0, delta1, nu = F(1,4), F(1,8), F(1,4), F(1,4)
    fee = 1/(eta*T)
    check("full_fee_exact", fee == F(4,T))
    check("partial_absorption", delta0+fee <= F(1,4))
    retention = 1-nu-F(1,s)
    check("conditional_total_absorption", delta0+fee+delta1 < retention)
    # Scalar component classification only; these are not spatial inputs.
    w = F(4)
    fixtures = (("low_endpoint", F(1), F(8), "low"),
                ("full_endpoint", F(1), F(4), "full"),
                ("fragment", F(1,2), F(2), "fragment"))
    for name, c, m, expected in fixtures:
        b = c/m
        actual = "low" if b <= delta0 else ("full" if c >= eta*w else "fragment")
        check(name+"_classification", actual == expected)
        if actual == "fragment":
            check("fragment_mass_cap", m < eta*w/delta0)
    # The contradictory endpoints are exact: threshold m>tau*a^n versus m<cap.
    cap, threshold = eta*w/delta0, F(8)
    check("empty_fragment_cap_at_threshold", cap == threshold)
    row = {"n": n, "sqrt_n": s, "T": T, "log_terms": terms,
           "H_interval": [str(Hlo),str(Hhi)], "full_fee": str(fee),
           "partial_coefficient": str(delta0+fee),
           "conditional_total_coefficient": str(delta0+fee+delta1),
           "original_retention": str(retention), "checks": checks,
           "check_count": len(checks), "status": "PASS"}
    rounds.append(row)
    total += len(checks)

results = {"scope": "rational scalar component guards only; no original spatial/gate or order evidence",
           "fragment_contract_proved": False, "rounds": rounds,
           "assertions": total, "status": "PASS"}
(HERE/(PREFIX+"_results.json")).write_text(json.dumps(results,ensure_ascii=False,indent=2)+"\n")
receipt = {"executed_at_utc": datetime.now(timezone.utc).isoformat(),
           "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           "registration_sha256": hashlib.sha256((HERE/(PREFIX+"_registration.md")).read_bytes()).hexdigest(),
           "assertions": total, "status": "PASS", "executions_authorized": 1,
           "scope": results["scope"]}
(HERE/(PREFIX+"_receipt.json")).write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+"\n")
print(json.dumps({"status":"PASS","rounds":len(rounds),"assertions":total,
                  "full_fees":[r["full_fee"] for r in rounds]}))
