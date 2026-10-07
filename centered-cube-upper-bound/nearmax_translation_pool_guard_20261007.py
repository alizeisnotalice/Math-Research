#!/usr/bin/env python3
"""Fresh exact complete-halo cell guard for whole-source translation pooling.
No old oracle/module import. Rational logs use a newly implemented 100-term
atanh series with binary range reduction and rigorous positive tail.
"""
from fractions import Fraction as F
from itertools import product
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json, gzip, time

HERE = Path(__file__).resolve().parent
P = "nearmax_translation_pool_guard_20261007"
J = [F(1), F(2)]
CHECKS = defaultdict(int)


def check(round_id, name, valid):
    CHECKS[round_id] += 1
    if not valid:
        raise AssertionError(name)


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def frac(value):
    value = F(value)
    return str(value.numerator) + "/" + str(value.denominator)


def cert(value):
    value = F(value)
    scale = 1 << 256
    floor = value.numerator * scale // value.denominator
    # Integer byte hash avoids large integer-to-decimal serialization.
    num, den = value.numerator, value.denominator
    payload = (b"-" if num < 0 else b"+") + abs(num).to_bytes(max(1, (abs(num).bit_length()+7)//8), "big")
    payload += b"/" + den.to_bytes(max(1, (den.bit_length()+7)//8), "big")
    out = {"fraction_byte_sha256": sha256(payload).hexdigest(),
           "dyadic256_lower": frac(F(floor, scale)),
           "dyadic256_upper": frac(F(floor + 1, scale)),
           "float_diagnostic_only": float(value)}
    if num.bit_length() + den.bit_length() < 400:
        out["exact_fraction"] = frac(value)
    return out


def atanh_log_bounds(r):
    z = (r - 1) / (r + 1)
    assert 0 <= z <= F(1, 3)
    lower = F(0)
    power = z
    for j in range(100):
        lower += 2 * power / (2*j + 1)
        power *= z*z
    remainder = 2 * power / (201*(1-z*z))
    return lower, lower + remainder


LOG2 = atanh_log_bounds(F(2))


@lru_cache(None)
def log_bounds(q):
    assert q > 0
    if q == 1:
        return F(0), F(0)
    r, k = q, 0
    while r >= 2:
        r /= 2
        k += 1
    while r < 1:
        r *= 2
        k -= 1
    lo, hi = atanh_log_bounds(r)
    if k >= 0:
        return lo+k*LOG2[0], hi+k*LOG2[1]
    return lo+k*LOG2[1], hi+k*LOG2[0]


def add_log(mapping, argument, coefficient):
    argument, coefficient = F(argument), F(coefficient)
    if argument == 1 or coefficient == 0:
        return
    mapping[argument] += coefficient
    if mapping[argument] == 0:
        del mapping[argument]


def plus_logs(*weighted_maps):
    out = defaultdict(F)
    for weight, mapping in weighted_maps:
        for q, c in mapping.items():
            add_log(out, q, weight*c)
    return out


def evaluate(mapping):
    low = high = F(0)
    for q, c in mapping.items():
        lo, hi = log_bounds(q)
        if c >= 0:
            low += c*lo
            high += c*hi
        else:
            low += c*hi
            high += c*lo
    return low, high


def describe(mapping):
    low, high = evaluate(mapping)
    return {
        "merged_rational_log_coefficients": [{"argument": frac(q), "coefficient": frac(c)}
                                             for q, c in sorted(mapping.items())],
        "rigorous_interval": {"lower": cert(low), "upper": cert(high)},
        "sign": "POSITIVE" if low > 0 else "NEGATIVE" if high < 0
                else "EXACT_ZERO" if not mapping else "UNRESOLVED_INTERVAL_CROSSES_ZERO",
    }


def response(source, x, R):
    return sum((m for y, m in source if all(2*abs(a-b) <= R for a, b in zip(x, y))), F(0)) / R**len(x)


def true_winner(source, x):
    vals = [response(source, x, R) for R in J]
    maximum = max(vals)
    return maximum, J[vals.index(maximum)], vals


def run_round(n):
    tau = F(3, 8*2**n)
    z = tuple(F(1, 4) for _ in range(n))
    source = [(tuple(F(0) for _ in range(n)), F(1, 3)),
              (tuple(F(3, 2) for _ in range(n)), F(2, 3))]
    plus = [(tuple(y_i+z_i for y_i, z_i in zip(y, z)), m) for y, m in source]
    minus = [(tuple(y_i-z_i for y_i, z_i in zip(y, z)), m) for y, m in source]
    pooled = [(y, m/2) for y, m in plus+minus]
    check(n, "complete original mass one", sum(m for _,m in source) == 1)
    check(n, "complete translated mass plus", sum(m for _,m in plus) == 1)
    check(n, "complete translated mass minus", sum(m for _,m in minus) == 1)
    check(n, "complete pooled mass one", sum(m for _,m in pooled) == 1)
    check(n, "pooled source equals half complete translations", pooled == [(y, m/2) for y,m in plus+minus])
    axes = []
    for i in range(n):
        cuts = sorted({y[i] + sign*R/2 for y,_ in source+plus+minus for R in J for sign in [-1,1]})
        axes.append(list(zip(cuts[:-1], cuts[1:])))
    phi = {name: defaultdict(F) for name in ["original","plus","minus","pooled"]}
    gain, regret, slack_positive = defaultdict(F), defaultdict(F), defaultdict(F)
    symdiff = intersection = f0volume = nonzero_regret_volume = F(0)
    integrals = {name:[F(0),F(0)] for name in phi}
    cells = []
    for bounds in product(*axes):
        x = tuple((left+right)/2 for left,right in bounds)
        volume = F(1)
        for left,right in bounds:
            volume *= right-left
        O, RO, vo = true_winner(source, x)
        A, RP, vp = true_winner(plus, x)
        B, RM, vm = true_winner(minus, x)
        MN, RN, vn = true_winner(pooled, x)
        d1 = B-response(minus, x, RP)
        d2 = A-response(plus, x, RM)
        d = min(d1,d2)
        beta = d/(A+B) if A+B else F(0)
        T0 = (A+B)/2
        eplus, eminus, f0 = A>tau, B>tau, T0>tau
        cross = eplus != eminus
        both = eplus and eminus
        check(n, "cell positive exact volume", volume > 0)
        check(n, "actual first cross regret nonnegative", d1 >= 0)
        check(n, "actual second cross regret nonnegative", d2 >= 0)
        check(n, "minimum regret <=min(A,B)", 0 <= d <= min(A,B))
        check(n, "beta in exact [0,1/2]", 0 <= beta <= F(1,2))
        check(n, "real pooled winner above two-winner lower", MN >= T0*(1-beta))
        check(n, "real pooled maximum below average maxima", MN <= T0)
        check(n, "pooled fixed responses average whole translated sources", vn == [(a+b)/2 for a,b in zip(vp,vm)])
        check(n, "winner tie rule original smallest R", RO == min(R for R,val in zip(J,vo) if val==O))
        check(n, "winner tie rule plus smallest R", RP == min(R for R,val in zip(J,vp) if val==A))
        check(n, "winner tie rule minus smallest R", RM == min(R for R,val in zip(J,vm) if val==B))
        check(n, "winner tie rule pooled smallest R", RN == min(R for R,val in zip(J,vn) if val==MN))
        for name, maxval, vals in [("original",O,vo),("plus",A,vp),("minus",B,vm),("pooled",MN,vn)]:
            if maxval > tau:
                add_log(phi[name], maxval/tau, tau*volume)
            for idx,val in enumerate(vals):
                integrals[name][idx] += volume*val
        amp_argument = (A+B)**2/(4*A*B) if both else F(1)
        check(n, "amplitude logcosh argument >=1", amp_argument >= 1)
        if both:
            intersection += volume
            add_log(gain, amp_argument, tau*volume/2)
        if cross:
            symdiff += volume
        if f0:
            f0volume += volume
            add_log(regret, 1/(1-beta), tau*volume)
            if beta > 0:
                nonzero_regret_volume += volume
        # Twice each cell's exact log residual, regrouped into one rational
        # argument. The positivity certificate avoids artificial cancellation.
        left_ratio = max(F(1),MN/tau)**2/(max(F(1),A/tau)*max(F(1),B/tau))
        rhs_ratio = amp_argument
        if cross:
            rhs_ratio /= 2
        if f0:
            rhs_ratio *= (1-beta)**2
        slack_argument = left_ratio/rhs_ratio
        check(n, "full pointwise pooling residual exponential >=1", slack_argument >= 1)
        add_log(slack_positive, slack_argument, tau*volume/2)
        cells.append({
            "bounds":[[frac(l),frac(r)] for l,r in bounds], "volume":frac(volume),
            "M_original":frac(O), "M_plus":frac(A), "M_minus":frac(B), "M_pooled":frac(MN),
            "winner_original":frac(RO), "winner_plus":frac(RP), "winner_minus":frac(RM),
            "winner_pooled":frac(RN), "d1":frac(d1), "d2":frac(d2), "beta":frac(beta),
            "strict_Eplus":eplus, "strict_Eminus":eminus, "strict_F0":f0,
            "slack_argument":frac(slack_argument),
        })
    check(n, "complete translated phi plus exact log map matches original", phi["plus"] == phi["original"])
    check(n, "complete translated phi minus exact log map matches original", phi["minus"] == phi["original"])
    for name, vals in integrals.items():
        for idx,value in enumerate(vals):
            check(n, name+" complete fixed-kernel source-once integral "+str(idx), value==1)
    delta = plus_logs((F(1),phi["pooled"]),(-F(1),phi["original"]))
    crossing = defaultdict(F)
    add_log(crossing, F(2), tau*symdiff/2)
    residual = plus_logs((F(1),delta),(-F(1),gain),(F(1),crossing),(F(1),regret))
    residual_interval = evaluate(residual)
    positive_interval = evaluate(slack_positive)
    check(n, "two exact log representations interval overlap", max(residual_interval[0],positive_interval[0]) <= min(residual_interval[1],positive_interval[1]))
    # Sign is certified from the independently positive rational arguments
    # even in an exactly zero case; undecided direct cancellation is retained.
    check(n, "positive residual representation all nonnegative terms", all(q>1 and c>0 for q,c in slack_positive.items()))
    check(n, "positive residual interval lower nonnegative", positive_interval[0] >= 0)
    raw = json.dumps(cells,ensure_ascii=False,separators=(",",":")).encode()
    path = HERE/(P+"_n"+str(n)+"_cells.json.gz")
    path.write_bytes(gzip.compress(raw,mtime=0))
    return {
        "n":n, "tau":frac(tau), "z":[frac(z_i) for z_i in z],
        "source":[{"y":[frac(a) for a in y],"mass":frac(m)} for y,m in source],
        "scales":[frac(R) for R in J], "cell_count":len(cells),
        "axis_interval_counts":[len(axis) for axis in axes], "checks":CHECKS[n],
        "strict_intersection_volume":frac(intersection), "strict_symmetric_difference_volume":frac(symdiff),
        "strict_F0_volume":frac(f0volume), "nonzero_regret_F0_volume":frac(nonzero_regret_volume),
        "phi_original":describe(phi["original"]), "phi_pooled":describe(phi["pooled"]),
        "phi_pooled_minus_original":describe(delta), "amplitude_gain":describe(gain),
        "threshold_crossing_cost":describe(crossing), "cross_winner_regret_cost":describe(regret),
        "direct_integrated_residual":describe(residual),
        "positive_recombined_residual":describe(slack_positive),
        "cells_file":path.name, "cells_sha256":sha(path),
        "complete_fixed_kernel_integrals":{name:[frac(x) for x in vals] for name,vals in integrals.items()},
        "certificate":"Full inequality holds from exact cellwise residual_argument>=1; totals also evaluated by rigorous log intervals.",
    }


def main():
    started = time.perf_counter()
    started_utc = datetime.now(timezone.utc).isoformat()
    registration_path = HERE/(P+"_registration.json")
    reg = json.loads(registration_path.read_text())
    check(0,"frozen script SHA",sha(Path(__file__))==reg["script_sha256"])
    check(0,"frozen analytic contract SHA",sha(HERE/"nearmax_smoothing_obstacle_20261007.md")==reg["proof_sha256"])
    check(0,"registered rounds match",[1,2,3]==reg["dimensions"])
    rounds = [run_round(n) for n in reg["dimensions"]]
    result = {"status":"PASS","seed":None,"dimensions":reg["dimensions"],
              "checks_total":sum(CHECKS.values()),"rounds":rounds,
              "script_sha256":sha(Path(__file__)),"registration_sha256":sha(registration_path),
              "proof_sha256":reg["proof_sha256"],"log_method":reg["log_method"],"scope":reg["scope"]}
    result_path = HERE/(P+"_results.json")
    result_path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    receipt = {"status":"PASS","started_utc":started_utc,"completed_utc":datetime.now(timezone.utc).isoformat(),
               "elapsed_seconds":time.perf_counter()-started,"checks_total":sum(CHECKS.values()),
               "round_check_counts":{str(n):CHECKS[n] for n in reg["dimensions"]},
               "registration_checks":CHECKS[0],"script_sha256":sha(Path(__file__)),
               "registration_sha256":sha(registration_path),"proof_sha256":reg["proof_sha256"],
               "results_sha256":sha(result_path),"seed":None,"execution_count":1,
               "old_oracle_invocations":0,"scope":reg["scope"]}
    (HERE/(P+"_receipt.json")).write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({"status":"PASS","checks_total":sum(CHECKS.values()),
                     "rounds":[{"n":r["n"],"cells":r["cell_count"],"checks":r["checks"]} for r in rounds],
                     "elapsed_seconds":receipt["elapsed_seconds"]}))


if __name__=="__main__":
    main()

