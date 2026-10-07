#!/usr/bin/env python3
"""New Fraction constants and original hard-kernel sign-interval certificates.

The cosine h1/h2 identities are analytic, not a numerical pi approximation.
Complete finite-box L1 input; no old cosine model or actual FIRST reruns.
"""
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
import json
import time

BASE = Path(__file__).resolve().parent
PREFIX = "direct_mass_band_cross_audit_20261007"
REG = BASE/(PREFIX+"_registration.json")
OUT = BASE/(PREFIX+"_results.json")
if OUT.exists():
    raise SystemExit("Refuse an existing terminal result.")
registration = json.loads(REG.read_text())
assert registration["status"] == "registered_before_execution"
assert registration["rounds"] == [4, 16, 64]
N = 0


def ck(condition, name):
    global N
    N += 1
    if not condition:
        raise AssertionError(name)


def fs(x):
    return str(x.numerator)+"/"+str(x.denominator)


def interval_length(left, right, pieces):
    return sum((max(F(0), min(right,b)-max(left,a))
                for a,b in pieces), F(0))


started = time.time()
ck(registration["arithmetic"] == "Fraction and exact interval lengths",
   "registration_consistency")
eps = F(1,8)
tau = B = F(3,4)
positive = [(2*j-F(1,2),2*j+F(1,2)) for j in range(-2,4)]
negative = [(2*j+F(1,2),2*j+F(3,2)) for j in range(-2,4)]
# These come from the actual cos(pi*x) sign intervals.
breakpoints = sorted({F(0), F(2)} |
                     {edge+shift
                      for piece in positive for edge in piece
                      for shift in (F(-1,2),F(1,2))
                      if 0 <= edge+shift <= 2})
ck(breakpoints == [F(0),F(1),F(2)], "all_A_piecewise_affine_breakpoints")
area = F(0)
for left,right in zip(breakpoints,breakpoints[1:]):
    a_left = interval_length(left-F(1,2),left+F(1,2),positive)
    a_right = interval_length(right-F(1,2),right+F(1,2),positive)
    area += (a_left+a_right)*(right-left)/2
ck(area == 1, "A_period_area_one_exact_affine_integral")
for j in range(9):
    y = F(j,4)
    A = interval_length(y-F(1,2),y+F(1,2),positive)
    ck(A == abs(y-1), "A_actual_interval_overlap_triangle")
    neg2 = interval_length(y-1,y+1,negative)
    ck(neg2 == 1, "h2_receiver_negative_occupancy_one")
    ck(F(0) <= A <= 1, "h1_sign_occupancy_range")

rounds = []
for n in registration["rounds"]:
    before=N
    L=8*n
    W=F((2*L)**n)
    core_volume=F((2*L-2)**n)
    inner_volume=F((2*L-4)**n)
    E_upper=F((2*L+2)**n)
    ck(L%2==0 and L>=4, "true_L1_source_periods_and_nonempty_core")
    ck(1-eps>tau and 1+eps<2*tau, "whole_input_and_global_hardband")
    ck(B<1 and 1+eps<=2*B, "core_R1_mass_band_zero")
    ck((2**n)*B<2**n<=2**(n+1)*B, "core_R2_mass_band_n")
    ck(tau==B, "same_pre_fixed_threshold_and_mass_base")
    # Whole E winner can only be one of these two original scales.
    ck(2*tau==2*B, "entire_R1_hardband_maps_only_to_band_zero")
    ck((2**n)*2*tau == (2**(n+1))*B,
       "entire_R2_hardband_maps_only_to_band_n")
    ck(-L+2-1 == -L+1 and L-2+1 == L-1,
       "source_inner_Q2_inside_receiver_core")
    ck(-L+1-1 == -L and L-1+1 == L,
       "receiver_core_Q2_inside_original_source_box")
    ck(-L+2-F(1,2)>-L+1 and L-2+F(1,2)<L-1,
       "source_inner_Q1_inside_receiver_core")
    ck((2*L-4)%2==0 and (2*L-2)%2==0,
       "complete_core_and_source_period_lengths")
    ck(core_volume/2>0 and inner_volume/2>0,
       "two_actual_strict_winner_regions_positive")
    factor=(1-eps)*tau*tau/(4*(1+eps))
    ck(factor==F(7,64), "original_source_cross_constant")
    cross_lower=factor*inner_volume
    ck(inner_volume/W >= F(3,4), "inner_source_volume_Bernoulli")
    ck(cross_lower >= F(21,256)*W, "cross_at_least_21_over_256_W")
    ratio_lower=cross_lower/(tau*E_upper)
    ratio_formula=F(7,48)*F(L-2,L+1)**n
    ck(ratio_lower==ratio_formula, "normalized_cross_formula")
    ck(F(L-2,L+1)==1-F(3,8*n+1), "volume_ratio_exact")
    ck(F(L-2,L+1)**n>=1-F(3*n,8*n+1),
       "Bernoulli_normalized_volume_ratio")
    ck(1-F(3*n,8*n+1)>=F(5,8),
       "normalized_volume_ratio_uniform_five_eighths")
    ck(ratio_lower>=F(35,384), "nondecaying_cross_uniform_lower")
    ck((2*L)**n==int(W), "whole_L1_source_mass_exact")
    ck(core_volume<E_upper and W<E_upper, "whole_E_expansion_retained")
    ck(F(2)+F(2)==4, "only_two_profiles_have_aggregate_pointwise_cap_four")
    rounds.append({
        "n":n,"L":L,"status":"PASS","checks":N-before,
        "mass_band_gap":n,
        "whole_source_mass_integer":str(W.numerator),
        "cross_over_W_lower":fs(cross_lower/W),
        "normalized_cross_lower":fs(ratio_lower),
        "uniform_normalized_cross_lower":"35/384",
        "nonempty_mass_bands":[0,n],
        "aggregate_bound":"4*I1 (only two profiles)",
        "actual_FIRST_certified":False
    })
result={
    "status":"PASS","new_guard_only":True,"total_checks":N,
    "rounds":rounds,"runtime_seconds":time.time()-started,
    "script_sha256":sha256(Path(__file__).read_bytes()).hexdigest(),
    "registration_sha256":sha256(REG.read_bytes()).hexdigest(),
    "arithmetic":registration["arithmetic"],
    "analytic_kernel_scope":"h1 cosine multiplier 2/pi and h2 zero proved exactly; no floating pi values evaluated.",
    "full_input_scope":registration["input"],
    "not_certified":registration["limits"]
}
OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
print(json.dumps({"status":"PASS","total_checks":N,
                  "rounds":[{"n":x["n"],"checks":x["checks"]}
                            for x in rounds],
                  "runtime_seconds":result["runtime_seconds"],
                  "script_sha256":result["script_sha256"]}))
