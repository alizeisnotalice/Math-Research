#!/usr/bin/env python3
"""New exact original-hard radial fragmentation components, no old reruns."""
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
import json
import time

ROOT = Path(__file__).resolve().parent
PREFIX = "scale_aware_capture_rearrangement_20261007"
REG = ROOT/(PREFIX+"_registration.json")
OUT = ROOT/(PREFIX+"_results.json")
if OUT.exists():
    raise SystemExit("Refuse to overwrite a terminal result.")
registration = json.loads(REG.read_text())
assert registration["status"] == "registered_before_execution"
assert registration["rounds"] == [4, 8, 12]
count = 0


def ck(v, label):
    global count
    count += 1
    if not v:
        raise AssertionError(label)


def fs(x):
    return str(x.numerator)+"/"+str(x.denominator)


def cert(x):
    raw = fs(x)
    return {"fraction_sha256": sha256(raw.encode()).hexdigest(),
            "numerator_bits": x.numerator.bit_length(),
            "denominator_bits": x.denominator.bit_length(),
            "preview_not_certificate": float(x)}


started = time.time()
ck(registration["arithmetic"] == "Fraction only", "registered_arithmetic")
rounds = []
for n in registration["rounds"]:
    before = count
    J = 2**n
    delta = F(1, 100*n*J*J)
    helper_mass = delta/(J+1)
    ts = [1+F(j, J) for j in range(J+1)]
    ck(helper_mass*(J+1) == delta, "entire_helper_source_mass")
    ck(1+delta+6*4**(n-1) > 0, "entire_background_source_mass")
    ck(delta < F(1, J), "precentral_signal_margin")
    ck((1+2*delta)**n > 1+delta, "later_arrival_uniform_signal_gap")
    ck(1-(1+delta)**(-n) <= n*delta, "power_loss_bound")
    sum_vol = F(0)
    largest_volume = F(0)
    min_width = None
    for j in range(J):
        left = ts[j]
        right_untrimmed = ts[j+1]
        right = right_untrimmed/(1+delta)
        width = right-left
        volume = (right**n-left**n)/(2*n)
        ideal = (right_untrimmed**n-left**n)/(2*n)
        m = 1+(j+1)*helper_mass
        coeff = volume/(1+m)
        ck(F(1) <= left < right < right_untrimmed <= F(2),
           "retained_radial_interval_positive")
        ck(width > 0 and volume > 0, "true_cone_positive_volume")
        ck(volume <= ideal <= F(1, 4), "each_block_low_volume")
        ck(ideal <= F(2**(n-2), J), "mean_value_volume_bound")
        ck(1 < m <= 1+delta < 2, "actual_prefix_capture_mass")
        ck(coeff >= volume/3, "actual_scale_aware_coefficient_lower")
        ck(ideal-volume <= delta*right_untrimmed**n/2,
           "exact_volume_loss")
        # New-source arrival at the upper retained boundary is already
        # (1+2 delta) times R0; the open interval gives strict inequality.
        ck(2*right_untrimmed-right == (1+2*delta)*right,
           "new_helper_arrival_ratio_at_retained_boundary")
        ck((1+delta)/(2*right_untrimmed-right)**n < 1/right**n,
           "true_continuous_future_candidate_lower_signal")
        # Prefix t_0..t_j is captured at R0, later points are not.
        ck(ts[j] <= left and ts[j+1] > right,
           "true_prefix_membership")
        # Complete central-source cube arrival is inside original scale range.
        mid = (left+right)/2
        ck(F(1) < mid < F(2), "central_winner_scale_domain")
        ck(1/mid**n > delta, "central_response_dominates_all_prearrival")
        sum_vol += volume
        largest_volume = max(largest_volume, volume)
        min_width = width if min_width is None else min(min_width, width)

    full_cone = F(J-1, 2*n)
    ck(sum_vol >= full_cone-F(1, 200*n), "summed_true_cone_loss_bound")
    ck(sum_vol >= F(J-1, 4*n), "retained_cone_at_least_half")
    # Do not build the enormous lcm of all distinct prefix-mass denominators.
    # Each coefficient is already checked >= volume/3; this is a rigorous
    # same-source lower certificate, not an exact value of the inflated proxy.
    common_source_coefficient_lower = sum_vol/3
    ck(common_source_coefficient_lower >= F(J-1, 12*n),
       "common_source_proxy_coefficient_lower")
    proxy_lower_from_actual_volume = common_source_coefficient_lower**2
    lower = F((J-1)**2, 144*n*n)
    ck(proxy_lower_from_actual_volume >= lower,
       "same_original_central_source_square_lower")
    global_volume_bound = F(2*J)
    ratio_lower = lower/global_volume_bound
    ck(ratio_lower >= F(J, 1152*n*n), "exponential_dimension_ratio_lower")
    # This construction does not refute a W-normalized estimate.
    W = 1+delta+6*4**(n-1)
    ck(W == 1+delta+F(6*4**(n-1)), "complete_input_normalization")
    ck(global_volume_bound == F(4*2**(n-1)),
       "full_receiver_support_bounding_box_volume")
    rounds.append({
        "n": n, "J": J, "status": "PASS", "checks": count-before,
        "delta": fs(delta), "whole_source_mass": fs(W),
        "minimum_retained_radial_width": fs(min_width),
        "largest_actual_cone_block_volume": cert(largest_volume),
        "summed_actual_cone_volume": cert(sum_vol),
        "common_source_proxy_lower_from_actual_volume":
            cert(proxy_lower_from_actual_volume),
        "proved_proxy_lower": fs(lower),
        "proved_global_E_volume_upper": fs(global_volume_bound),
        "proved_ratio_lower": fs(ratio_lower),
        "weaker_simple_ratio_lower": fs(F(J, 1152*n*n)),
        "actual_FIRST_certified": False,
        "W_counterexample_certified": False,
    })

result = {
    "status": "PASS", "new_guard_only": True, "total_checks": count,
    "rounds": rounds, "runtime_seconds": time.time()-started,
    "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    "registration_sha256": sha256(REG.read_bytes()).hexdigest(),
    "scope": registration["new_scope"],
    "not_certified": registration["not_certified"],
    "L1_status": "Analytic stable compact-subcone packet transfer only; no numerical packet simulation.",
}
OUT.write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n")
print(json.dumps({"status": result["status"], "total_checks": count,
                  "rounds": [{"n": r["n"], "J": r["J"],
                              "checks": r["checks"],
                              "ratio_lower": r["proved_ratio_lower"]}
                             for r in rounds],
                  "runtime_seconds": result["runtime_seconds"],
                  "script_sha256": result["script_sha256"]}))
