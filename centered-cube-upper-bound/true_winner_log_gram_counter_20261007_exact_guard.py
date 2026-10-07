#!/usr/bin/env python3
"""New exact certificates for the continuous true-winner log-Gram counter.

These are original hard-kernel component checks with an exact L1 packet lift,
not actual FIRST/history samples. No floating logarithm is used as a proof.
"""
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
import json
import time

BASE = Path(__file__).resolve().parent
PREFIX = "true_winner_log_gram_counter_20261007"
REG = BASE / (PREFIX + "_registration.json")
OUT = BASE / (PREFIX + "_results.json")
if OUT.exists():
    raise SystemExit("Refusing to overwrite a terminal result.")
registration = json.loads(REG.read_text())
assert registration["status"] == "registered_before_execution"
assert registration["rounds"] == [16, 64, 256]

total_checks = 0


def check(condition, name):
    global total_checks
    total_checks += 1
    if not condition:
        raise AssertionError(name)


def small(x):
    return str(x.numerator) + "/" + str(x.denominator)


def receipt(x):
    raw = small(x)
    return {
        "exact_fraction_sha256": sha256(raw.encode()).hexdigest(),
        "decimal_preview_not_certificate": float(x),
        "numerator_bits": x.numerator.bit_length(),
        "denominator_bits": x.denominator.bit_length(),
    }


started = time.time()
check(registration["arithmetic"] == "fractions.Fraction, all deterministic",
      "registration_arithmetic")
rounds = []
for J in registration["rounds"]:
    before = total_checks
    h = F(1, 100 * J)
    eps = F(1, 10**6 * J**2)
    ss = [i * h for i in range(J + 1)]
    zz = [(F(9, 10) * s - 15 * s*s) / (1+s) for s in ss]
    mm = [1 / (J * (1+s)) for s in ss]
    aa = [mm[i] - mm[i+1] for i in range(J)] + [mm[J]]
    suffix_masses = [F(0)]*(J+1)
    running_mass = F(0)
    for i in range(J, -1, -1):
        running_mass += aa[i]
        suffix_masses[i] = running_mass
    check(sum(aa, F(0)) == F(1, J), "source_mass_telescope")
    check(zz[0] == 0 and zz[J] == F(3, 404), "position_endpoints")
    check(F(2) + zz[J] < F(201, 100), "complete_atomic_E_volume")
    check(F(6) + sum(aa, F(0)) == F(6) + F(1, J), "full_background_mass")
    check(F(1197, 2000) > F(29, 50) * F(101, 100)**2,
          "analytic_derivative_uniform_residual")
    check(-2-eps > -3 and 2+zz[J]+eps < 3,
          "all_queries_on_complete_E_have_full_background")
    check(F(2)+zz[J]+2*eps < F(201, 100),
          "complete_L1_E_volume")
    check(F(3, 5)-eps-zz[J]-eps > F(1, 2),
          "all_L1_partial_bands_start_above_one")
    check(F(9, 10)-eps-zz[0]+eps < 1,
          "all_L1_partial_bands_end_below_two")

    for i in range(J + 1):
        s, z, m, alpha = ss[i], zz[i], mm[i], aa[i]
        derivative = (F(9, 10)-30*s-15*s*s)/(1+s)**2
        check(derivative >= F(29, 50), "grid_derivative_lower")
        check(alpha > 0, "positive_mass")
        check(suffix_masses[i] == m, "each_suffix_mass")
        check(alpha >= F(1, 103*J**2), "packet_mass_minimum")
        check(alpha/(4*eps) > F(1, J),
              "partial_band_strict_signal_derivative")
        check(F(1) < 2*(F(3, 5)-z) < F(2),
              "atomic_arrival_left_region")
        check(F(1) < 2*(F(9, 10)-z) < F(2),
              "atomic_arrival_right_region")
        if i < J:
            check(zz[i+1]-zz[i] >= F(29, 50)*h,
                  "integrated_derivative_gap")
            check(zz[i+1]-zz[i] > 2*eps, "disjoint_packets")
            check(aa[i] == h/(J*(1+ss[i])*(1+ss[i+1])),
                  "exact_packet_mass_identity")

    def inv(k, x):
        # Reciprocal of the excess over the true background level.
        return 2*(x-zz[k])/mm[k]

    intervals = []
    for i in range(1, J):
        left = F(9, 10)-30*(ss[i]+h/2)
        right = F(9, 10)-30*(ss[i]-h/2)
        mid = F(9, 10)-30*ss[i]
        v = right-left
        check(v == F(3, 10*J), "unique_winner_interval_volume")
        check(F(3, 5) < left < mid < right < F(9, 10),
              "internal_interval_location")
        check(v/mm[i] >= F(3, 10), "log_argument_lower")
        check(v/(v+mm[i]) >= F(3, 13),
              "analytic_log_lower_via_t_over_one_plus_t")
        check(v/(1+mm[i]) < v, "exact_coefficient_keeps_scale_floor")
        check(inv(i, left) == inv(i+1, left), "left_adjacent_tie")
        check(inv(i, right) == inv(i-1, right), "right_adjacent_tie")
        x_l1 = mid-eps
        selected_l1 = 2*(x_l1-zz[i]+eps)
        check(selected_l1 == 2*(mid-zz[i]),
              "L1_complete_capture_scale_same_envelope")
        check(F(1) < selected_l1 < F(2), "L1_selected_scale_domain")
        for k in range(J + 1):
            quadratic = 2*J*(mid+(mid-F(9, 10))*ss[k]+15*ss[k]**2)
            check(inv(k, mid) == quadratic, "inverse_quadratic_identity")
            check(k == i or inv(k, mid) > inv(i, mid),
                  "strict_internal_true_winner_all_candidates")
            check(inv(k, left) >= inv(i, left), "left_tie_global_envelope")
            check(inv(k, right) >= inv(i, right), "right_tie_global_envelope")
            if k not in (i, i+1):
                check(inv(k, left) > inv(i, left), "left_only_two_ties")
            if k not in (i, i-1):
                check(inv(k, right) > inv(i, right), "right_only_two_ties")
            check(2*(x_l1-zz[k]+eps)/mm[k] == inv(k, mid),
                  "all_L1_complete_capture_candidates_same_envelope")
            # At selected end all packets i..J are complete, earlier absent.
            edge = x_l1-selected_l1/2
            if k >= i:
                check(edge <= zz[k]-eps, "selected_suffix_packet_complete")
            else:
                check(edge >= zz[k]+eps, "earlier_packet_uncaptured")
        intervals.append((left-eps, right-eps))

    lower_log = F(3, 13)
    proxy_lower = F(0)
    for i in range(1, J):
        for j in range(1, J):
            common = suffix_masses[max(i, j)]
            check(common == mm[max(i, j)], "actual_suffix_overlap_mass")
            check(common >= aa[J], "last_packet_common_to_all_blocks")
            proxy_lower += common * lower_log**2
    common_source_lower = aa[J] * ((J-1)*lower_log)**2
    check(common_source_lower == F(900*(J-1)**2, 17069*J),
          "explicit_linear_counter_lower")
    check(proxy_lower >= common_source_lower, "full_overlap_contains_common_source")
    check(common_source_lower >= F(900, 4*17069)*J,
          "uniform_J_linear_lower_for_J_ge_two")
    check(sum((r-l for l, r in intervals), F(0)) == F(3*(J-1), 10*J),
          "disjoint_true_receiver_block_total")
    check(sum((r-l for l, r in intervals), F(0)) < F(2)+zz[J]+2*eps,
          "blocks_inside_complete_E")
    rounds.append({
        "J": J,
        "checks": total_checks-before,
        "status": "PASS",
        "epsilon": small(eps),
        "complete_source_mass": small(F(6)+F(1, J)),
        "atomic_E_volume": small(F(2)+zz[J]),
        "L1_E_volume": small(F(2)+zz[J]+2*eps),
        "true_internal_block_volume_each": small(F(3, 10*J)),
        "common_source_proxy_lower": small(common_source_lower),
        "common_source_proxy_lower_preview": float(common_source_lower),
        "all_suffix_rational_proxy_lower": receipt(proxy_lower),
        "minimum_packet_derivative_margin": small(min(aa)/(4*eps)-F(1, J)),
        "finite_scale_family_certified": False,
        "actual_FIRST_certified": False,
    })

result = {
    "status": "PASS",
    "scope": registration["scope"],
    "new_guard_only": True,
    "random_seed": None,
    "total_checks": total_checks,
    "rounds": rounds,
    "runtime_seconds": time.time()-started,
    "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    "registration_sha256": sha256(REG.read_bytes()).hexdigest(),
    "true_L1_transfer": "Exact partial-band monotonicity and complete-packet endpoint, not an atomic-limit-only claim.",
    "logarithm_certificate": "log(13/10)>=3/13 is proved analytically; rational checks do not evaluate logarithms.",
    "limits": registration["not_certified"],
}
OUT.write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n")
print(json.dumps({"status": result["status"], "total_checks": total_checks,
                  "rounds": [{"J": x["J"], "checks": x["checks"],
                              "proxy_lower": x["common_source_proxy_lower"]}
                             for x in rounds],
                  "runtime_seconds": result["runtime_seconds"],
                  "script_sha256": result["script_sha256"]}, ensure_ascii=False))
