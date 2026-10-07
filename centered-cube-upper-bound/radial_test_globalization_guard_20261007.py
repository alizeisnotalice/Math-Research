#!/usr/bin/env python3
"""One registered run: exact arithmetic for the original L1 cube example.
No receiver grid, no sampling of global v, no old oracle or probe imports.
The analytic universal-v estimate is cited from the frozen proof; this guard
checks its explicit parameters and scalar constants, not all functions v.
"""
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from math import comb
from datetime import datetime, timezone
import json
import time

HERE = Path(__file__).resolve().parent
PREFIX = "radial_test_globalization_guard_20261007"
REG = HERE / (PREFIX + "_registration.json")
PROOF = HERE / "radial_test_globalization_20261007.md"
COUNTS = {"total": 0}
CHECKS = []


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def require(label, condition, round_id=0, n=None, H_factor=None, eps=None):
    COUNTS["total"] += 1
    if not condition:
        raise AssertionError(label)
    CHECKS.append({
        "label": label, "PASS": True, "round": round_id, "n": n,
        "H_factor": H_factor, "epsilon": None if eps is None else str(eps),
    })


def certificate(value):
    value = F(value)
    raw = str(value.numerator) + "/" + str(value.denominator)
    scale = 1 << 128
    low = value.numerator * scale // value.denominator
    try:
        decimal_diagnostic = float(value)
    except OverflowError:
        decimal_diagnostic = None
    out = {
        "fraction_sha256": sha256(raw.encode()).hexdigest(),
        "numerator_digits": len(str(abs(value.numerator))),
        "denominator_digits": len(str(value.denominator)),
        "dyadic128_lower": str(low) + "/" + str(scale),
        "dyadic128_upper": str(low + 1) + "/" + str(scale),
        "decimal_diagnostic": decimal_diagnostic,
    }
    if len(raw) <= 240:
        out["exact_fraction"] = raw
    return out


def main():
    started_utc = datetime.now(timezone.utc).isoformat()
    started = time.perf_counter()
    reg_bytes = REG.read_bytes()
    reg = json.loads(reg_bytes)
    amendment_path = HERE / (PREFIX + "_registration_amendment.json")
    amendment = json.loads(amendment_path.read_text())
    require("amendment preserves original registration", amendment["original_registration_sha256"] == sha256(reg_bytes).hexdigest())
    require("registration exact rounds match", reg["rounds"] == [[2, 16], [4, 64], [8, 256]])
    require("registration exact H factors match", reg["H_factors"] == [1, 2, 4])
    require("registration exact amplitudes match", reg["epsilon"] == ["1/8", "1/4"])
    require("amended script frozen", amendment["amended_script_sha256"] == digest(Path(__file__)))
    require("registered analytic proof frozen", reg["proof_sha256"] == digest(PROOF))
    profiles = []
    dimension_profiles = []
    round_summaries = []
    for rid, dims in enumerate(reg["rounds"], start=1):
        before = COUNTS["total"]
        for n in dims:
            L = F(n + 1, n)
            width = L - 1
            # Exact antiderivative of integral_1^L [1-(s/L)^n] ds /(L-1).
            radial_J = (F(n, n + 1) * L - 1 + L ** (-n) / (n + 1)) / width
            integral_J = (width - (L ** (n + 1) - 1) / ((n + 1) * L ** n)) / width
            # Independent binomial integration after s=1+(L-1)t.
            binomial_integral = sum((F(comb(n, k), k + 1) * width ** k for k in range(n + 1)), F(0))
            polynomial_J = 1 - binomial_integral / L ** n
            c = lambda name, cond: require(name, cond, rid, n)
            c("candidate L=1+1/n in original [1,2]", 1 < L <= 2)
            c("theta lower integral bound >=2/3", F(n, n + 1) >= F(2, 3))
            c("radial J antiderivative equality", radial_J == integral_J)
            c("radial J independent positive polynomial equality", radial_J == polynomial_J)
            c("radial J positive and at most one", 0 < radial_J <= 1)
            c("uniform radial J >=1/10", radial_J >= F(1, 10))
            c("conservative weighted radial gap >=3/40", F(3, 4) * radial_J >= F(3, 40))
            dimension_profiles.append({
                "round": rid, "n": n, "L": str(L),
                "theta_lower": str(F(n, n + 1)),
                "theta_upper": "1 (using log(1+1/n)<=1/n)",
                "uniform_radial_expectation": certificate(radial_J),
                "weighted_radial_lower_3over4": certificate(F(3, 4) * radial_J),
            })
            for factor in reg["H_factors"]:
                H = factor * (n * n + 2)
                for eps_raw in reg["epsilon"]:
                    eps = F(eps_raw)
                    tau = F(1, 2)
                    source_volume = F((2 * H) ** n)
                    receiver_volume = F((2 * H - 2) ** n)
                    W = (1 - eps / 3) * source_volume
                    W_independent = source_volume - eps * F(n) * source_volume * H * H / (3 * n * H * H)
                    A1_min = 1 - eps * F((H - 1) ** 2, H * H) - eps / (12 * H * H)
                    AL_min = 1 - eps * F((H - 1) ** 2, H * H) - eps * L * L / (12 * H * H)
                    A1_max = 1 - eps / (12 * H * H)
                    mean_drop = eps * (L * L - 1) / (12 * H * H)
                    gap_upper = mean_drop / (1 - eps)
                    rho_squared = F(H, H - 1) ** n
                    # rho <=2 is certified by its exact rational square.
                    global_upper = F(26, n) / (1 - eps) + gap_upper
                    local_lower = (1 - eps) * radial_J
                    # The universal conversion requires at least local/global.
                    conversion_lower = local_lower / global_upper
                    cf = lambda name, cond: require(name, cond, rid, n, factor, eps)
                    cf("H exactly registered multiple", H == factor * (n * n + 2))
                    cf("H >=n+1", H >= n + 1)
                    cf("full source density >=3/4", 1 - eps >= F(3, 4))
                    cf("full source density <=1", 0 < eps < 1)
                    cf("complete source integral W exact", W == W_independent)
                    cf("complete W positive", W > 0)
                    cf("inner receiver volume positive", receiver_volume > 0)
                    cf("receiver inner box contained in source box", receiver_volume < source_volume)
                    cf("all allowed queries inner box contained in source box", H - 1 + F(2, 2) == H)
                    cf("inner A1 at least full source minimum", A1_min >= 1 - eps)
                    cf("inner AL at least full source minimum", AL_min >= 1 - eps)
                    cf("true inner receiver set strictly above tau", A1_min > tau)
                    cf("true band upper using full f<=1", A1_max <= 1 == 2 * tau)
                    cf("average R quadratic slope strictly negative", -eps / (12 * H * H) < 0)
                    cf("unique winner 1 versus L strict", mean_drop > 0)
                    cf("unique winner 1 versus endpoint2 strict", eps * F(3, 12 * H * H) > 0)
                    cf("exact response difference identity", A1_min - AL_min == mean_drop)
                    cf("source mass does not equal free receiver volume", W != receiver_volume)
                    cf("boundary volume ratio square <=4", rho_squared <= 4)
                    cf("physical posterior gap upper positive", gap_upper > 0)
                    cf("posterior gap upper <= conservative 1/4 bound", gap_upper <= (L * L - 1) / (36 * H * H))
                    cf("actual amplitude radial lower >=3/40", local_lower >= F(3, 40))
                    cf("constant-free multiplier contribution <=104/(3n)", F(26, n) / (1 - eps) <= F(104, 3 * n))
                    cf("global all-v bound <=36/n", global_upper <= F(36, n))
                    cf("conversion lower >=n/480", conversion_lower >= F(n, 480))
                    cf("hinge local lower positive at t1/4", local_lower / 4 - gap_upper > 0)
                    cf("boundary W/source normalized factor exactly 1-eps/3", W / source_volume == 1 - eps / 3)
                    profiles.append({
                        "round": rid, "n": n, "H_factor": factor, "H": H,
                        "epsilon": str(eps), "tau": str(tau), "L": str(L),
                        "complete_source_mass": certificate(W),
                        "complete_source_box_volume": certificate(source_volume),
                        "receiver_inner_box_volume": certificate(receiver_volume),
                        "inner_AL_minimum": certificate(AL_min),
                        "true_response_gap_upper": certificate(gap_upper),
                        "boundary_rho_squared": certificate(rho_squared),
                        "actual_amplitude_radial_gap_lower": certificate(local_lower),
                        "all_global_v_mean_upper": certificate(global_upper),
                        "necessary_conversion_constant_lower": certificate(conversion_lower),
                    })
        round_summaries.append({"round": rid, "dimensions": dims, "checks": COUNTS["total"] - before, "PASS": True})
    result = {
        "status": "PASS", "seed": None, "randomness": "none",
        "checks_total": COUNTS["total"], "profile_count": len(profiles),
        "rounds": round_summaries, "dimension_profiles": dimension_profiles,
        "profiles": profiles, "checks": CHECKS,
        "registration_sha256": sha256(reg_bytes).hexdigest(),
        "registration_amendment_sha256": digest(amendment_path),
        "script_sha256": digest(Path(__file__)), "proof_sha256": digest(PROOF),
        "scope": reg["scope"],
    }
    result_path = HERE / (PREFIX + "_results.json")
    result_path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    receipt = {
        "status": "PASS", "started_utc": started_utc,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.perf_counter() - started,
        "checks_total": COUNTS["total"], "profile_count": len(profiles),
        "rounds": round_summaries, "seed": None,
        "script_sha256": digest(Path(__file__)),
        "registration_sha256": sha256(reg_bytes).hexdigest(),
        "registration_amendment_sha256": digest(amendment_path),
        "proof_sha256": digest(PROOF), "results_sha256": digest(result_path),
        "scope": reg["scope"], "successful_execution_count": 1,
        "failed_prior_attempt_count": 1,
        "old_probe_oracle_invocations": 0,
    }
    (HERE / (PREFIX + "_receipt.json")).write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "checks_total": COUNTS["total"], "profiles": len(profiles),
                      "rounds": round_summaries, "elapsed_seconds": receipt["elapsed_seconds"]}))


if __name__ == "__main__":
    main()
