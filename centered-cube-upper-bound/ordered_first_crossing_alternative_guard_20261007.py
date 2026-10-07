#!/usr/bin/env python3
"""New first-receiver-crossing scalar receipts and original-symbol diagnostics."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import time
import numpy as np

BASE = Path(__file__).resolve().parent
REG = BASE / "ordered_first_crossing_alternative_registration_20261007.json"
OUT = BASE / "ordered_first_crossing_alternative_results_20261007.json"
reg = json.loads(REG.read_text())
assert reg["status"] == "registered_before_execution"
start = time.time()


def scalar_round(n):
    count = 0
    names = []
    for j in range(n):
        v = F(j, n)
        for u in [F(1), 1 + F(1, n), F(2), F(8)]:
            r = (1+u)/(1+v)
            for name, result in [
                ("high_ratio_ge_one", r >= 1),
                ("log_linear_coefficient_exact", (1+v)*(r-1) == u-v),
                ("hinge_cross_positive", (int(u >= 1)-int(v >= 1))*(u-v) >= 0),
                ("firsthit_test_cap", 1+v < 2),
            ]:
                assert result, (n, j, u, name)
                count += 1
                names.append(name)
    for u in [F(0), F(1,4), F(1), F(3,2), F(8)]:
        for v in [F(0), F(1,4), F(1), F(3,2), F(8)]:
            value = (int(u >= 1)-int(v >= 1))*(u-v)
            assert value >= 0
            count += 1
            names.append("hinge_pair_including_threshold_equality")
    assert F(4,5)-F(7,20)-F(2,5) == F(1,20)
    count += 1
    return {"n": n, "exact_checks": count, "status": "exact_PASS",
            "check_name_sha256": hashlib.sha256("\n".join(names).encode()).hexdigest()}


def B(v):
    v = np.asarray(v, dtype=np.float64)
    ans = np.zeros_like(v)
    mask = v > 0
    ans[mask] = v[mask]/np.log1p(v[mask])-1
    return ans


def phi(v):
    return v-np.log1p(v)


def model_round(cfg):
    n, N, panels = cfg["n"], cfg["phase_points_per_axis"], cfg["time_panels"]
    phase = 2*np.pi*np.arange(N)/N
    theta, phase_fast = np.meshgrid(phase, phase, indexing="ij")
    A = (7/20)*np.cos(theta)
    D = (2/5)*np.cos(phase_fast)
    gs, gf = np.log(17)/16, np.log(4097)/4096

    def response(z):
        e = np.exp(-z)
        slow = gs+(1-gs)*e
        fast = gf+(1-gf)*e
        value = 4/5 + A*slow-D*fast**n
        derivative = e*(-A*(1-gs)+D*n*(1-gf)*fast**(n-1))
        return value, derivative

    f, _ = response(0.0)
    terminal, _ = response(1.0)
    peak_value = np.maximum(f, terminal)
    peak_time = np.where(terminal > f, 1.0, 0.0)
    # All possible interior maxima are found analytically. If A,D<0 the
    # stationary point is a minimum; the other sign cases are monotone.
    stationary_max = (A > 0) & (D > 0)
    critical_power = np.ones_like(A)
    critical_power[stationary_max] = (
        A[stationary_max]*(1-gs)/(D[stationary_max]*n*(1-gf)))
    critical_fast = critical_power**(1/(n-1))
    critical_e = (critical_fast-gf)/(1-gf)
    valid = stationary_max & (critical_e > math.exp(-1)) & (critical_e < 1)
    safe_e = np.clip(critical_e, math.exp(-1), 1)
    critical_time = -np.log(safe_e)
    critical_value, critical_derivative = response(critical_time)
    use = valid & (critical_value > peak_value)
    peak_value = np.where(use, critical_value, peak_value)
    peak_time = np.where(use, critical_time, peak_time)

    initial_high = f >= 1
    hit = (~initial_high) & (peak_value >= 1)
    # This bracket is the increasing interval ending at the exact maximum.
    # Equality at a tangent maximum or at Z belongs to closed first-touch.
    lo, hi = np.zeros_like(A), peak_time.copy()
    for _ in range(reg["root_bisection_iterations"]):
        mid = (lo+hi)/2
        vmid, _ = response(mid)
        lo = np.where(hit & (vmid < 1), mid, lo)
        hi = np.where(hit & (vmid >= 1), mid, hi)
    tau = np.where(initial_high, 0.0, np.where(hit, hi, np.inf))
    nonhit = (~initial_high) & (~hit)
    touch_value, _ = response(np.where(hit, tau, 0))
    touch_residual = float(np.max(np.abs(touch_value[hit]-1))) if np.any(hit) else 0.
    derivative_residual = float(np.max(np.abs(critical_derivative[valid]))) if np.any(valid) else 0.
    closed_last_time = int(np.count_nonzero(hit & (tau == 1)))
    tangency_grid = int(np.count_nonzero(hit & valid & (peak_value == 1)))

    freq = np.fft.fftfreq(N)*N
    bj = B((4*freq[:,None]+64*freq[None,:])**2)
    bl = B((64*freq[None,:])**2)
    keys = ["flow", "entropy_rate", "hinge_flux", "C", "C_high",
            "C_rec_signed", "C_alive", "R_hist", "Q_refill",
            "B_rec", "weighted_bregman", "active_fraction",
            "inactive_low_fraction"]
    streams = {k: [] for k in keys}
    max_pre_hit_violation = 0.
    max_flow_generator_residual = 0.
    max_J_constant_residual = 0.
    # The nonlinear sharp-gate FFT is a diagnostic, not a positive discrete
    # Markov replacement. The multiplier is the actual original kernel.
    for step in range(panels+1):
        z = step/panels
        v, vp = response(z)
        d = math.exp(-z)
        jm = 1/(1+d*bj)-1 + (n-1)*(1/(1+d*bl)-1)

        def J(a):
            return np.fft.ifft2(np.fft.fft2(a)*jm).real

        H = ((~initial_high) & (z < tau)).astype(float)
        test = H*(1+v)
        logv = np.log1p(v)
        high = (v >= 1).astype(float)
        rec = (1-H)*(v < 1)
        Jlog = J(logv)
        C = float(np.mean(test*Jlog))
        C_high = float(np.mean(test*(J(logv*high)-logv*J(high))))
        C_rec = float(np.mean(test*(J(logv*rec)-logv*J(rec))))
        R_hist = float(np.mean(H*(J(v*rec)-v*J(rec))))
        Q_refill = float(np.mean(rec*(J(v*high)-v*J(high))))
        values = {
            "flow": float(np.mean(H*vp)),
            "entropy_rate": float(-np.mean(v/(1+v)*vp)),
            "hinge_flux": float(-np.mean(high*vp)),
            "C": C, "C_high": C_high, "C_rec_signed": C_rec,
            "C_alive": C-C_high-C_rec, "R_hist": R_hist,
            "Q_refill": Q_refill, "B_rec": R_hist-C_rec,
            "weighted_bregman": float(np.mean(test*(vp/(1+v)-Jlog))),
            "active_fraction": float(np.mean(H)),
            "inactive_low_fraction": float(np.mean(rec)),
        }
        for k in keys:
            streams[k].append(values[k])
        if np.any(H):
            max_pre_hit_violation = max(max_pre_hit_violation, float(np.max(v[H > 0]-1)))
        max_flow_generator_residual = max(max_flow_generator_residual,
                                         float(np.max(np.abs(J(v)-vp))))
        max_J_constant_residual = max(max_J_constant_residual,
                                     float(np.max(np.abs(J(np.ones_like(v))))))

    weights = np.ones(panels+1)
    weights[1:-1:2] = 4
    weights[2:-1:2] = 2
    integral = {k: float(np.dot(weights, x)/(3*panels))
                for k, x in streams.items()}
    initial_low_mass = float(np.mean(f*(~initial_high)))
    nonhit_terminal_mass = float(np.mean(terminal*nonhit))
    stopped_wealth = float(np.mean(hit)+nonhit_terminal_mass-initial_low_mass)
    entropy_drop = float(np.mean(phi(f)-phi(terminal)))
    hinge_drop = float(np.mean(np.maximum(f-1,0)-np.maximum(terminal-1,0)))
    deficit_terminal = float(np.mean((tau <= 1)*np.maximum(1-terminal,0)))
    Evolume = float(np.mean(peak_value > 1))
    W = float(np.mean(f))
    residuals = {
        "flow_minus_exact_stopped_wealth": integral["flow"]-stopped_wealth,
        "entropy_integral_minus_exact_drop": integral["entropy_rate"]-entropy_drop,
        "hinge_integral_minus_exact_drop": integral["hinge_flux"]-hinge_drop,
        "flow_minus_bregman_and_commutator": integral["flow"]-integral["weighted_bregman"]-integral["C"],
        "hist_minus_refill_minus_terminal_deficit": integral["R_hist"]-integral["Q_refill"]-deficit_terminal,
        "original_J_v_minus_exact_derivative_max": max_flow_generator_residual,
        "J_constant_max": max_J_constant_residual,
        "first_touch_root_residual_max": touch_residual,
        "stationary_derivative_residual_max": derivative_residual,
        "prehit_response_minus_eta_max": max_pre_hit_violation,
    }
    tol = reg["numeric_screen_absolute_tolerance"]
    screens = {
        "all_finite": all(math.isfinite(x) for x in integral.values()),
        "source_positive": float(f.min()) >= 1/20-1e-12,
        "source_mass_4_over_5": abs(W-4/5) < 1e-12,
        "flow_identity": abs(residuals["flow_minus_exact_stopped_wealth"]) <= tol,
        "entropy_identity": abs(residuals["entropy_integral_minus_exact_drop"]) <= tol,
        "hinge_identity": abs(residuals["hinge_integral_minus_exact_drop"]) <= tol,
        "bregman_identity": abs(residuals["flow_minus_bregman_and_commutator"]) <= tol,
        "history_deficit_identity": abs(residuals["hist_minus_refill_minus_terminal_deficit"]) <= tol,
        "generator_exact_modes": max_flow_generator_residual < 1e-9,
        "first_touch_peak_bracket": touch_residual < 1e-9,
        "stationary_point": derivative_residual < 1e-9,
        "prehit_below_eta": max_pre_hit_violation < 1e-9,
        "hinge_fee_le_W": -tol <= integral["hinge_flux"] <= W+tol,
        "high_branch_le_hinge": -tol <= integral["C_high"] <= integral["hinge_flux"]+tol,
        "alive_internal_nonpositive": integral["C_alive"] <= tol,
        "bregman_positive": integral["weighted_bregman"] >= -tol,
        "rec_bregman_positive": integral["B_rec"] >= -tol,
        "terminal_weak_interface": Evolume <= 4*W+integral["C_rec_signed"]+tol,
    }
    return {
        **cfg, "scope": "original-symbol 2-phase torus diagnostic, not interval",
        "W": W, "Evolume": Evolume,
        "initial_high_fraction": float(np.mean(initial_high)),
        "new_hit_fraction": float(np.mean(hit)),
        "nonhit_fraction": float(np.mean(nonhit)),
        "closed_terminal_hit_grid_count": closed_last_time,
        "tangent_hit_grid_count": tangency_grid,
        "first_touch_by_exact_critical_point_then_bisection": True,
        "first_touch_time_grid_used_for_peak_search": False,
        "terminal_deficit": deficit_terminal,
        "exact_phase_grid_entropy_drop": entropy_drop,
        "exact_phase_grid_hinge_drop": hinge_drop,
        "integrals": integral, "residuals": residuals,
        "screens": screens,
        "screens_all_pass": all(screens.values()),
        "sample_times": {str(j/panels): {k: streams[k][j] for k in keys}
                         for j in [0, panels//4, panels//2, 3*panels//4, panels]},
    }


exact = [scalar_round(n) for n in reg["exact_scalar_rounds"]]
numeric = [model_round(cfg) for cfg in reg["models"]]
payload = {
    "status": "terminal_recorded",
    "exact_status": "all_exact_PASS",
    "exact_checks": sum(r["exact_checks"] for r in exact),
    "numeric_screens_all_pass": all(r["screens_all_pass"] for r in numeric),
    "elapsed_seconds": time.time()-start,
    "random_seed": None,
    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "registration_sha256": hashlib.sha256(REG.read_bytes()).hexdigest(),
    "old_data_rerun": False,
    "exact": exact, "numeric": numeric,
    "limitations": ["sharp gate FFT and Simpson not interval-certified",
                    "periodic 2-phase input is not Rn L1",
                    "no actual FIRST/CP/GP/history claim",
                    "general proof is analytic, no new weak endpoint"],
}
OUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2)+"\n")
print(json.dumps({"status": payload["status"],
                  "exact_checks": payload["exact_checks"],
                  "numeric_screens_all_pass": payload["numeric_screens_all_pass"],
                  "elapsed_seconds": payload["elapsed_seconds"],
                  "output": str(OUT)}, ensure_ascii=False))
