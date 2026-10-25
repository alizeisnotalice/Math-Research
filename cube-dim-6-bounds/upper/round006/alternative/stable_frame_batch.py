"""Finite alpha=3/2 body32 diagnostics. No certified integration or C6 bound.

One model, 16 spatial parameters, two 4096-point scrambled Sobol replicates.
Time integration is analytic using independent positive-stable Gaussian mixing.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import time

import numpy as np
import scipy
from scipy.special import gamma
from scipy.stats import levy_stable, qmc
import sympy as sp

HERE = Path(__file__).resolve().parent
START = time.monotonic()
ALPHA = 1.5
BETA = ALPHA / 2
V = np.array([(1,) + s for s in itertools.product((-1, 1), repeat=5)], dtype=float)
BASE_IDS = [0, 3, 5, 14, 22, 25]
FACE = np.column_stack((np.ones(8), qmc.Sobol(5, scramble=True, seed=6006003).random_base2(3)))
ANCHORS = np.array([
    [1, 0, 0, 0, 0, 0],
    [1, 1, 0, 0, 0, 0],
    [1, 1, 1, 0, 0, 0],
    [1, .8, .6, .4, .2, 0],
    [1, .99, .9, .8, .45, .1],
    [1, .9, .9, .9, .9, .9],
    [1, .99, .99, .99, .99, .9],
], dtype=float)
POINTS = np.vstack((FACE, ANCHORS))


def json_hash(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def positive_stable(uniforms):
    """Kanter formula, Laplace transform exp(-lambda**beta).

    All inputs are open-interval uniforms; no clipping is silently performed.
    """
    assert np.all((uniforms > 0) & (uniforms < 1))
    u = np.pi * uniforms[:, :32]
    e = -np.log(uniforms[:, 32:])
    return (np.sin(BETA * u) / np.sin(u) ** (1 / BETA)
            * (np.sin((1 - BETA) * u) / e) ** ((1 - BETA) / BETA))


def batch(seed):
    uniforms = qmc.Sobol(64, scramble=True, seed=seed).random_base2(12)
    stable = positive_stable(uniforms)
    covariance = np.einsum("nv,vi,vj->nij", stable, V, V)
    eigenvalues = np.linalg.eigvalsh(covariance)
    assert np.all(eigenvalues > 0)
    sign, logdet = np.linalg.slogdet(covariance)
    assert np.all(sign == 1)
    inverse = np.linalg.inv(covariance)
    quadratic = np.einsum("pi,nij,pj->np", POINTS, inverse, POINTS)
    constant = ALPHA * gamma((6 - ALPHA) / 2) / (2 ** (ALPHA + 1) * np.pi ** 3)
    values = constant * np.exp(-logdet[:, None] / 2) * quadratic ** (-(6 - ALPHA) / 2)
    assert np.all(np.isfinite(values)) and np.all(values > 0)
    means = values.mean(axis=0)
    sample_g = means.min()
    face_mean_power = np.mean(means[:8] ** (4 / 3))
    cost = 4 ** (4 / 3) * face_mean_power / sample_g ** (4 / 3)
    # These checks test normalization of the sampling implementation only.
    laplace = {str(lam): {"sample":float(np.mean(np.exp(-lam * stable))),
                         "exact":float(np.exp(-lam ** BETA))}
               for lam in [.25, 1., 4.]}
    negative_moment = {str(r): {"sample":float(np.mean(stable ** (-r))),
                               "exact":float(gamma(1 + r / BETA) / gamma(1 + r))}
                       for r in [1., 2.]}
    return {
        "seed":seed, "N":len(stable),
        "green_estimates":means.tolist(),
        "sample_standard_error_not_qmc_error_bound":(values.std(axis=0, ddof=1) / math.sqrt(len(values))).tolist(),
        "maximum_single_sample_fraction":(values.max(axis=0) / values.sum(axis=0)).tolist(),
        "empirical_minimum":float(sample_g), "empirical_minimum_index":int(means.argmin()),
        "eight_face_sample_mean_power":float(face_mean_power),
        "formal_cost_no_certificate":float(cost),
        "laplace_checks":laplace, "negative_moment_checks":negative_moment,
        "minimum_eigenvalue":float(eigenvalues.min()),
        "maximum_spectral_condition_number":float((eigenvalues[:, -1] / eigenvalues[:, 0]).max()),
        "uniform_input_sha256":hashlib.sha256(uniforms.tobytes()).hexdigest(),
    }


def exact_checks():
    alpha = sp.Rational(3, 2)
    d = sp.Integer(6)
    q = d / (d - alpha)
    c, g = sp.symbols("c g", positive=True)
    cutoff_cost = (d - alpha) / alpha * c ** (-alpha / (d - alpha)) / (g - c)
    optimum = alpha * g / d
    critical_derivative = sp.simplify(sp.diff(cutoff_cost, c).subs(c, optimum))
    optimized_factor = sp.simplify(cutoff_cost.subs(c, optimum) * g ** q)
    assert critical_derivative == 0
    assert sp.simplify(optimized_factor - (d / alpha) ** q) == 0
    determinant = sp.Matrix(V[BASE_IDS].astype(int)).det()
    assert abs(determinant) == 160
    assert np.array_equal(V.T @ V, 32 * np.eye(6))
    assert all(sp.Matrix(V[[i, j]].astype(int)).rank() == 2
               for i, j in itertools.combinations(range(32), 2))
    density_sup = gamma(1 / ALPHA) / (ALPHA * np.pi)
    # Proven in report.md by contour rotation; not fitted to stable.pdf.
    density_tail_envelope = gamma(1 + ALPHA) / (np.pi * np.sin(np.pi/(2*ALPHA)) ** (1+ALPHA))
    one_big = density_tail_envelope * density_sup ** 5 * 64 ** (1+ALPHA)
    two_big = (31 * density_tail_envelope ** 2 * density_sup ** 4
               * 64 ** (1+ALPHA) * (62 * np.sqrt(6)) ** (1+ALPHA))
    green_line_envelope = 3 * two_big + .75 * one_big
    return {"alpha":str(alpha), "green_degree":str(alpha-d), "q":str(q),
            "cutoff_over_g":str(alpha/d), "optimized_cost_factor":str(optimized_factor),
            "base_indices":BASE_IDS, "signed_determinant":int(determinant),
            "pair_count_rank_two":496,
            "short_time_exponents":{str(k):str(k-(d-k)/alpha) for k in range(1,7)},
            "critical_alpha_rank_one":str((d-1)/2),
            "critical_alpha_rank_two":str((d-2)/3),
            "explicit_analytic_constants_not_interval_evaluated":{
                "M_density_sup":float(density_sup),
                "A_density_tail_envelope":float(density_tail_envelope),
                "C1_one_big":float(one_big), "C2_two_big":float(two_big),
                "D_G_le_D_rho_minus2":float(green_line_envelope),
                "positive_corner_cap_I_le_constant_delta_7over3":
                    float((6*green_line_envelope)**(4/3)*np.pi**2/28)}}


def tail_and_stratum_checks():
    # Density of the six-vector invertible subframe at points with k nonzero
    # frame coordinates: exact positive product, but not the full 32-frame law.
    rows = []
    tail_constant = gamma(1 + ALPHA) * np.sin(np.pi * ALPHA / 2) / np.pi
    for k in [1,2]:
        times = np.array([1e-2, 1e-3, 1e-4, 1e-5])
        pdfs = []
        for t in times:
            scale = t ** (1 / ALPHA)
            at_one = levy_stable.pdf(1 / scale, ALPHA, 0) / scale
            at_zero = gamma(1 / ALPHA) / (ALPHA * np.pi * scale)
            pdfs.append(float(at_one ** k * at_zero ** (6-k) / 160))
        slopes = np.diff(np.log(pdfs)) / np.diff(np.log(times))
        rows.append({"rank":k,"times":times.tolist(), "subframe_heat_density":pdfs,
                     "log_slopes":slopes.tolist(), "predicted_exponent":k-(6-k)/ALPHA})
    ratios = [{"x":x, "density_tail_ratio":float(levy_stable.pdf(x,ALPHA,0)*x**(1+ALPHA)/tail_constant)}
              for x in [10.,30.,100.,300.]]
    return {"one_dimensional_tail_constant":float(tail_constant), "tail_ratios":ratios,
            "subframe_strata":rows}


def main():
    inputs = {"model":"body32 equal-weight symmetric stable", "alpha":"3/2", "V":V.astype(int).tolist(),
              "spatial_parameters":POINTS.tolist()+[[1,1,1,1,1,1]], "replicate_seeds":[6006151,6006152],
              "N_each":4096, "point_15_status":"exact line singularity G=+infinity, not sampled"}
    result = {
        "status":"finite numerical feasibility only; no continuous g/I certificate; no C6 new bound",
        "inputs":inputs, "input_sha256":json_hash(inputs),
        "registration_sha256":hashlib.sha256((HERE/"registration.json").read_bytes()).hexdigest(),
        "source_code_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "versions":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"sympy":sp.__version__},
        "exact_checks":exact_checks(), "tail_and_stratum_checks":tail_and_stratum_checks(),
    }
    result["replicates"] = [batch(seed) for seed in inputs["replicate_seeds"]]
    a,b = [np.array(row["green_estimates"]) for row in result["replicates"]]
    result["replicate_relative_difference"] = (2*abs(a-b)/(a+b)).tolist()
    result["max_replicate_relative_difference"] = float(max(result["replicate_relative_difference"]))
    result["wall_seconds"] = time.monotonic()-START
    assert result["wall_seconds"] < 1200
    (HERE/"stable_frame_results.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"input_sha256":result["input_sha256"],
                      "costs":[r["formal_cost_no_certificate"] for r in result["replicates"]],
                      "empirical_minimum":[r["empirical_minimum"] for r in result["replicates"]],
                      "max_replicate_relative_difference":result["max_replicate_relative_difference"],
                      "wall_seconds":result["wall_seconds"]},indent=2))


if __name__ == "__main__":
    main()
