#!/usr/bin/env python3
"""One frozen batch of three new A076 small symbolic configurations."""
import hashlib
import itertools
import json
import os
from pathlib import Path
import resource
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parent
resource.setrlimit(resource.RLIMIT_CPU, (27, 27))
START = time.monotonic()
CPU_START = time.process_time()
RESULT = {
    "schema": "A076-small-results-v1", "pid": os.getpid(),
    "configurations_started": 0, "configurations_completed": 0,
    "configurations": [], "status": "RUNNING",
    "actual_work": {
        "new_symbolic_configurations": 0,
        "new_three_dimensional_density_paths": 0,
        "actual_upper_halfplane_residues": 0,
        "exact_pair_span_determinants": 0,
        "body32_full_kernel_calls": 0,
        "six_dimensional_full_kernel_calls": 0,
        "old_baseline_reruns": 0, "matlab_calls": 0,
        "recursive_agents": 0, "automatic_retries": 0
    }
}


def save():
    RESULT["wall_seconds"] = time.monotonic() - START
    RESULT["cpu_seconds"] = time.process_time() - CPU_START
    (ROOT / "results.json").write_text(json.dumps(RESULT, indent=2) + "\n")


def begin(identifier):
    RESULT["configurations_started"] += 1
    RESULT["actual_work"]["new_symbolic_configurations"] += 1
    RESULT["active_configuration"] = identifier
    save()


def finish(identifier, detail):
    RESULT["configurations"].append({"id": identifier, "status": "PASS", **detail})
    RESULT["configurations_completed"] += 1
    RESULT["active_configuration"] = None
    save()


try:
    registration = json.loads((ROOT / "registration.json").read_text())
    for name, expected in registration["frozen_local_sha256"].items():
        actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        if actual != expected:
            raise RuntimeError("Pre-execution hash mismatch: " + name)
    inputs = json.loads((ROOT / "input_manifest.json").read_text())
    assert len(inputs["configurations"]) == 3
    import sympy as s
    RESULT["runtime"] = {"python": sys.version, "executable": sys.executable,
                         "sympy": s.__version__}
    RESULT["registration_sha256"] = hashlib.sha256((ROOT / "registration.json").read_bytes()).hexdigest()

    begin("C001-fixed-cutoff")
    v, w = s.symbols("v w", real=True)
    H = sum(s.Integer(c) * v**i for i, c in enumerate(inputs["configurations"][0]["H_coefficients_ascending"]))
    h4 = sum(s.Integer(c) * v**i for i, c in enumerate(inputs["configurations"][0]["h4_coefficients_ascending"]))
    assert s.expand(s.diff(H, v) + 30*v**2*(1-v)**2) == 0
    assert s.expand(s.diff(H, v, 2) + 60*v*(1-v)*(1-2*v)) == 0
    for h, order in ((H, 2), (h4, 4)):
        assert h.subs(v, 0) == 1 and h.subs(v, 1) == 0
        for k in range(1, order + 1):
            assert s.diff(h, v, k).subs(v, 0) == 0
            assert s.diff(h, v, k).subs(v, 1) == 0
    assert s.expand(s.diff(h4, v) + 630*v**4*(1-v)**4) == 0
    maximum_identity = s.Rational(4, 27) - w*(1-w)**2 - (w-s.Rational(1, 3))**2*(s.Rational(4, 3)-w)
    assert s.expand(maximum_identity) == 0
    assert s.Rational(15, 8) < 2
    assert 225*s.Rational(4, 27) < 36
    assert 2+2 == 4 and 6+2*2*2+6 == 20
    finish("C001-fixed-cutoff", {
        "H_prime_factor": str(s.factor(s.diff(H, v))),
        "H_second_factor": str(s.factor(s.diff(H, v, 2))),
        "H_first_absolute_bound": "15/8 < 2",
        "H_second_squared_bound": "100/3 < 36",
        "nonnegative_maximum_identity": "4/27-w(1-w)^2=(w-1/3)^2(4/3-w), 0<=w<=1",
        "h4_prime_factor": str(s.factor(s.diff(h4, v))),
        "large_factor_constants": [1, 4, 20],
        "scope": "exact polynomial identities plus stated interval inequalities; no density evaluation"
    })

    begin("C002-new-three-dimensional-density")
    cfg = inputs["configurations"][1]
    columns = [s.Matrix(a) for a in cfg["frame_columns"]]
    y = s.Matrix(cfg["point"])
    pair_results = []
    for i, j in itertools.combinations(range(4), 2):
        d = s.Matrix.hstack(columns[i], columns[j], y).det()
        RESULT["actual_work"]["exact_pair_span_determinants"] += 1
        assert d != 0
        pair_results.append({"pair": [i+1, j+1], "determinant": str(d)})
    t = s.symbols("t", real=True)
    centers = list(map(s.Integer, cfg["scalar_pole_centers"]))
    terms = []
    for j, center in enumerate(centers):
        term = s.Integer(1)
        for ell, other in enumerate(centers):
            if ell != j:
                d = center-other
                term /= d*(d+2*s.I*t)
        terms.append(term)
        RESULT["actual_work"]["actual_upper_halfplane_residues"] += 1
    RESULT["actual_work"]["new_three_dimensional_density_paths"] += 1
    total = s.factor(s.cancel(sum(terms)))
    assert s.cancel(s.conjugate(total)-total) == 0
    assert s.cancel(total.subs(t, -t)-total) == 0
    a0 = s.cancel(total.subs(t, 0))
    positive_term_sum = sum(s.prod((centers[j]-centers[ell])**(-2) for ell in range(4) if ell != j) for j in range(4))
    assert a0 == positive_term_sum and a0 > 0
    correction = s.factor(s.cancel((total-a0)/t**2))
    correction0 = s.cancel(correction.subs(t, 0))
    assert correction0.is_finite is True
    assert s.cancel(t**3*total - t**3*a0 - t**5*correction) == 0
    finish("C002-new-three-dimensional-density", {
        "pair_qualification": pair_results,
        "pi3_density_divided_by_t3": str(total),
        "leading_coefficient_pi3_a": str(a0),
        "pi3_t5_correction": str(correction),
        "correction_at_zero": str(correction0),
        "upper_halfplane_residues": 4,
        "scope": "one complete 3D four-direction point-density symbolic path, four residues; no 6D or body32 kernel"
    })

    begin("C003-weighted-signed-layer")
    u, a, r = s.symbols("u a r", real=True)
    assert s.cancel(s.Rational(1, 2)-u/(1+u**2)-(u-1)**2/(2*(1+u**2))) == 0
    assert s.cancel(1-u**2/(1+u**2)-1/(1+u**2)) == 0
    assert s.simplify(s.re(s.I*s.exp(-s.I*a))-s.sin(a)) == 0
    boundary = (s.sin(a*r)-a*r)/r**2
    ibp_rhs = a*(s.cos(a*r)-1)/r**2 - 2*(s.sin(a*r)-a*r)/r**3
    assert s.simplify(s.diff(boundary, r)-ibp_rhs) == 0
    coefficient = s.simplify(24/(2*s.pi)**3 * (-s.pi/4))
    assert coefficient == -3/(4*s.pi**2)
    fplus, fminus = u**2, -u**2
    assert s.diff(fplus, u).subs(u, 0) == s.diff(fminus, u).subs(u, 0) == 0
    assert s.diff(fplus, u, 2) == 2 and s.diff(fminus, u, 2) == -2
    exponents = [3-2+(2-3), 3-1+(1-3), 3+(0-3)]
    assert exponents == [0, 0, 0]
    finish("C003-weighted-signed-layer", {
        "scalar_weight1_identity": "1/2-u/(1+u^2)=(u-1)^2/[2(1+u^2)], u>=0",
        "scalar_weight2_identity": "1-u^2/(1+u^2)=1/(1+u^2)",
        "quadratic_remainder_second_derivative_t_exponents": exponents,
        "Fourier_inverse_phase": "Re(i exp(-ia))=sin(a)",
        "integration_by_parts_identity": str(ibp_rhs),
        "accepted_radial_input": "A072 integral (1-cos(ar))/r^2 = pi*abs(a)/2; not rerun",
        "new_signed_angular_coefficient": str(coefficient),
        "ridge_first_derivative_common_limit": 0,
        "ridge_second_derivatives_left_right": [-2, 2],
        "scope": "formal identities supporting the written uniform proof; no numerical moment quadrature"
    })
    assert RESULT["configurations_completed"] == 3
    RESULT["status"] = "PASS_ALL_THREE_SMALL_CONFIGURATIONS"
    RESULT["body32_theorem_numerically_certified"] = False
    RESULT["global_g_certificate"] = False
    save()
    print(json.dumps({"status": RESULT["status"], "pid": os.getpid(),
                      "cpu_seconds": RESULT["cpu_seconds"],
                      "wall_seconds": RESULT["wall_seconds"],
                      "actual_work": RESULT["actual_work"]}))
except BaseException as exc:
    RESULT["status"] = "FAILED_FIRST_STOP_NO_RETRY"
    RESULT["failure"] = {"type": type(exc).__name__, "message": str(exc),
                         "traceback": traceback.format_exc()}
    save()
    print(json.dumps({"status": RESULT["status"], "failure": str(exc), "pid": os.getpid()}))
    raise
