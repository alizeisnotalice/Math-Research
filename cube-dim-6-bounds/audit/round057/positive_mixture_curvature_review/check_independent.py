#!/usr/bin/env python3
"""Fixed A057 independent arithmetic audit; synthetic jets only, no kernels."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import copy
import hashlib
import importlib.util
import json
import resource
import signal
import sys
import time

resource.setrlimit(resource.RLIMIT_CPU, (30, 30))
signal.alarm(900)
START_CPU = time.process_time()
START_WALL = time.monotonic()

import sympy as s
import flint

ROOT = Path(__file__).resolve().parent
REGISTRATION = json.loads((ROOT / "input_manifest.json").read_text())
OUT = ROOT / "results.json"
if OUT.exists():
    raise RuntimeError("Completed audit batch cannot be overwritten")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


for path, expected in REGISTRATION["pre_execution_sha256"].items():
    assert digest(path) == expected, path
SOURCE = Path(REGISTRATION["author_directory"])
spec = importlib.util.spec_from_file_location("a057_source_box", SOURCE / "box_certificate.py")
source = importlib.util.module_from_spec(spec)
spec.loader.exec_module(source)

checks = []
call_count = 0


def call(data, precision=64):
    global call_count
    call_count += 1
    return source.certify(data, precision_bits=precision)


def rejected(base, cases):
    rows = []
    for name, mutate in cases:
        trial = copy.deepcopy(base)
        mutate(trial)
        try:
            call(trial)
        except ValueError as exc:
            rows.append({"case": name, "rejected": True, "message": str(exc)})
        else:
            raise AssertionError(f"Invalid input accepted: {name}")
    return rows


def exp_rational_bracket(value):
    """Exact positive Taylor sum and geometric tail; reciprocal for negative x."""
    x = abs(value)
    n = 200
    assert x < n + 2
    term = F(1)
    partial = term
    for k in range(1, n + 1):
        term *= x / k
        partial += term
    next_term = term * x / (n + 1)
    upper = partial + next_term / (1 - x / (n + 2))
    if value >= 0:
        return partial, upper
    return 1 / upper, 1 / partial


# 1: distinct star basis, affine fiber reconstruction, and one null perturbation.
star = [(1, 1, 1, 1, 1, 1)]
for j in range(1, 6):
    v = [1] * 6
    v[j] = -1
    star.append(tuple(v))
remaining = [v for tail in product((-1, 1), repeat=5)
             if (v := (1,) + tail) not in star]
N = s.Matrix.hstack(*map(s.Matrix, star))
M = s.Matrix.hstack(*map(s.Matrix, remaining))
A = N.row_join(M)
R = A.T / 32
K = (-N.inv() * M).col_join(s.eye(26))
assert N.det() == -32
assert A * K == s.zeros(6, 26)
# For all x,w, N^-1(x-M(R_M x+w)) = R_N x + K_N w.
assert N.inv() * (s.eye(6) - M * R[6:, :]) == R[:6, :]
assert -N.inv() * M == K[:6, :]
E = s.zeros(26, 6)
E[0, 0], E[3, 2] = s.Rational(1, 7), s.Rational(-2, 9)
D = K * E
assert D != s.zeros(32, 6)
assert A * (R + D) == s.eye(6)
assert (R + D).T * (R + D) - R.T * R == D.T * D
checks.append({"configuration": 1, "star_basis_determinant": str(N.det()),
               "fiber_affine_identity": True, "AK_zero": True,
               "nonzero_nullspace_perturbation_decomposition": True, "pass": True})

# 2: quotient/variance algebra with arbitrary positive component masses.
w = s.symbols("w0:3", positive=True)
score = s.symbols("score0:3", real=True)
tau = s.symbols("tau0:3", real=True)
mass = sum(w)
mean = sum(wi * si for wi, si in zip(w, score)) / mass
quotient = sum(wi * (ti + si ** 2) for wi, ti, si in zip(w, tau, score)) / mass - mean ** 2
variance_form = sum(wi * ti for wi, ti in zip(w, tau)) / mass + sum(
    wi * (si - mean) ** 2 for wi, si in zip(w, score)) / mass
assert s.factor(quotient - variance_form) == 0
checks.append({"configuration": 2, "three_mass_quotient_variance_identity": True,
               "continuous_mixture_proof_is_analytic": True, "pass": True})

window = {"kind": "B", "a": "3/2", "T": "5"}
center = ["1", "0", "0", "0", "0", "0"]
base = {
    "model": source.MODEL, "model_definition": source.MODEL_DEFINITION,
    "window": copy.deepcopy(window), "center": list(center),
    "face_box": [["-1", "12"], ["0", "7"], ["0", "0"], ["0", "0"], ["0", "0"]],
    "quadratic": "frame32", "target": "1",
    "jet": {"model": source.MODEL, "model_definition": source.MODEL_DEFINITION,
            "window": copy.deepcopy(window), "center": list(center),
            "F_interval": ["2", "5"],
            "face_gradient_intervals": [["1/4", "1/3"], ["1", "2"],
                                        ["0", "0"], ["0", "0"], ["0", "0"]],
            "source_proof": "SYNTHETIC_AUDIT_ONLY_NOT_A_BODY32_JET"}}


def verify_corners(result, expected):
    distinct = {}
    for row in result["corners"]:
        delta = tuple(map(F, row["displacement"]))
        key = (delta[0], delta[1])
        dot, f, exponent = map(F, expected[key])
        assert F(row["gradient_dot_lower"]) == dot
        assert F(row["minimizing_center_value"]) == f
        assert F(row["exponent"]) == exponent
        el, eu = exp_rational_bracket(exponent)
        lo, hi = map(F, [row["certificate_lower"], row["certificate_upper"]])
        assert lo <= f * el <= f * eu <= hi
        distinct[key] = {"displacement_2d": list(map(str, key)),
                         "dot": str(dot), "minimizing_F": str(f), "exponent": str(exponent),
                         "arb_lower": str(lo), "arb_upper": str(hi),
                         "contains_independent_exact_Taylor_bracket": True}
    assert len(result["corners"]) == 32 and len(distinct) == 4
    assert F(result["tangent_bound_lower"]) == min(F(row["certificate_lower"]) for row in result["corners"])
    assert F(result["tangent_bound_upper"]) == min(F(row["certificate_upper"]) for row in result["corners"])
    return list(distinct.values())


# 3: asymmetric band, all three clamp branches; independent rational exp proof.
frame_expected = {
    (F(-1), F(0)): ("-1/3", "2", "-13/72"),
    (F(12), F(0)): ("3", "3", "-1"),
    (F(-1), F(7)): ("20/3", "5", "23/36"),
    (F(12), F(7)): ("10", "5", "-49/72"),
}
frame = call(base)
assert frame["status"] == "PASS_CONDITIONAL_ON_MATCHING_CENTER_JET"
frame_rows = verify_corners(frame, frame_expected)
assert F(frame["tangent_bound_upper"]) < F(1104, 1000)
assert F(frame["tangent_bound_lower"]) > F(1103, 1000)
checks.append({"configuration": 3, "synthetic_only": True, "status": frame["status"],
               "unique_corners": frame_rows, "winning_clamp_is_interior": True,
               "tangent_bound_lower": frame["tangent_bound_lower"],
               "tangent_bound_upper": frame["tangent_bound_upper"], "pass": True})

# 4: same independent input, alternate proved fixed quadratic, no source Q calls.
basis = copy.deepcopy(base)
basis["quadratic"], basis["target"] = "basis6", "1/10"
basis_expected = {
    (F(-1), F(0)): ("-1/3", "2", "-19/72"),
    (F(12), F(0)): ("3", "3", "-13"),
    (F(-1), F(7)): ("20/3", "5", "-67/18"),
    (F(12), F(7)): ("10", "5", "-1039/72"),
}
basis_result = call(basis)
assert basis_result["status"] == "UNRESOLVED"
basis_rows = verify_corners(basis_result, basis_expected)
checks.append({"configuration": 4, "synthetic_only": True, "status": basis_result["status"],
               "unique_corners": basis_rows, "pass": True})

# 5: eight fixed inexact/bool input cases.
inexact = [
    ("float_root_a", lambda d: d["window"].update(a=1.5)),
    ("float_root_T", lambda d: d["window"].update(T=5.0)),
    ("float_center", lambda d: d["center"].__setitem__(1, 0.0)),
    ("float_box", lambda d: d["face_box"][0].__setitem__(0, -1.0)),
    ("float_F", lambda d: d["jet"]["F_interval"].__setitem__(0, 2.0)),
    ("float_gradient", lambda d: d["jet"]["face_gradient_intervals"][0].__setitem__(0, 0.25)),
    ("float_target", lambda d: d.update(target=1.0)),
    ("boolean_target", lambda d: d.update(target=True)),
]
checks.append({"configuration": 5, "cases": rejected(base, inexact), "pass": True})

# 6: eight fixed metadata cases.
metadata = [
    ("jet_window_kind", lambda d: d["jet"]["window"].update(kind="H")),
    ("jet_window_a", lambda d: d["jet"]["window"].update(a="2")),
    ("jet_center", lambda d: d["jet"]["center"].__setitem__(2, "1/7")),
    ("jet_model", lambda d: d["jet"].update(model="sparse4")),
    ("jet_normalization", lambda d: d["jet"].update(model_definition="32 factors with different weights")),
    ("empty_provenance", lambda d: d["jet"].update(source_proof="")),
    ("nonincreasing_band", lambda d: d["window"].update(T="3/2")),
    ("H_with_finite_T", lambda d: d["window"].update(kind="H")),
]
checks.append({"configuration": 6, "cases": rejected(base, metadata), "pass": True})

# 7: exact equality and just-above-target behavior at a point box.
point = copy.deepcopy(base)
point["face_box"] = [["0", "0"] for _ in range(5)]
point["target"] = "2"
point_pass = call(point)
assert point_pass["status"] == "PASS_CONDITIONAL_ON_MATCHING_CENTER_JET"
assert point_pass["tangent_bound_lower"] == point_pass["tangent_bound_upper"] == "2"
point["target"] = "20000000000000000001/10000000000000000000"
point_unresolved = call(point)
assert point_unresolved["status"] == "UNRESOLVED"
checks.append({"configuration": 7, "synthetic_only": True, "exact_point_bound": "2",
               "target_equal": point_pass["status"], "target_above": point_unresolved["status"],
               "pass": True})

# 8: four fixed geometric/shape cases.
geometry = [
    ("reversed_face_interval", lambda d: d["face_box"].__setitem__(0, ["12", "-1"])),
    ("center_outside", lambda d: d["face_box"].__setitem__(0, ["1", "12"])),
    ("four_gradients", lambda d: d["jet"]["face_gradient_intervals"].pop()),
    ("unsupported_quadratic", lambda d: d.update(quadratic="pointwise_min_bases")),
]
checks.append({"configuration": 8, "cases": rejected(base, geometry), "pass": True})

assert len(checks) == 8 and call_count == 24
for path, expected in REGISTRATION["pre_execution_sha256"].items():
    assert digest(path) == expected, path
result = {
    "status": "PASS_EIGHT_INDEPENDENT_AUDIT_CONFIGURATIONS",
    "configurations": 8, "checks": checks, "source_certify_calls": call_count,
    "accepted_input_interface_calls": 4, "rejected_input_interface_calls": 20,
    "production_G_H_p_gradient_Hessian_kernel_calls": 0, "old_formal_calls": 0,
    "author_fixed_batch_runs": 0, "matlab_calls": 0, "fullfan_calls": 0,
    "synthetic_jets_are_kernel_evidence": False, "global_g_certificate": False,
    "new_C6_bound": False, "precision_bits": 64, "Taylor_terms_through_degree": 200,
    "python": sys.version, "sympy": s.__version__, "python_flint": flint.__version__,
    "cpu_seconds": time.process_time() - START_CPU,
    "wall_seconds": time.monotonic() - START_WALL,
    "pre_and_post_source_hashes_match": True,
    "pre_execution_sha256": REGISTRATION["pre_execution_sha256"],
}
OUT.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({k: result[k] for k in ["status", "configurations", "source_certify_calls",
                                       "production_G_H_p_gradient_Hessian_kernel_calls",
                                       "cpu_seconds", "wall_seconds"]}, indent=2))
