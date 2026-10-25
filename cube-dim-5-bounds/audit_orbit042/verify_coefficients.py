#!/usr/bin/env python3
"""One exact, finite coefficient audit; no SPD or Gamma-budget computation."""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import signal
import sys
import time

import sympy as sp

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--output-dir", type=Path, default=HERE)
args = parser.parse_args()
OUT = args.output_dir.resolve()
OUT.mkdir(parents=True, exist_ok=True)
started = time.monotonic()

def checkpoint(stage):
    data = {"pid": os.getpid(), "stage": stage,
            "elapsed_seconds": time.monotonic() - started,
            "deadline_seconds": 120,
            "recoverable_command": f"OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 {sys.executable} -B {__file__} --output-dir {OUT}"}
    (OUT / "checkpoint.json").write_text(json.dumps(data, indent=2) + "\n")

def timed_out(_signum, _frame):
    checkpoint("timeout_no_certificate")
    raise TimeoutError("Registered 120-second recoverable deadline")

signal.signal(signal.SIGALRM, timed_out)
signal.alarm(120)
checkpoint("source_hash")
reg = json.loads((HERE / "registration.json").read_text())
source = Path(reg["source"]["path"])
actual_source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
assert actual_source_hash == reg["source"]["sha256"]
checks = {}

def identity(name, lhs, rhs):
    residual = sp.simplify(lhs - rhs)
    assert residual == 0, (name, residual)
    checks[name] = {"lhs": str(lhs), "rhs": str(rhs), "residual": "0"}

checkpoint("pair_and_holder")
q, t = sp.symbols("q t", positive=True)
identity("log_sinh_Hessian", sp.diff(sp.log(sp.sinh(q*t)), q, 2),
         -t**2 / sp.sinh(q*t)**2)
F = sp.cosh(q*t) - 1 - 2**(q-1)*(sp.cosh(t)**q - 1)
identity("pair_gap_derivative", sp.diff(F, t),
         q*(sp.sinh(q*t) - 2**(q-1)*sp.cosh(t)**(q-1)*sp.sinh(t)))
a, b = 2 - 2/q, 2/q - 1
identity("Holder_weights_sum", a+b, 1)
identity("Holder_moment_exponent", a*q/2+b*q, 1)
identity("Holder_after_raising_first", q*a, 2*q-2)
identity("Holder_after_raising_second", q*b, 2-q)
identity("pair_normalization", 2/2**q + (2**q-2)/2**q, 1)
identity("pair_gap_constant", (2**q-2)/2**q, 1-2**(1-q))
Q = sp.Rational(5,4)
identity("sphere_face_coefficient_without_pi_squared", sp.Rational(8,3)/160, sp.Rational(1,60))
identity("sphere_homogeneity", -4*Q, -5)
identity("determinant_norm_exponent", sp.Rational(-1,8)/Q, sp.Rational(-1,10))
identity("orbit_norm_energy_exponent", (Q-1)/Q, sp.Rational(1,5))

checkpoint("Schur_and_Beta")
s, k, z, y = sp.symbols("s k z y", nonnegative=True)
identity("covariance_precision_parameter", 1/(1-s)-1, s/(1-s))
identity("single_flip_Schur_ratio", (1+k-z*z*k), 1+(1-z*z)*k)
identity("single_flip_u", 1-(2*t-1)**2, 4*t*(1-t))
identity("single_gap_Hessian", sp.diff(1-(1+k*y)**sp.Rational(-1,2), y, 2),
         -sp.Rational(3,4)*k*k*(1+k*y)**sp.Rational(-5,2))
r1, r2, u = sp.symbols("tau1 tau2 u")
kap1, kap2 = r1/(1-r1), r2/(1-r2)
identity("split_rank_two_determinant_ratio",
         ((1-(1-u)*r1)*(1-(1-u)*r2))/((1-r1)*(1-r2)),
         (1+u*kap1)*(1+u*kap2))
identity("split_D", (1+kap1)*(1+kap2), 1/((1-r1)*(1-r2)))
identity("split_K_trace", 1/(1-r1)+1/(1-r2)-2, kap1+kap2)
identity("split_product_minus_linear", (1+u*kap1)*(1+u*kap2)-1-u*(kap1+kap2), u*u*kap1*kap2)
identity("Beta_Y_mean", Q/(Q+sp.Rational(1,2)), sp.Rational(5,7))
identity("C2_from_Y_Beta_recurrence", (Q-sp.Rational(1,2))/(Q-1), 3)
identity("C2_from_T_Beta_recurrence", sp.Rational(1,4)*16*sp.Rational(3,2)*sp.Rational(1,2), 3)
checks["Beta_rank_formula"] = {"T_form": "C_k=4^(-k/2) B(q-k/2,q-k/2)/B(q,q)",
    "Y_form": "C_k=B(q-k/2,1/2)/B(q,1/2)",
    "integrability_condition": "k < 2q=5/2; ranks 0,1,2 are finite"}

checkpoint("tilted_moments_and_coarse_coefficients")
beta, N = sp.symbols("beta N", positive=True)
cross = 2*beta/(N*N*(1+beta))
extra = (1-beta)/(N*(1+beta))
diag = cross+extra
identity("tilted_cross_Gamma_cancel", (2*beta/N**3)*(N/(1+beta)), cross)
identity("tilted_diag_extra_Gamma_cancel", ((1-beta)/N**2)*(N/(1+beta)), extra)
identity("mass_second_moment", N*diag+N*(N-1)*cross, 1)
identity("Cauchy_cross", cross.subs({beta:sp.Rational(1,2), N:16}), sp.Rational(1,384))
identity("Cauchy_diagonal", diag.subs({beta:sp.Rational(1,2), N:16}), sp.Rational(3,128))
V = [(1,)+tail for tail in itertools.product((-1,1), repeat=4)]
balances = [sum(v[i]*v[j] for v in V) for i in range(5) for j in range(i+1,5)]
assert len(V) == 16 and len(balances) == 10 and balances == [0]*10
checks["body_pair_sign_balance"] = {"directions": len(V), "coordinate_pairs": 10, "signed_sums": balances}
identity("E_Mij_squared", (N*(diag-cross)).subs(beta, sp.Rational(1,2)), sp.Rational(1,3))
identity("E_w", 10*sp.Rational(1,3), sp.Rational(10,3))
identity("five_flip_q_gap", sp.Rational(1,5)*sp.Rational(5,14)*sp.Rational(1,4)*2, sp.Rational(1,28))
identity("norm_gap", sp.Rational(1,28)/Q, sp.Rational(1,35))
identity("J_subtraction", sp.Rational(1,35)*sp.Rational(10,3), sp.Rational(2,21))
assert sp.Rational(41,42)**2 > sp.Rational(19,20)
checks["coarse_relative_budget_floor"] = {"rational_residual": str(sp.Rational(41,42)**2-sp.Rational(19,20)), "positive": True}

results = {"status": "all_exact_checks_passed", "source_sha256": actual_source_hash,
           "symbolic_or_rational_checks": checks,
           "new_I_or_special_function_numerics": False,
           "scope": "Necessary identities only; report contains the analytic inequality and measure arguments."}
(OUT / "results.json").write_text(json.dumps(results, indent=2) + "\n")
signal.alarm(0)
checkpoint("complete")
receipt = {"pid": os.getpid(), "elapsed_seconds": time.monotonic()-started,
           "exit_code": 0, "python": sys.version.split()[0], "sympy": sp.__version__,
           "checks": len(checks), "threads": 1}
(OUT / "run_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
print(json.dumps(receipt))
