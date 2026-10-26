"""Finite exact checks for the new anchor-time curvature formula only.

No kernel, H, gradient, log, old baseline, or transcendental evaluation.
This file refuses to overwrite a completed receipt.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import time

import sympy as sp

HERE = Path(__file__).resolve().parent
REG = HERE / "registration.json"
OUT = HERE / "results.json"
if OUT.exists():
    raise SystemExit("Completed receipt exists; no rerun authorized by this entry point.")
started = time.monotonic()
reg_raw = REG.read_bytes()
reg = json.loads(reg_raw)
script_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
assert reg["script_sha256"] == script_sha
for source in reg["sources"]:
    assert hashlib.sha256(Path(source["path"]).read_bytes()).hexdigest() == source["sha256"]

checks = {}
H0, H1, L, U = sp.symbols("H0 H1 L U", positive=True)
a, b = sp.symbols("a b", positive=True)
q_stair = (a * H0 - (a - b) * H1) / (16 * H0)
q_end = (a - (a - b) * L / U) / 16
endpoint_residual = (a - b) * ((H1 - L) * U + L * (U - H0)) / (16 * H0 * U)
assert sp.cancel(q_end - q_stair - endpoint_residual) == 0
checks["endpoint_direction_factorization"] = str(endpoint_residual)

special = sp.cancel(q_end.subs({a: 128**2, b: 64**2}))
assert sp.cancel(special - 1024 * (1 - sp.Rational(3, 4) * L / U)) == 0
checks["epsilon128_epsilon64_coefficient"] = str(special)
assert sp.cancel(special - 256 - 768 * (U - L) / U) == 0
assert sp.cancel(1024 - special - 768 * L / U) == 0
checks["positive_lower_and_strict_improvement"] = {
    "qU_minus_256": str(768 * (U - L) / U),
    "1024_minus_qU": str(768 * L / U),
}

# Three cutoffs, represented by the nonnegative masses in the three bands.
# H_j is the sum of all masses at/after cutoff j, including the infinite tail.
a0, a1, a2 = sp.symbols("a0 a1 a2", positive=True)
m0, m1, m2 = sp.symbols("m0 m1 m2", nonnegative=True)
tail_expression = a0 * (m0 + m1 + m2) - (a0-a1) * (m1+m2) - (a1-a2) * m2
assert sp.expand(tail_expression - (a0*m0 + a1*m1 + a2*m2)) == 0
checks["multicutoff_telescoping"] = str(sp.expand(tail_expression))

receipt = {
    "status": "passed_finite_new_symbolic_identities_only",
    "pid": os.getpid(),
    "runtime_seconds": time.monotonic() - started,
    "python": os.sys.version,
    "sympy": sp.__version__,
    "registration_sha256": hashlib.sha256(reg_raw).hexdigest(),
    "script_sha256": script_sha,
    "checks": checks,
    "kernel_evaluations": 0,
    "H_evaluations": 0,
    "gradient_evaluations": 0,
    "log_evaluations": 0,
    "old_baseline_runs": 0,
    "matlab_calls": 0,
}
OUT.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n")
print(json.dumps(receipt, ensure_ascii=False, indent=2))
