"""One-shot exact scalar guard; no original spatial/history simulation."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
REGISTRY = HERE / "actual_fullfuture_spatial_dual_registration_20261007.json"
RESULT = HERE / "actual_fullfuture_spatial_dual_parameter_results_20261007.json"
registration = json.loads(REGISTRY.read_text(encoding="utf-8"))
if RESULT.exists():
    raise RuntimeError("One-shot guard already has a receipt; do not overwrite it.")

checks = []


def check(label, condition):
    if not condition:
        raise AssertionError(label)
    checks.append(label)


def exact(q):
    return {"numerator": str(q.numerator), "denominator": str(q.denominator)}


rows = []
previous = None
for case in registration["rounds"]:
    n, m, terms = case["n"], case["m"], case["atanh_terms"]
    z = F(1, 2)
    lo = 2 * sum((z ** (2*j+1) / (2*j+1) for j in range(terms)), F(0))
    tail = 2 * z ** (2*terms+1) / ((2*terms+1) * (1-z*z))
    hi = lo + tail
    ratio_lo = 4*lo - 2*hi*hi/n
    ratio_hi = 4*hi
    sigma_lo, sigma_hi = F(m, n)*ratio_lo, F(m, n)*ratio_hi
    error_hi = 2*hi*hi/n
    width = ratio_hi-ratio_lo
    prefix = f"n={n}"
    check(prefix+": sixth-power dimension", n == m**6)
    check(prefix+": log3 positive enclosure", F(1) < lo < hi < F(3, 2))
    check(prefix+": exponential argument below one", hi/n < 1)
    check(prefix+": nonempty ratio interval", 0 < ratio_lo < ratio_hi)
    check(prefix+": sigma strictly inside (0,1)", 0 < sigma_lo < sigma_hi < 1)
    check(prefix+": n sigma strictly exceeds m+1", m*ratio_lo > m+1)
    check(prefix+": exact enclosure-width identity", width == 4*tail+error_hi)
    check(prefix+": vanishing error envelope below 3/n", 0 < error_hi < F(3, n))
    if previous is not None:
        check(prefix+": ratio enclosure narrows", width < previous["width"])
        check(prefix+": asymptotic error bound decreases", error_hi < previous["error_hi"])
    rows.append({
        **case,
        "log3_lower": exact(lo), "log3_upper": exact(hi),
        "ratio_lower": exact(ratio_lo), "ratio_upper": exact(ratio_hi),
        "sigma_lower": exact(sigma_lo), "sigma_upper": exact(sigma_hi),
        "asymptotic_error_upper": exact(error_hi),
        "simple_error_upper": exact(F(3, n)),
        "enclosure_width": exact(width),
        "actual_ledger_dimension": n >= 512
    })
    previous = {"width": width, "error_hi": error_hi}

if len(checks) != registration["expected_checks"]:
    raise AssertionError("Unexpected assertion count")
result = {
    "id": registration["id"], "status": "PASS_EXACT_SCALAR_PARAMETER_GUARD",
    "checks_passed": len(checks), "rounds": rows, "checks": checks,
    "scope": registration["scope"], "ledger": registration["ledger"],
    "registration_sha256": hashlib.sha256(REGISTRY.read_bytes()).hexdigest(),
    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "runs": 1, "random_seed": None,
    "spatial_integrals_computed": False, "actual_gate_samples": False
}
RESULT.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
print(json.dumps({"status": result["status"], "checks": len(checks),
                  "dimensions": [r["n"] for r in rows],
                  "result": str(RESULT)}, ensure_ascii=False))
