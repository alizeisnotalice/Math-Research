"""Reproduce the finite Fraction check originally executed as a tool snippet.

This is a new run. It neither rewrites the historical JSON nor attributes a
process ID or timestamp to that earlier execution. It makes no new geometric
or asymptotic claim.
"""
from fractions import Fraction as F
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import sys

HERE = Path(__file__).resolve().parent
OLD = HERE / "micro_increment_fraction.json"
OUTPUT = HERE / "micro_increment_fraction_recheck.json"

# This core is the same finite rational calculation as the executed snippet.
k0 = F(16, 3)
k1 = -4*k0 + F(640, 81)
l0 = k0**3/32
l1 = (3*k0**2*k1 + 2*k0**3)/32 - 5*k0**4/768
f0dd = 2 - 3*k0/4
f1d = -8*k0 + 5*k0**3/48
assert f0dd*k1 + f1d == 0
micro8 = 16*l0/256
micro9 = (16*l1 - 96*l0 - 5*k0**4/48)/512
protect9 = -(5*k0**4/48)/512
assert micro8 == F(8, 27)
assert micro9 == -F(496, 243)
assert protect9 == -F(40, 243)
result = {
    "scope": "R073 finite rational normalization only; host and actual construction proved separately",
    "kappa0": str(k0),
    "kappa1": str(k1),
    "lambda_macro_constant": str(l0),
    "lambda_macro_linear": str(l1),
    "kappa2_multiplier_at_order9": str(f0dd*k1 + f1d),
    "micro_increment_av8": str(micro8),
    "micro_increment_av9": str(micro9),
    "new_first_protection_only_av9": str(protect9),
    "assertions": "passed",
}

old_bytes = OLD.read_bytes()
old_result = json.loads(old_bytes)
assert result == old_result, (result, old_result)
output = {
    "run_kind": "new reproducibility recheck; not the historical execution",
    "python": sys.version,
    "pid": os.getpid(),
    "executed_at_utc": datetime.now(timezone.utc).isoformat(),
    "input_script": Path(__file__).name,
    "input_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "historical_json": OLD.name,
    "historical_json_sha256": hashlib.sha256(old_bytes).hexdigest(),
    "historical_payload_matches_exactly": True,
    "result": result,
}
# Exclusive creation prevents an accidental rerun from overwriting this output.
with OUTPUT.open("x", encoding="utf-8") as stream:
    json.dump(output, stream, ensure_ascii=False, indent=2)
    stream.write("\n")
print(json.dumps(output, ensure_ascii=False, indent=2))
