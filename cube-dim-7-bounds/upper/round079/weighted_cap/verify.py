#!/usr/bin/env python3
"""New LOCAL weighted-cap arithmetic, conditional on geometric guards."""
from fractions import Fraction as F
from hashlib import sha256
from math import factorial
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[1] / "round073/cantelli_frontier/PROOF.md"
SOURCE_HASH = "50d48b85f02a14f50614f0eea5cedbd7ccd68e9b2f7e1564fe8b41fcd2af9c4b"
assert sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_HASH
L0 = F(4096)
N0 = L0**64/F(factorial(64))
x0 = F(231, 256)*L0
eN = eB = F(1, 2**19)
eJ = eRad = F(1, 2**20)
qmax = F(11569, 12800)
eta = (1-qmax)/(2*qmax)
nu = F(2, 3)-F(5, 3)*eN
count_floor = F(333, 500)
cfixed = count_floor*eta/(1+eJ)
P0 = 1-eB-F(6, 768)
Uhi = (1+eN)/(P0*(1-eJ))
Mdegree = 2*(1+eRad)*Uhi
exp100_term4 = F(100)**4/F(factorial(4))
checks = {
    "source_sha256_matches": True,
    "positive_band_eta": 0 < qmax < 1 and eta == F(1231, 23138),
    "N0_exceeds_million": N0 > 10**6,
    "count_omission_budget": nu-F(1, 10**6) > count_floor > 0,
    "degree_bound_lt_three": 0 < Mdegree < 3,
    "new_cap_exponent_gt_100": cfixed*x0 > 100,
    "positive_exp_term_pays_eS": exp100_term4 > 3*2**20,
}
assert all(checks.values()), checks
values = dict(x0=x0, nu=nu, eta=eta, count_floor=count_floor,
              cfixed=cfixed, cap_exponent_at_4096=cfixed*x0,
              Mdegree=Mdegree, exp100_term4=exp100_term4)
receipt = {
    "status": "PASS", "scope": "7 new LOCAL cap gates only; no full C_n theorem",
    "conditions": "Existing full-matching count, J, geometry and radial-width gates must hold",
    "source_sha256": SOURCE_HASH,
    "checks": checks,
    "exact_values": {k: str(value) for k, value in values.items()},
    "display_only": {k: float(value) for k, value in values.items()},
}
(HERE/"verification.json").write_text(json.dumps(receipt, indent=2)+"\n")
print(json.dumps({"status": "PASS", "local_new_gates": len(checks),
                  "cap_exponent_lower": float(cfixed*x0),
                  "scope": "weighted cap only"}))
