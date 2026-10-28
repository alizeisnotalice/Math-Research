#!/usr/bin/env python3
"""Bounded exact intervals and uniform constants; no parameter or dimension grid."""
from fractions import Fraction as F
from hashlib import sha256
from math import factorial
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent

def log20(y):
    assert y > 1
    z = (y-1)/(y+1)
    lo = 2*sum((z**(2*j+1)/F(2*j+1) for j in range(20)), F(0))
    return lo, lo+2*z**41/(41*(1-z*z))

# Odd/even partial sums enclose exp(-1/2), by the alternating-series sign.
alo = sum((F(-1, 2)**k/F(factorial(k)) for k in range(22)), F(0))
ahi = sum((F(-1, 2)**k/F(factorial(k)) for k in range(21)), F(0))
slo, shi = 1/(2*ahi-1), 1/(2*alo-1)
blo, bhi = 3-4*ahi, 3-4*alo
clo, _ = log20(1/bhi)
_, chi = log20(1/blo)
lg5lo, _ = log20(F(5, 3))
_, lg4hi = log20(F(8, 5))
A_prime = 2*F(38, 15)*5+F(2, 3)*(5+6)+1
J_prime = (A_prime+17*5)/2
J1_prime = F(9, 5)+5*F(3, 5)*(F(15, 2)+F(25, 3))+9
checks = {
    "uniform_ell_floor": (1+F(1, 5))/2 == F(3, 5),
    "uniform_Q_bound": 1+F(3, 5)+F(3, 5)**2/2 > F(5, 3),
    "uniform_Q_derivative_bounds": F(4, 5)/(2-F(4, 5)) == F(2, 3) and 1/(2-F(4, 5)) == F(5, 6),
    "uniform_J_derivative_arithmetic": A_prime == F(101, 3) and J_prime == F(178, 3) < 60 and 2*60 == 120 and F(17, 2) < 9,
    "uniform_J1_derivative_arithmetic": J1_prime == F(583, 10) < 60 and 5*F(3, 5)**2 == F(9, 5) < 2,
    "uniform_J2_bound": 5*2*F(3, 5)**3/3 == F(18, 25) < 1,
    "uniform_log_moment_rounding": 9**2 == 81 and 4*9**2 == 324 and 2*9 == 18,
    "uniform_count_at_error_limit": F(1, 5)-F(6, 5)*F(1, 12) == F(1, 10) > 0,
    "uniform_face_variance_rounding": 5/(1-F(1, 12)) < 10,
    "uniform_transport_constants": 2*120 == 240 and 2*60 == 120,
    "exp_half_elementary_upper": F(13, 8)+F(1, 48)/(1-F(1, 8)) == F(277, 168) < F(5, 3),
    "exp_half_elementary_lower": F(13, 8) > F(8, 5),
    "critical_exp_interval": F(3, 5) < alo < ahi < F(5, 8),
    "critical_s_inside_family": 4 < slo < shi < 5,
    "derivative_endpoint_signs": lg4hi < F(1, 2) < lg5lo,
    "critical_log_positive_domain": 0 < blo < bhi < 1,
    "critical_s_decimal_enclosure": F(469348449, 10**8) < slo < shi < F(469348450, 10**8),
    "critical_c_decimal_enclosure": F(55533956, 10**8) < clo < chi < F(55533957, 10**8),
}
assert all(checks.values()), checks
values = dict(exp_minus_half_lower=alo, exp_minus_half_upper=ahi,
              sstar_lower=slo, sstar_upper=shi,
              cstar_lower=clo, cstar_upper=chi,
              uniform_J_prime_bound=J_prime, uniform_J1_prime_bound=J1_prime)
receipt = {
    "status": "PASS", "registration_sha256": sha256((HERE/"REGISTRATION.json").read_bytes()).hexdigest(),
    "scope": "18 exact constant/interval checks only; no grid, finite threshold or operator acceptance",
    "checks": checks,
    "exact_values": {k: str(v) for k, v in values.items()},
    "display_only": {k: float(v) for k, v in values.items()},
}
(HERE/"verification.json").write_text(json.dumps(receipt, indent=2)+"\n")
print(json.dumps({"status": "PASS", "exact_checks": len(checks),
                  "sstar_interval": [float(slo), float(shi)],
                  "cstar_interval": [float(clo), float(chi)]}))
