"""Independent narrow arithmetic preflight; not a proof of the operator claim."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import hashlib, json, os

L = Q(8192)
x = Q(91, 100) * L
h = 1 / (256 * L)
trim = Q(1, 2**19)
eN = Q(1, 2**19)
eJ = Q(1, 2**20)
rlo, rhi = 1 + trim, Q(273, 200) - trim
# A different n lower envelope from the author's 64-term Taylor gate.
nlo = L**32 / factorial(32)
p = 1 - Q(1, 2**19) - Q(3, 128)
gap = Q(9, 200) - Q(2211, 50000)/(1-eJ) - Q(183, 1000)/x
v = Q(531, 100000)
ulo = p*(1-eN)/(1+eJ)
uhi = (1+eN)/(p*(1-eJ))
M = 2*(1+eJ)*uhi
b = 2*Q(11,25)*(1-eJ)*ulo-eJ-Q(1,2**19)
Gamma = Q(311,1000)*b/M
checks = {
    "trim_over_h_is_four": trim/h == 4,
    "matching_baseline_and_active_domain": trim > 2*h and rlo-2*h > 1,
    "count_drift_with_independent_n_envelope": 2*h+1/nlo < 3*eN/4,
    "sharp_logproduct_shift": 6*x*h < Q(3,128),
    "first_profile_even_old_derivative_bridge": 16*h < Q(1,100000),
    "second_profile_coarse_bridge": 2*h < Q(1,1000000),
    "pointwise_J": x/(2*nlo) < eJ,
    "radial_width": 1/(1-(h/rlo)**2) < 1+eJ,
    "positive_gap": gap > Q(1,1330),
    "cantelli_probability_direct": v/(v+x*gap*gap) < Q(14,25),
    "cantelli_conservative_gap": v/(v+x/Q(1330**2)) < Q(14,25),
    "conditional_energy_budget": (1-eJ)*Gamma > Q(13,100),
}
assert all(checks.values()), checks
out = {"status":"PASS", "pid":os.getpid(), "checks":checks,
       "scope":"tuple arithmetic and decreasing-width preflight only; no frozen full proof yet",
       "gap":str(gap), "conservative_tail":str(v/(v+x/Q(1330**2))),
       "conditional_Gamma":str(Gamma),
       "infinite_tail_note":"For L>=8192, h and 1/(L^32/32!) decrease; x grows; xh is constant. Every displayed adverse margin is therefore bounded by its base value.",
       "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name("preflight.json").write_text(json.dumps(out, indent=2)+"\n")
print(json.dumps({"status":"PASS", "pid":out["pid"], "checks":len(checks)}))
