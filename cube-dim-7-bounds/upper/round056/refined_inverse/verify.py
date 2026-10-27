"""Exact fixed-tuple budget and finite-gate checks; no search or simulation."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import json

L = Q(2**20)
d, sigma, chi = Q(3, 10), Q(7, 100), Q(7, 200)
x = (1 - sigma)*L
eps0, trim = Q(1, 2**24), Q(1, 2**20)
e, epsA = Q(1, 1024), Q(1, 64)
qmin, qmax = Q(31, 40), Q(37, 40)
rlo, rhi = Q(43, 40) + trim, Q(49, 40) - trim
h = eps0/(L*d)
N = L**8/factorial(8)
P0 = 1 - (Q(1, 256) + 1280*eps0)
assert P0 == Q(65275, 65536)
Ulo, Uhi = P0*(1-e)/(1+e), (1+e)/(P0*(1-e))
M = Q(21, 8)
zA = (rhi-rlo)/rlo
mass = (1-epsA)*(zA-zA**2/2+zA**3/3-zA**4/4)
rowgap = (1-e)*2*Ulo/qmax - 2*e
gamma = (1-e)*mass*rowgap/M
craw = gamma - sigma/(1-sigma)
u = Q(500, 3083)
Qpin = u-u*u/2+u**3/3
fac = Q(43, 31)/(1-e)
ML = 6*d**3*L**2/(qmin**2*N) + d*L/(qmin*N)
VL = 12*d**3*L**2/(qmin**2*N)
cap_bound = factorial(8)*d*d/(eps0*Q(93, 20000)**8*L**7)
tau_bound = factorial(8)*d/(eps0*Q(49, 20000)**8*L**6*Q(2**10))
scalar_bound = (Q(factorial(8))/(sigma*L)**8
                +2*Q(2**10)*factorial(8)/(chi*L)**8
                +Q(factorial(8))/L**8)
checks = {
    "domain": Q(0)<zA<1 and d<=Q(1,3) and 0<gamma<1,
    "n_lower": N>=75,
    "x_lower": x>=1000,
    "t_lower": N/L>=4,
    "source_radius": qmin*N/L>=4*rhi,
    "h_domain": h<=min(trim,d/8,Q(1,10000)),
    "zeta_guard": 1/(x*N*N)<=h/2,
    "g_R_domain": d<=(Q(43,40)-Q(1,10000))/2,
    "J_individual": 2*d*L/(qmin*N)<=e,
    "baseline_mean": ML<=Q(1,512),
    "baseline_variance": VL<=Q(1,2**26),
    "count_probability": 4/(e*e*N)<=epsA/4,
    "count_calibration": 3*h+1/N<=e/2,
    "target_cap_probability": Q(3,4)**74<=epsA/4,
    "target_cap_side": Q(9,8)-Q(20,49)<Q(3,4),
    "pin_log_bound": Qpin<Q(1,6),
    "mean_rate": fac*Qpin**2<Q(4,125),
    "variance_rate": Q(2,3)*fac*Qpin**3<Q(1,300),
    "pin_budget": Qpin/x<=Q(1,2000),
    "Chebyshev_gap": chi-Q(4,125)-Q(1,2000)==Q(1,400),
    "retention_probability": Q(1600,3)/x<=e,
    "row_max": 2*Uhi/qmin<M,
    "source_cap_exponent": Q(74,75)*Q(3,74)==Q(1,25),
    "source_cap_gap": Q(1,25)-chi==Q(1,200),
    "source_cap_payment": cap_bound<=e,
    "tau_payment": tau_bound<=e,
    "positive_energy": rowgap>0,
    "net_margin": craw>Q(1,40),
    "full_scalar_payment": scalar_bound<craw-Q(1,40),
}
assert all(checks.values()), [name for name, ok in checks.items() if not ok]
result = {
    "status": "PASS",
    "scope": "exact rational fixed-tuple and L0 bases; analytic proof supplies all-L extension and original-operator bridge",
    "L0": str(L),
    "checks": checks,
    "P0": str(P0), "Ulo": str(Ulo), "Uhi": str(Uhi),
    "Qpin_upper": str(Qpin), "mass_lower": str(mass),
    "gamma": str(gamma), "craw_lower": str(craw),
    "cap_ratio_upper": str(cap_bound), "tau_ratio_upper": str(tau_bound),
    "scalar_positive_upper": str(scalar_bound),
    "display": {"gamma": float(gamma), "craw": float(craw),
                "row_max": float(2*Uhi/qmin), "cap_ratio": float(cap_bound),
                "tau_ratio": float(tau_bound), "scalar_positive": float(scalar_bound)},
}
Path(__file__).with_name("verification.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({"status": "PASS", "checks": len(checks), "display": result["display"]}))
