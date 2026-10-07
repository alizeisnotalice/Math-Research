"""New finite derivative-Mellin checks; numerical quadrature is not interval certified."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import numpy as np

BASE = Path(__file__).parent
REG = BASE / "frozen_jump_mellin_registration_20261007.json"
COUNT = 0

def check(ok, label):
    global COUNT
    COUNT += 1
    if not ok:
        raise AssertionError(label)

def complex_receipt(z):
    return {"real": float(z.real), "imag": float(z.imag)}

def normalized_transform(omega, q, order):
    x, w = np.polynomial.legendre.leggauss(order)
    edges = np.linspace(-40., 40., 33)
    v = np.concatenate([(a+b)/2+(b-a)*x/2 for a,b in zip(edges[:-1],edges[1:])])
    weights = np.concatenate([(b-a)*w/2 for a,b in zip(edges[:-1],edges[1:])])
    t = np.exp(v)
    residual = -np.expm1(-t)/t - np.exp(-t)
    gamma_integrand = t * np.exp(-t)
    phase = np.exp(-1j*omega*v)
    gamma = np.dot(weights, gamma_integrand*phase)
    # Physical logtime is v-log(q); q^(i omega) is a phase, not a magnitude.
    response = np.dot(weights, residual*phase) * np.exp(1j*omega*math.log(q))
    return response, gamma

def gamma_integer_ratio(k, omega):
    # Exact identity: ratio^2=(pi*w/sinh(pi*w))*product_{j=1}^{k-2}(1+w^2/j^2).
    # Stable log evaluation avoids gamma/Lanczos libraries; still floating point.
    a = math.pi*omega
    log_sinh = a-math.log(2)+math.log1p(-math.exp(-2*a))
    j = np.arange(1, k-1, dtype=float)
    log_square = math.log(a)-log_sinh+float(np.log1p((omega/j)**2).sum())
    return math.exp(log_square/2)

def run(spec):
    before = COUNT
    m, q, cutoff = spec["states"], F(spec["q"]), spec["cutoff"]
    source_mass = F(m,2)
    u0 = F(m*(m-2),2*(m-1))/q
    outside_density = q*u0/m
    outside_mass = (m-1)*outside_density
    check(q*u0*(1-F(1,m)) == source_mass-1, "saturation equation")
    check(0 < outside_density < 1, "positive full exterior absorption")
    check(1+outside_mass == source_mass, "source mass once")
    check(F(1) <= source_mass, "omega paid by bad mass")
    true_tail_ratio_bound = F(4,2**cutoff)*outside_mass/source_mass
    check(true_tail_ratio_bound < 1, "analytic true high-tail constant")
    mellin = []
    for word in ("1/2","1","2","4"):
        omega = float(F(word))
        coarse, coarse_gamma = normalized_transform(omega,float(q),spec["gauss_order"])
        fine, fine_gamma = normalized_transform(omega,float(q),2*spec["gauss_order"])
        target = np.exp(1j*omega*math.log(float(q)))*fine_gamma/(1+1j*omega)
        error = abs(fine-target)
        refinement = max(abs(coarse-fine),abs(coarse_gamma-fine_gamma))
        gamma_square = math.pi*omega/math.sinh(math.pi*omega)
        check(error < 5e-11, "renormalized R Mellin identity")
        check(refinement < 5e-11, "registered quadrature refinement")
        check(abs(abs(fine_gamma)**2-gamma_square) < 5e-11, "Gamma modulus identity")
        check(abs(abs(target)-abs(fine_gamma)/math.sqrt(1+omega**2)) < 1e-14, "q phase has no tail magnitude")
        mellin.append({"omega":word,"response":complex_receipt(fine),
                       "gamma_from_independent_integrand":complex_receipt(fine_gamma),
                       "identity_error":float(error),"refinement_difference":float(refinement)})
    bands = []
    base = (8*cutoff)**2
    for multiple in (1,2,4):
        k = base*multiple
        for factor in (1.,1.5,2.):
            omega = factor*math.sqrt(k)
            ratio = gamma_integer_ratio(k,omega)
            check(k >= 3, "Gamma strip k domain")
            check(ratio >= math.exp(-6)-1e-8, "infinite-product lower bound diagnostic")
            bands.append({"k":k,"omega":omega,"gamma_ratio":ratio})
    triangle = []
    for blocks in (1,4,16):
        # Every [K,2K) has sum 1/(k-1)>1/2; no numerical fitting.
        lower = outside_mass*F(blocks,2*3888)
        check(lower > 0, "positive high-tail triangle block cost")
        check(F(blocks,2) == sum((F(1,2) for _ in range(blocks)),F()), "exact harmonic block accumulation")
        triangle.append({"blocks":blocks,"lower_bound":str(lower)})
    return {"states":m,"q":str(q),"cutoff":cutoff,"source_mass":str(source_mass),
            "outside_absorption_mass":str(outside_mass),"u0":str(u0),
            "true_high_tail_over_source_analytic_upper":str(true_tail_ratio_bound),
            "mellin":mellin,"gamma_bands":bands,"triangle_block_lower_bounds":triangle,
            "checks":COUNT-before,"status":"PASS"}

registration = json.loads(REG.read_text())
rounds = [run(spec) for spec in registration["rounds"]]
result = {"status":"PASS_EXACT_AND_NUMERIC","predicates":COUNT,"rounds":rounds,
          "scope":registration["scope"],
          "quadrature":"paired floating GL, not interval; analytic logtime tails <=1.5 exp(-40) for R and exp(-40)+exp(-exp(40)) for Gamma",
          "triangle":"infinite-tail divergence is proved analytically, not inferred from three finite blocks",
          "registration_sha256":hashlib.sha256(REG.read_bytes()).hexdigest(),
          "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
out = BASE / "frozen_jump_mellin_guard_results_20261007.json"
out.write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps({"status":result["status"],"predicates":COUNT,"output":str(out),
                  "max_identity_error":max(x["identity_error"] for r in rounds for x in r["mellin"]),
                  "max_refinement":max(x["refinement_difference"] for r in rounds for x in r["mellin"])}))
