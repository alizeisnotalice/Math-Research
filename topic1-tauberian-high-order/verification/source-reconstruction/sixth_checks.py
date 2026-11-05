"""Small exact algebra checks for the R063 independent reconstruction.

These identities are supporting evidence only.  They do not check any covering,
measurability, witness, compactness, incidence, or supremum argument.
"""
from pathlib import Path
import json
import sympy as s

v, a, V, M, S, Av, c, b, A0 = s.symbols(
    "v a V M S Av c b A0", positive=True)
f = s.symbols("f0:7")
tau = b*v/c
checks = {}
lhs = sum(tau**(6-h)*(c*f[h]-b*v*f[h-1]) for h in range(1, 7))
checks["six_layer_telescope"] = s.expand(
    lhs-(c*f[6]-b*v*tau**5*f[0])) == 0
checks["six_layer_initial_endpoint"] = s.simplify(
    (b*v*tau**5*(A0*v*M)/a)
    - (A0*b**6/c**5)*v**7*M/a) == 0

q = s.symbols("q", real=True)
zeta = q/(2*q-1)
checks["holder_exponent_factorization"] = s.factor(
    q+6*zeta-6) == (q-2)*(2*q-3)/(2*q-1)
checks["holder_small_set_power"] = s.simplify(
    1-zeta-(q-1)/(2*q-1)) == 0
checks["higher_support_minimum_gap_factor"] = s.simplify(
    q+6*zeta-6-s.Rational(3,5)
    -(q-3)*(10*q-11)/(5*(2*q-1))) == 0

L1, L2, z = s.symbols("L1 L2 z", real=True)
poly = (L1+z)*(L2+z)*(L1+L2+2*z)
checks["positive_branch_convexity"] = s.expand(
    s.diff(poly, z, 2)-6*(L1+L2+2*z)) == 0

n, beta, u = s.symbols("n beta u", positive=True)
k = n**2/(27*beta**2)
t = 3*beta*u/n
checks["scalar_maximum_factorization"] = s.factor(
    k-u**2*(1-2*beta*u/n)-k*(t-1)**2*(2*t+1)) == 0

MI, SI = s.symbols("MI SI", real=True)
checks["interior_identity"] = s.expand(
    (1-a)*SI+MI-Av*MI/a
    - ((1-a)*(SI-MI/a)-(Av-1)*MI/a)) == 0

rho, N = s.symbols("rho N", positive=True)
delta = rho/(N*(1+2*rho))
radius = 1/(N*(1+2*rho))
checks["periodic_neighbor_touch"] = s.simplify(radius+2*delta-1/N) == 0
checks["periodic_boundary_slack"] = s.simplify(
    1/N-radius-delta-rho/(N*(1+2*rho))) == 0

result = {
    "sympy_version": s.__version__,
    "checks": checks,
    "all_passed": all(checks.values()),
    "scope": "Exact algebra only; all-dimensional geometry is proved in rebuild.md."
}
assert result["all_passed"], result
Path(__file__).with_name("algebra_results.json").write_text(
    json.dumps(result, ensure_ascii=False, indent=2)+"\n")
print(json.dumps(result, ensure_ascii=False, indent=2))
