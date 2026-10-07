#!/usr/bin/env python3
"""New original-symbol guard. No reset kernels, no old data rerun.

Floating GL pairs diagnose analytic inequalities and contour signs. They are
not interval certificates. Infinite tails have separately stated bounds.
"""
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
PREFIX = "frozen_jump_allvisited_hormander"
REG = HERE / (PREFIX + "_registration_20261007.json")
OUT = HERE / (PREFIX + "_guard_results_20261007.json")
registration = json.loads(REG.read_text())
count = 0

def check(test, name):
    global count
    count += 1
    if not bool(test):
        raise AssertionError(name)

def B(v):
    v = np.asarray(v, dtype=float)
    # Stable near zero; series through v^4. Grid never reaches below 2^-24.
    small = v < 1e-3
    safe = np.where(small, 1.0, v)
    exact = safe / np.log1p(safe) - 1.0
    series = v/2 - v*v/12 + v**3/24 - 19*v**4/720
    return np.where(small, series, exact)

def elasticity(v):
    if v < 1e-3:
        b = float(B(v))
        derivative = .5-v/6+v*v/8-19*v**3/180
        return v*derivative/b
    log = math.log1p(v)
    derivative = (log-v/(1+v))/(log*log)
    return v*derivative/float(B(v))

def invB(y):
    lo, hi = 0.0, max(1.0, 4*y)
    while float(B(hi)) < y:
        hi *= 2
    for _ in range(100):
        mid = (lo+hi)/2
        if float(B(mid)) < y:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2

def GL(fun, left, right, order, panels=1):
    z, w = np.polynomial.legendre.leggauss(order)
    result = 0.0
    for j in range(panels):
        a = left+(right-left)*j/panels
        b = left+(right-left)*(j+1)/panels
        x = (a+b)/2+(b-a)*z/2
        result += (b-a)/2*np.dot(w, fun(x))
    return result

def peak(c, a, order):
    t = c*a
    bt, bc = invB(1/t), invB(1/c)
    ratio = math.sqrt(bc/bt)
    upper = 65536*ratio
    def integrand(u):
        bv = B(bt*u*u)
        g = 1/(1+c*bv)
        return np.exp(-a*(1-g))*(-np.expm1(-a*g))
    partial = GL(integrand, 0, 1, order)
    partial += GL(lambda w: np.exp(w)*integrand(np.exp(w)),
                  0, math.log(upper), order, 16)
    partial /= math.pi
    tail = 3*a*math.exp(-a/2)*(bc/bt)**(2/3)*upper**(-1/3)/math.pi
    return partial, tail, bt, bc

def direction(n, kind):
    base = (1+np.arange(n)/(2*n))/math.sqrt(n)
    if kind == "balanced":
        return base
    if kind == "one_large_coordinate":
        base[0] = math.sqrt(n)
        return base
    return 2.0**(np.arange(n)%5-2)/math.sqrt(n)

def contour(n, c, zeta, omega, order, shifted):
    beta = 1/256
    psi = B(zeta*zeta)/(1+c*B(zeta*zeta))
    sumpsi = float(psi.sum())
    low = math.log(c)-40/n
    upper_t = (40+n*math.log(2))/(math.cos(beta)*sumpsi)
    high = math.log(upper_t)
    angle = -beta if shifted else 0.0
    def fun(v):
        t = np.exp(v+1j*angle)
        # Difference represented without catastrophic cancellation.
        factors = np.exp(-t[:,None]*psi[None,:]) * (
            -np.expm1(-t[:,None]*(1/c-psi)[None,:]))
        whole = np.prod(factors, axis=1)
        return np.exp(-1j*omega*v+omega*angle)*whole
    val = GL(fun, low, high, order, 32)
    val *= 1j*omega/(1+1j*omega)
    small_tail = math.exp(-40)/n
    big_arg = math.cos(beta)*upper_t*sumpsi
    big_tail = math.exp(n*math.log(2)-big_arg)/big_arg
    return complex(val), small_tail+big_tail, [low,high]

rounds = []
for rr in registration["rounds"]:
    n = rr["n"]
    q0, q1 = rr["gl_orders"]
    symbol = []
    for exponent in registration["symbol_lambda_log2_grid"]:
        v = 2.0**exponent
        e = elasticity(v)
        check(2/3-1e-10 <= e <= 1+1e-10, "original B elasticity")
        for scale in (2,16,256):
            ratio = float(B(scale*v)/B(v))
            check(scale**(2/3)*(1-1e-9) <= ratio <= scale*(1+1e-9),
                  "original B scaling")
        symbol.append({"lambda":v,"elasticity":e})
    density = []
    for c in registration["c"]:
        for a in registration["density_t_over_c"]:
            coarse, tail, bt, bc = peak(c,a,q0)
            fine, tail2, _, _ = peak(c,a,q1)
            check(abs(tail-tail2) < 1e-15, "fixed analytic tail")
            check(0 <= fine <= fine+tail <= 5, "density peak analytic bound")
            check(abs(fine-coarse) < 1e-7, "peak quadrature refinement")
            check(math.sqrt(bc/bt) <= a**.75*(1+1e-12), "inverse scaling")
            # Arbitrary h tested relative to the actual inverse clock.
            for relative in (1/16,.25,1,4,16):
                r = relative/math.sqrt(bt)
                A = float(B(r**-2))/(1+c*float(B(r**-2)))
                s = c*a*A
                if s >= 1:
                    check(r*math.sqrt(bt) <= s**-.5*(1+1e-10),
                          "inverse concentration direction")
            density.append({"c":c,"t_over_c":a,"normalized_peak_coarse":coarse,
                "normalized_peak_fine":fine,"analytic_positive_tail_upper":tail,
                "quadrature_pair_difference":abs(fine-coarse),"beta_t":bt,"beta_c":bc})
    contours = []
    for c in registration["c"]:
        for kind in registration["fourier_directions"]:
            zeta = direction(n,kind)
            for omega in registration["mellin_omega"]:
                real0, tails, limits = contour(n,c,zeta,omega,q0,False)
                real1, _, _ = contour(n,c,zeta,omega,q1,False)
                shift0, _, _ = contour(n,c,zeta,omega,q0,True)
                shift1, _, _ = contour(n,c,zeta,omega,q1,True)
                check(abs(real1-real0) < 1e-7, "real contour refinement")
                check(abs(shift1-shift0) < 1e-7, "shift contour refinement")
                check(abs(real1-shift1) < 1e-7+2*tails, "full contour shift sign")
                contours.append({"c":c,"direction":kind,"omega":omega,
                    "real":[real1.real,real1.imag],"shifted":[shift1.real,shift1.imag],
                    "contour_discrepancy":abs(real1-shift1),
                    "real_GL_pair_difference":abs(real1-real0),
                    "shift_GL_pair_difference":abs(shift1-shift0),
                    "analytic_two_endpoint_tail_upper":tails,"logtime_limits":limits})
    # Exact coarse exponent margins: valid for all n>=1, checked at saved n.
    hormander_exponent = Fraction(12+3*n-24*n*n)
    check(hormander_exponent <= -9, "H cutoff exponent")
    check(Fraction(1,2**9) >= Fraction(1,2**(-hormander_exponent)), "H < 1/512")
    # n^(3/2)>=n gives conservative second L2 exponent.
    l2_exponent = Fraction(2)+Fraction(7*n,2)-Fraction(3072*n,5)
    check(l2_exponent < -10, "L2 cutoff exponent")
    check(Fraction(8,2**(1024*n*n-n)) < Fraction(1,1024), "first L2 cutoff")
    check(Fraction(1,2**10)+Fraction(1,2**10) == Fraction(1,512), "L2 sum budget")
    theta = Fraction(1,128)
    check(24*theta < Fraction(1,2), "Taylor sector safe")
    check(Fraction(1,256) == theta/2, "contour beta")
    check(Fraction(4)*(Fraction(1,512)+Fraction(1,512)) == Fraction(1,64), "CZ constant")
    rounds.append({"n":n,"gl_orders":rr["gl_orders"],"symbol":symbol,
        "density":density,"contours":contours,
        "cutoff":1024*n*n,"exact_H_exponent_upper":str(hormander_exponent),
        "exact_L2_exponent_upper":str(l2_exponent)})

result = {"status":"PASS_EXACT_AND_NUMERIC","predicates":count,
    "registration_sha256":hashlib.sha256(REG.read_bytes()).hexdigest(),
    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "scope":"Original B/psi_c symbol, real-density bounds and complete full-product Mellin contour; floating quadrature is not an interval proof. Uniform c, spatial Hörmander and fixed maximal conclusions are analytic.",
    "rounds":rounds}
OUT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
print(result["status"],count,"predicates")
