#!/usr/bin/env python3
"""Deterministic kernel-component diagnostics and exact constant guards.

This does not sample G_c, sources, FIRST histories, or actual residual outputs.
The mathematical spectrum and mixture proof remains in the companion note.
"""
from fractions import Fraction as F
from math import comb, exp, isqrt, log, sqrt
from pathlib import Path
import json


def sqrt_upper(q, digits=6):
    scale = 10**digits
    k = isqrt(q.numerator * scale * scale // q.denominator)
    out = F(k + 1, scale)
    assert out * out >= q
    return out


def simpson(fun, left, right, panels):
    step = (right - left) / panels
    return step / 3 * (fun(left) + fun(right) + sum(
        (4 if j % 2 else 2) * fun(left + j * step)
        for j in range(1, panels)))


def component(u, panels):
    # Integrate in distance z=u(1/2-x) from the inner boundary.
    def fields(z):
        x = 0.5 - z / u
        left, right = exp(-z), exp(-u + z)
        b = 1 - (left + right) / 2
        bp = -u * (left - right) / 2
        return x, b, bp
    def integral(kind):
        def integrand(z):
            x, b, bp = fields(z)
            if kind == "J":
                return x*x*bp*bp/b
            if kind == "I":
                return (-1-x*bp/b)**2*b
            if kind == "mean":
                return -b-x*bp
            return b
        return 2/u * simpson(integrand, 0, u/2, panels)
    tail = 1-exp(-u)
    j_out = tail * (u/4+1+2/u)
    i_out = tail * (u/4+1/u)
    return {"J": integral("J")+j_out,
            "I": integral("I")+i_out,
            "mean": integral("mean")+tail/2,
            "mass": integral("mass")+tail/u,
            "J_in": integral("J")}


def scalar_diagnostic(u):
    coarse, fine = component(u, 8000), component(u, 16000)
    identities = abs(fine["mass"]-1) < 1e-9 and abs(fine["mean"]) < 1e-9
    identities &= abs(fine["I"]-(fine["J"]-1)) < 1e-9
    inequalities = fine["J_in"] <= u/4+1e-9 and fine["I"] <= u/2+2+1e-9
    stable = max(abs(fine[k]-coarse[k]) for k in fine) < 1e-8
    # Independently evaluate the convolution as a Laplace CDF difference.
    def cdf(y):
        return 1-exp(-u*y)/2 if y >= 0 else exp(u*y)/2
    errors = []
    for x in [0, 0.125, 0.499, 0.5, 0.501, 0.75, 1]:
        b = (1-(exp(-u*(0.5-x))+exp(-u*(0.5+x)))/2
             if x <= 0.5 else (exp(-u*(x-0.5))-exp(-u*(x+0.5)))/2)
        errors.append(abs(b-(cdf(x+0.5)-cdf(x-0.5))))
    assert identities and inequalities and stable and max(errors) < 1e-14
    return {"u": u, "fine": fine,
            "max_grid_difference": max(abs(fine[k]-coarse[k]) for k in fine),
            "convolution_cdf_max_error": max(errors),
            "status": "PASS", "certificate": False,
            "scope": "Floating-point scalar b_u diagnostic; no G_c or FIRST sample."}


def exact_round(n, rates):
    root = isqrt(n)
    ceilroot = root if root*root == n else root+1
    # ln(n+2)<log_2(n+2)<=(n+1).bit_length(), without floating logs.
    t_upper = F(12+(n+1).bit_length(), n)
    assert t_upper < F(1, 16)
    # e > sum_{j=0}^4 1/j! = 65/24 > 8/3.
    e_lower = F(65, 24)
    assert e_lower > F(8, 3)
    phi_half_lower = F(1, 2)-1/(3*e_lower)
    assert phi_half_lower > F(3, 8)
    b_lower = (F(3, 8)-t_upper)/(1-t_upper)
    assert b_lower >= F(5, 16)
    assert F(32, 5) < 7 and F(32, 5)**2 < 64
    assert 2+F(16, 2) == 10
    mask_checks = []
    for m in sorted({0, 1, ceilroot, n//2, n}):
        k = n-m
        # Triangle inequality proves the symbolic bound exactly; guard its
        # cross-term and retain rational upper radicals in the receipt.
        left_squared = m*m+10*k
        upper = F(m)+sqrt_upper(F(10*k))
        assert upper*upper >= left_squared
        mask_checks.append({"m": m, "k": k,
                            "TV_upper_rational": str(F(m)+sqrt_upper(F(left_squared))),
                            "triangle_status": "PASS"})
    parameter_checks = []
    for requested_L in [F(1), F(16), F(ceilroot)]:
        requested_r = requested_L/F(ceilroot)
        r0 = min(F(1, 2), requested_r)
        delta = (1-t_upper)*r0
        # Actual r=(1-sigma)/(1-t), with sigma>=1-delta.
        for t in [F(0), t_upper/2, t_upper]:
            assert delta/(1-t) <= r0
        variance = r0*(1-r0)+(1-r0)*10
        assert variance <= 11
        sn = sqrt_upper(F(n))
        s_regular = sqrt_upper(n*n*r0*r0+11*n)
        mixed = 8*sn*s_regular+8*n+8*n*r0*sn
        # log(b/a)<=log 2<1; this rational receipt bounds (16), not a
        # numerically measured maximal integral.
        J_safe = 1+8*r0*sn+sqrt_upper(F(10*n))+r0*mixed
        assert J_safe > 0
        parameter_checks.append({"requested_L_over_ceil_sqrt_n": str(requested_L),
                                 "r0_after_cap": str(r0), "delta_safe": str(delta),
                                 "J_safe_rational": str(J_safe),
                                 "actual_mapping_status": "PASS"})
    mode_checks = []
    for m in sorted({1, ceilroot, n//2}):
        # Exact binomial modal mass at r=m/n.  Verify the lower constant
        # without evaluating a radical or decimal probability.
        probability = F(comb(n,m)*m**m*(n-m)**(n-m), n**n)
        assert probability*probability*m >= F(9, 400)
        mode_checks.append({"m": m, "mode_lower_squared_status": "PASS"})
    return {"n": n, "T_upper_rational": str(t_upper),
            "phi_half_lower_rational": str(phi_half_lower),
            "b_inner_lower_rational": str(b_lower),
            "exact_constant_status": "PASS", "mask_checks": mask_checks,
            "parameter_checks": parameter_checks, "mode_checks": mode_checks,
            "scalar_diagnostics": [scalar_diagnostic(u) for u in rates]}


def main():
    result = {"name": "far_large_softness_fisher_guard_20261007",
              "seed": None, "randomness": "none",
              "scope": "Exact rational constant/domain/Bernoulli-mode guards plus deterministic scalar b_u diagnostics. No actual input/FIRST sample; no G_c spectral quadrature; no claim that the whole large-softness branch is paid.",
              "rounds": [exact_round(512, [1,2,8]),
                         exact_round(1024, [1.5,4,16]),
                         exact_round(4096, [3,32,64])],
              "status": "PASS"}
    output = Path(__file__).with_name("far_large_softness_fisher_guard_20261007_results.json")
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n")
    print("PASS: 3 rounds, 9 scalar diagnostics, exact domain/constant/TV/mapping/mode guards.")
    print(output)


if __name__ == "__main__":
    main()
