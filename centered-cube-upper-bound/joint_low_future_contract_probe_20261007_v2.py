#!/usr/bin/env python3
"""Three fixed true-kernel interval probes; no actual history/volume claim."""
from decimal import Decimal as D, Context, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import hashlib
import json

BASE = Path(__file__).resolve().parent
OUT = BASE / 'joint_low_future_contract_probe_20261007_v2.json'
LO = Context(prec=60, rounding=ROUND_FLOOR)
HI = Context(prec=60, rounding=ROUND_CEILING)


def dec(x, ctx):
    return ctx.divide(D(x.numerator), D(x.denominator))


def exp_bounds(l, u):
    # Decimal.exp is correctly rounded (half-even); adjacent 60-digit numbers
    # enclose the exact exponential regardless of the context rounding setting.
    return LO.next_minus(LO.exp(l)), HI.next_plus(HI.exp(u))


def phi_interval(t, steps):
    assert F(0) <= t <= F(1, 2)
    cs = [F(1, 2) + t, F(1, 2) - t]
    coeff = [(dec(c, LO), dec(c, HI)) for c in cs]
    left, right = D(0), D(0)
    for j in range(1, steps + 1):
        s = dec(F(j, steps), LO)  # steps is a power of two: exact decimal.
        fl, fu = D(0), D(0)
        for cl, cu in coeff:
            el, eu = exp_bounds(LO.divide(cu.copy_negate(), s), HI.divide(cl.copy_negate(), s))
            vl = max(D(0), LO.multiply(s, LO.subtract(D(1), eu)))
            vu = HI.multiply(s, HI.subtract(D(1), el))
            fl, fu = LO.add(fl, vl), HI.add(fu, vu)
        if j < steps:
            left = LO.add(left, fl)
        right = HI.add(right, fu)
    # s(1-exp(-c/s)) is increasing on [0,1], with value zero at zero.
    return F(LO.divide(left, D(steps))), F(HI.divide(right, D(steps)))


def log_interval(x, terms=32):
    assert x > 0
    shift = 0
    while x >= 2:
        x /= 2
        shift += 1
    while x < 1:
        x *= 2
        shift -= 1
    def core(z):
        low = 2 * sum((z ** (2*j+1)/(2*j+1) for j in range(terms)), F(0))
        high = low + 2*z**(2*terms+1)/((2*terms+1)*(1-z*z))
        return low, high
    low, high = core((x-1)/(x+1))
    l2, h2 = core(F(1, 3))
    return ((low+shift*l2, high+shift*h2) if shift >= 0
            else (low+shift*h2, high+shift*l2))


def ratio(x):
    return f'{x.numerator}/{x.denominator}'


def digest_fraction(x):
    def digest(k):
        raw = k.to_bytes(max(1, (k.bit_length()+7)//8), 'big')
        return dict(bits=k.bit_length(), sha256=hashlib.sha256(raw).hexdigest())
    return dict(numerator=digest(x.numerator), denominator=digest(x.denominator))


def main():
    assert not OUT.exists(), 'Frozen result refuses overwrite'
    rounds = []
    ln3lo, ln3hi = log_interval(F(3))
    for n, steps in [(512, 1024), (1024, 4096), (4096, 16384)]:
        k = (n + 1).bit_length()  # ceil(log2(n+2)); n+2 is not a power of two.
        alpha = isqrt(n)
        if alpha*alpha < n:
            alpha += 1
        radius, spacing = F(n+16, n), F(k, n)
        offsets = [F(1, 2), F(1, 2)-F(k, n+16)]
        phi = [phi_interval(t, steps) for t in offsets]
        avglo = sum(z[0] for z in phi)/2
        avghi = sum(z[1] for z in phi)/2
        assert F(3,8) < avglo <= avghi < F(3,4)
        Alo, Ahi = 1-avghi, 1-avglo
        rootlo, roothi = exp_bounds(dec(-ln3hi/n, LO), dec(-ln3lo/n, HI))
        slo, shi = (1-F(roothi))/Ahi, (1-F(rootlo))/Alo
        assert F(1,n) < slo <= shi < F(1,16)

        # Exact finite inequalities supporting the full-continuum analytic proof.
        below_radius_bound = (F(3,4)*radius)**n
        future_large_s_bound = F(7,8)**n
        assert below_radius_bound < F(1,3)
        assert future_large_s_bound < F(1,3)
        assert radius**n < 2**n
        beta_lower = F(16*n, n+16)
        assert beta_lower > 11  # v*=log(65536/49)<log(2^11)<11.
        assert radius-2*spacing < 1

        atom = F(1,2**n)
        eta_lower = F(49,65536*n)  # <= eta0/sqrt(n).
        greedy_threshold_lower = 2*eta_lower/(3*radius**n)
        assert atom < eta_lower
        assert atom < greedy_threshold_lower
        lnlo, lnhi = log_interval(F(n+2))
        llhi = log_interval(lnhi)[1]
        assert 2*llhi < k  # old h=2 a loglog(n+2)/n is smaller than d.
        lowZ = F(1,2**alpha) + F(9,10)**n
        epsilon_lower = 1/lnhi**4
        assert lowZ < epsilon_lower
        assert n*lowZ < F(1,16)
        rounds.append(dict(n=n, k=k, alpha=alpha, riemann_steps=steps,
            radius_over_a=ratio(radius), spacing_over_a=ratio(spacing),
            offsets_over_radius=[ratio(t) for t in offsets],
            original_phi_intervals=[[ratio(z[0]),ratio(z[1])] for z in phi],
            avg_phi_interval=[ratio(avglo),ratio(avghi)],
            sigma_interval=[ratio(slo),ratio(shi)],
            n_sigma_display_only=[float(n*slo),float(n*shi)],
            beta_certified_lower=ratio(beta_lower),
            single_atom_mass=digest_fraction(atom),
            beta_average_uniform_upper=digest_fraction(lowZ),
            epsilon_certified_lower=ratio(epsilon_lower),
            below_R0_future_relative_upper=digest_fraction(below_radius_bound),
            large_s_future_relative_upper=digest_fraction(future_large_s_bound),
            contracts=dict(hard_unique_winner=True, hardband=True,
                complete_continuous_future_cap=True, rightmost_FIRST_kernel_definition=True,
                pointwise_low_Z_all_original_sources=True, complete_output_low_average=True,
                microbox_nonconcentration=True, beta_above_vstar=True,
                hard_short_shell_mass='1-2^(-n)', old_loglog_microbox_greedy_not_triggered=True),
            limitations='One receiver point, finite positive atomic measure. No positive-volume, birth, CPGP, forest, or complete actual history certificate.'))
    result = dict(status='PASS_THREE_FIXED_ORIGINAL_KERNEL_JOINT_CONTRACT_PROBES_V2',
        rounds=rounds, arithmetic='Directed 60-digit Decimal with adjacent correctly rounded exp bounds; exact Fraction predicates; monotone Riemann enclosure. Display-only floats never used to pass.',
        preregistration='Only revised R0=a(1+16/n), d=ceil(log2(n+2))*a/n; no old winner=a or d=3a/n run.',
        supersedes='v1 decimal unary negation rounded exp input endpoints at the default 28-digit context. V1 interval results are invalid; frozen v1 files retained. V2 uses exact Decimal.copy_negate() at the same three parameters.',
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    with OUT.open('x') as f:
        json.dump(result,f,ensure_ascii=False,indent=2)
    print(json.dumps(dict(status=result['status'], rounds=[dict(n=r['n'],k=r['k'],n_sigma=r['n_sigma_display_only']) for r in rounds])))


if __name__ == '__main__':
    main()
