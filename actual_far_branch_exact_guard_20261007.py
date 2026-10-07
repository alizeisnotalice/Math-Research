"""Exact constants and discrete coefficients for the actual-far branch audit.

No G_c evaluation, Monte Carlo, random seed, or actual FIRST sample occurs.
Binomial tails use the integer common denominator D**n, p0=1/D.
Large integers are represented by binary SHA256 and diagnostic log10 metadata.
"""
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt, log10
from pathlib import Path
import json


def rational(x):
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def integer_log10(x):
    if x == 0:
        return None
    shift = max(0, x.bit_length()-53)
    return log10(x >> shift)+shift*log10(2)


def integer_receipt(x):
    assert x >= 0
    data = x.to_bytes(max(1, (x.bit_length()+7)//8), "big")
    return dict(bit_length=x.bit_length(), unsigned_big_endian_sha256=sha256(data).hexdigest(),
                log10_approx=integer_log10(x))


def fraction_receipt(num, den):
    assert den > 0
    return dict(numerator=integer_receipt(num), denominator=integer_receipt(den),
                reduced=False, log10_approx=integer_log10(num)-integer_log10(den) if num else None,
                proof_arithmetic="integer numerator and stated common denominator; log10 is diagnostic only")


def tail_with_integer_denominator(n, D, K):
    """P(Bin(n,1/D)>K), exactly via its finite lower complement."""
    den = D**n
    term = (D-1)**n  # k=0 numerator
    low = term
    for k in range(K):
        term, rem = divmod(term*(n-k), (k+1)*(D-1))
        assert rem == 0  # binomial recurrence is integral
        low += term
    num = den-low
    assert 0 < num < den
    return num, den


def round_checks(n):
    s = isqrt(n)
    assert s*s == n and n >= 512
    D, K = 16*s, s-2
    p0, S0 = F(1, D), F(1, 16*s)
    eta0 = F(49, 65536)
    eta_n = eta0/s
    assert F(1, n) < S0 < 1
    assert K+1 > n*p0
    num, den = tail_with_integer_denominator(n, D, K)

    # tau <= 2**(-(K+1)) (1+p0)**n, using one common D**n.
    mgf_num = (D+1)**n
    mgf_den = den*2**(K+1)
    assert num*2**(K+1) <= mgf_num
    # Target many-active coefficient: n*tau <= 4*sqrt(n).
    assert n*num <= 4*s*den
    # Vertex budget, all comparisons use the original source normalization.
    M_checks = []
    for M in (3, s-1, K+2):
        assert M*eta_n <= eta0
        M_checks.append(dict(M=M, eta_n=rational(eta_n), M_eta_n=rational(M*eta_n),
                             eta0=rational(eta0), exact_pass=True))
    assert 4*s*eta_n == F(49, 16384)

    # Sample coefficient monotonicity in p without enormous Fraction gcds:
    # p0/4,p0/2,p0 share denominator E=4D.
    coefficient_checks = []
    E = 4*D
    common_den = E**n
    for k in (K+1, 2*K, n):
        assert k*D > n  # derivative sign for every 0<p<=p0
        values = [j**k*(E-j)**(n-k) for j in (1, 2, 4)]
        assert values[0] <= values[1] <= values[2]
        coefficient_checks.append(dict(k=k, positive_derivative_integer_margin=k*D-n,
                                       p_values=[rational(p0/4), rational(p0/2), rational(p0)],
                                       common_denominator=integer_receipt(common_den),
                                       numerator_receipts=[integer_receipt(x) for x in values], exact_pass=True))

    # Tiny-mark algebra from the stated analytic density bound G_(c,b)<=2/b.
    # Prior first-event total mass <=1; reference comparison is factor 2.
    tiny_checks = []
    for a_over_b in (F(1), F(3, 4), F(1, 2)):
        mark_bound = 4*a_over_b/n
        reference_source_mass_bound = 2*mark_bound
        assert reference_source_mass_bound <= F(8, n)
        tiny_checks.append(dict(a_over_b=rational(a_over_b), mark_mass_bound=rational(mark_bound),
                                reference_comparison_factor=2,
                                reference_source_mass_per_W=rational(reference_source_mass_bound),
                                common_upper_bound=rational(F(8, n)), exact_pass=True))
    Ch = F(16, 3)
    # Substitute only the proved bounds N_h<=n,D_n<=2sqrt(n).
    tiny_fee_per_W = F(8, n)*Ch*n*(2*s)
    assert tiny_fee_per_W == 16*Ch*s
    many_fee_per_W = F(2)*Ch*n*F(num, den)
    assert many_fee_per_W <= 8*Ch*s
    return dict(n=n, sqrt_n=s, p0=rational(p0), K=K,
                tail_tau=fraction_receipt(num, den),
                tail_method="D**n minus lower binomial numerator sum k=0..K; exact integral recurrence",
                mgf_bound=fraction_receipt(mgf_num, mgf_den),
                exact_tau_le_mgf=True, exact_n_tau_le_4_sqrt_n=True,
                n_tau=fraction_receipt(n*num, den),
                target_n_tau_bound=4*s,
                coefficient_monotonicity_samples=coefficient_checks,
                sigma_domain=dict(lower=rational(F(1, n)), upper=rational(S0),
                                  exact_lower_lt_upper=True, necessary_nonempty_condition_n_gt_256=n > 256,
                                  scope="Necessary parameter-range check only; no actual FIRST existence claim"),
                tiny_source_mass_checks=tiny_checks,
                tiny_space_fee_per_W=rational(tiny_fee_per_W),
                many_space_fee_per_W=fraction_receipt(2*Ch.numerator*n*num, Ch.denominator*den),
                many_fee_upper_bound_per_W=rational(8*Ch*s),
                vertex_checks=M_checks,
                exact_vertex_absorption_coefficient=rational(4*s*eta_n),
                pass_=True)


def main():
    rounds = [round_checks(n) for n in (1024, 4096, 16384)]
    result = dict(status="PASS_EXACT_ACTUAL_FAR_BRANCH_COMPONENTS", rounds=rounds,
                  random_seed="none: deterministic exact arithmetic; no random sampling",
                  scope="Analytic-constant substitutions and finite Bernoulli coefficient certificates only. No evaluation or simulation of G_c, no actual FIRST/GOOD/CP/GP instance, no spatial-envelope theorem established by numerics, and no certification of the remaining actual far traffic.",
                  large_integer_encoding="Unsigned big-endian byte SHA256 plus bit length and diagnostic log10; exact PASS comparisons run on full Python integers.",
                  output="actual_far_branch_exact_guard_20261007_results.json")
    path = Path(__file__).with_name(result['output'])
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps(dict(status=result['status'], n=[r['n'] for r in rounds],
                          coefficient_samples=sum(len(r['coefficient_monotonicity_samples']) for r in rounds),
                          output=str(path))))


if __name__ == "__main__":
    main()
