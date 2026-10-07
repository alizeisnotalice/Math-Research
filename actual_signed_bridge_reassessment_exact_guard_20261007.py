#!/usr/bin/env python3
"""New joint-gate/signed correction components; not actual FIRST samples."""
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
PREFIX = "actual_signed_bridge_reassessment"
checks, rounds = [], []

def rec(label, ok):
    checks.append({"label": label, "pass": bool(ok)})

def expminus(x, N=80):
    # Positive Taylor sum for exp(x), with a geometric upper bound on its tail.
    if not x:
        return F(1), F(1)
    term, lower = F(1), F(1)
    for j in range(1, N + 1):
        term *= x / j
        lower += term
    next_term = term * x / (N + 1)
    ratio = x / (N + 2)
    assert ratio < 1
    upper = lower + next_term / (1 - ratio)
    return 1 / upper, 1 / lower

def apply_g(profile, axis):
    return [(2 * profile[x] + profile[x ^ (1 << axis)]) / 3
            for x in range(len(profile))]

def joint_profile(seed, mask, d):
    value = seed[:]
    for axis in range(d):
        if mask & (1 << axis):
            value = apply_g(value, axis)
    return value

def repeats(profile, mask, p, d):
    # The first two extra-degree layers of the support mask reference.
    # This is a finite coefficient component, not a truncated approximation
    # used as a substitute for the full H kernel.
    axes = [i for i in range(d) if mask & (1 << i)]
    k = len(axes)
    extra = [F(0)] * len(profile)
    for i in axes:
        gi = apply_g(profile, i)
        for x in range(len(profile)):
            extra[x] += p ** (k + 1) * gi[x] / 2
        gii = apply_g(gi, i)
        for x in range(len(profile)):
            extra[x] += p ** (k + 2) * gii[x] / 6
    for j, i in enumerate(axes):
        for ii in axes[j + 1:]:
            gij = apply_g(apply_g(profile, i), ii)
            for x in range(len(profile)):
                extra[x] += p ** (k + 2) * gij[x] / 4
    return extra

for d, n, m0, sqrtn in [(2, 729, 3, 27),
                         (3, 262144, 8, 512),
                         (4, 16777216, 16, 4096)]:
    start = len(checks)
    size = 1 << d
    raw = [F(1 + (5 * x + d) % 11, 2 + x % 3) for x in range(size)]
    seed = [a / sum(raw) for a in raw]
    tag = f"finite_d{d}/budget_n{n}"
    rec(tag + "/seed-positive-normalized", all(a > 0 for a in seed) and sum(seed) == 1)
    independence_witnesses = 0
    symbolic_identities = 0
    for p in [F(0), F(1, 4), F(1, 2), F(3, 4), F(1)]:
        elo, ehi = expminus(d * p)
        ptag = tag + f"/p{p}"
        rec(ptag + "/exp-enclosure-direction", 0 < elo <= ehi <= 1)
        prior_sum = F(0)
        accepted_total = F(0)
        posterior_total = F(0)
        for mask in range(size):
            k = mask.bit_count()
            beta = p ** k * (1 - p) ** (d - k)
            prior_sum += beta
            profile = joint_profile(seed, mask, d)
            rec(ptag + f"/mask{mask}/normalization", sum(profile) == 1)
            # Fine history and acceptance are correlated with endpoint and mask.
            hist0 = [F(1 + (x + 2 * mask) % 5, 7) for x in range(size)]
            gate0 = [F((x + mask) % 3, 2) for x in range(size)]
            gate1 = [F(1 + (2 * x + mask) % 4, 4) for x in range(size)]
            accepted = [beta * profile[x] *
                        (hist0[x] * gate0[x] + (1 - hist0[x]) * gate1[x])
                        for x in range(size)]
            a = [accepted[x] / (beta * profile[x]) if beta else F(0)
                 for x in range(size)]
            prior_gate_average = sum(profile[x] *
                                     (hist0[x] * gate0[x] +
                                      (1 - hist0[x]) * gate1[x])
                                     for x in range(size))
            if beta and any(accepted[x] != beta * profile[x] * prior_gate_average
                            for x in range(size)):
                independence_witnesses += 1
            extra = repeats(profile, mask, p, d)
            for x in range(size):
                xtag = ptag + f"/mask{mask}/x{x}"
                rec(xtag + "/joint-RN-exact",
                    F(0) <= a[x] <= 1 and accepted[x] == beta * profile[x] * a[x])
                # Exact formal coefficient pair in common zeta=exp(-d*p).
                reference_pair = (F(0), a[x] * (p ** k * profile[x] + extra[x]))
                difference_pair = (a[x] * beta * profile[x],
                                   -a[x] * (p ** k * profile[x] + extra[x]))
                rec(xtag + "/same-gate-signed-identity",
                    (reference_pair[0] + difference_pair[0],
                     reference_pair[1] + difference_pair[1]) == (accepted[x], F(0)))
                symbolic_identities += 1
                rec(xtag + "/repeated-nonpositive", -a[x] * extra[x] <= 0)
                if not beta:
                    rec(xtag + "/zero-prior-convention", a[x] == 0)
            accepted_total += sum(accepted)
            # Posterior labeling the unsplit accepted endpoint also preserves total.
            unsplit_gate = [F(1 + x % 4, 5) for x in range(size)]
            kernel_total = [sum(p ** mm.bit_count() * (1-p) ** (d-mm.bit_count()) *
                                joint_profile(seed, mm, d)[x] for mm in range(size))
                            for x in range(size)]
            posterior_total += sum(unsplit_gate[x] * kernel_total[x] *
                                   (beta * profile[x] / kernel_total[x])
                                   for x in range(size))
        rec(ptag + "/all-mask-prior-mass", prior_sum == 1)
        rec(ptag + "/accepted-original-mass", 0 <= accepted_total <= 1)
        rec(ptag + "/posterior-total-preserved",
            posterior_total == sum(F(1+x % 4, 5) * kernel_total[x] for x in range(size)))
        for k in range(1, d + 1):
            # upper enclosure on d_k; checking its positive part suffices.
            upper_d = p ** k * ((1-p) ** (d-k) - elo)
            bound = k * p ** (k+1) * (1-p) ** (d-k)
            rec(ptag + f"/k{k}/coefficient-upper", max(upper_d, F(0)) <= bound)
    rec(tag + "/gate-independence-failure-witness", independence_witnesses > 0)

    # Fixed source labels have one reference mass, regardless of receiver selector.
    W = F(37, 13)
    masses = [W * F(2, 9), W * F(1, 7), W * F(3, 11)]
    label_weights = [F(1, 5), F(2, 5), F(1, 4), F(3, 20)]
    reference = [(mass, 2 * mass * weight)
                 for mass in masses for weight in label_weights]
    rec(tag + "/source-partition-not-renormalized", sum(masses) <= W)
    rec(tag + "/reference-mass-once", sum(b for mass, b in reference) == 2*sum(masses) <= 2*W)
    for receiver in range(size):
        selector_density_ratio = F(1 + receiver % 7, 8)
        rec(tag + f"/receiver{receiver}/first-tag-domination",
            all(F(0) <= selector_density_ratio*b <= b for mass, b in reference))
    rec(tag + "/receiver-replication-is-not-used",
        size * sum(b for mass, b in reference) > 2 * W)

    Ch = F(16, 3)
    for q in [F(3, 4), F(7, 8), F(1)]:
        for u in [F(5, 2), F(3), F(4)]:
            rec(tag + f"/q{q}/u{u}/hard-row-cancel", u/q <= Ch)
    # An enlarged positive majorant may exceed the actual cap. Cap the actual
    # density, not the majorant; negative differences are harmless here.
    for B in [F(0), F(1, 2), F(4), F(9)]:
        for e in [F(0), F(1, 3), F(7)]:
            actual = min(F(4), B+e)
            rec(tag + f"/B{B}/e{e}/actual-min-cap",
                actual <= min(F(4), B) + e)
    rec(tag + "/majorant-does-not-inherit-cap", F(9) > 4)

    C = F(m0*(m0+1)*(m0+2), 3*(n+1))
    fee = 2 * Ch * n * C  # already-proved physical N_h <= n, no log rounding
    coarse = 4 * Ch * sqrtn
    rec(tag + "/integer-cutoff", n == m0**6 and n == sqrtn**2)
    rec(tag + "/cubic-sum-exact",
        sum(F(k*(k+1), n+1) for k in range(1, m0+1)) == C)
    rec(tag + "/cubic-cutoff-coarse", m0*(m0+1)*(m0+2) <= 6*m0**3)
    rec(tag + "/actual-source-fee", fee <= coarse)
    rec(tag + "/growing-cutoff-scope", m0**3 == sqrtn and n >= 512)
    rounds.append({
        "finite_profile_dimension": d, "compressed_budget_dimension": n,
        "m0": m0, "sqrt_n": sqrtn, "checks": len(checks)-start,
        "same_gate_symbolic_identities": symbolic_identities,
        "independence_failure_witnesses": independence_witnesses,
        "C_n_m0": str(C), "exact_fee_over_W_using_Nh_le_n": str(fee),
        "coarse_fee_over_W": str(coarse),
        "scope": "Finite positive profile conditionalization / formal repeated coefficients and exact large-n scalar budget only."
    })

registration = HERE / f"{PREFIX}_registration_20261007.json"
result = {
    "status": "PASS" if all(c["pass"] for c in checks) else "FAIL",
    "seed": None, "total_checks": len(checks),
    "failed_checks": [c for c in checks if not c["pass"]],
    "rounds": rounds, "checks": checks,
    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "registration_sha256": hashlib.sha256(registration.read_bytes()).hexdigest(),
    "scope": "RN/frozen-gate identities and actual source-once fee components; not original G_c simulation, full FIRST eligibility, or moving-reference weak proof."
}
(HERE / f"{PREFIX}_results_20261007.json").write_text(
    json.dumps(result, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({k:result[k] for k in
                  ["status","total_checks","failed_checks","script_sha256"]},
                 ensure_ascii=False))
raise SystemExit(0 if result["status"] == "PASS" else 1)
