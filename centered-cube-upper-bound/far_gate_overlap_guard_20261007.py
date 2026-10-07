#!/usr/bin/env python3
"""Exact finite-row LP guards; not actual FIRST or spatial-cost experiments."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import hashlib
import json
import time

OUT = Path(__file__).resolve().parent
PREFIX = "far_gate_overlap_guard_20261007"
EPS = F(49, 4096)
SITES = [(0, 0), (1, 0), (0, 1), (1, 1)]


def fq(x):
    return str(x.numerator) + "/" + str(x.denominator)


def solve(rows, rhs):
    a = [[F(x) for x in row] + [F(b)] for row, b in zip(rows, rhs)]
    for j in range(4):
        pivot = next((i for i in range(j, 4) if a[i][j]), None)
        if pivot is None:
            return None
        a[j], a[pivot] = a[pivot], a[j]
        divisor = a[j][j]
        a[j] = [x / divisor for x in a[j]]
        for i in range(4):
            if i != j:
                coefficient = a[i][j]
                a[i] = [x - coefficient * y for x, y in zip(a[i], a[j])]
    return [a[i][4] for i in range(4)]


def vertex_minimum(rho, t):
    # Two equalities; choose two active inequalities from Ee lower bound and
    # four probability nonnegativity faces. Enumerates every bounded LP vertex.
    basic = [[1, 1, 1, 1], [0, 1, 0, 1]]
    active = [([0, 0, 1, 1], 1 - rho)]
    for j in range(4):
        active.append(([int(i == j) for i in range(4)], 0))
    feasible = []
    for first, second in combinations(active, 2):
        p = solve(basic + [first[0], second[0]], [1, t, first[1], second[1]])
        if p is not None and min(p) >= 0 and p[2] + p[3] >= 1 - rho:
            feasible.append(p)
    assert feasible
    best = min(feasible, key=lambda p: p[3])
    return best[3], best, len(feasible)


def check_row(rho, t):
    opt = max(F(0), t - rho)
    p = [rho - t + opt, t - opt, 1 - rho - opt, opt]
    assert min(p) >= 0 and sum(p) == 1
    assert p[1] + p[3] == t and p[2] + p[3] == 1 - rho
    dual = [-1, 1, 1] if t >= rho else [0, 0, 0]
    slacks = [F(g * e) - dual[0] - dual[1] * g - dual[2] * e for g, e in SITES]
    assert dual[2] >= 0 and min(slacks) >= 0
    dual_value = dual[0] + dual[1] * t + dual[2] * (1 - rho)
    assert dual_value == opt
    enum, best, vertex_count = vertex_minimum(rho, t)
    assert enum == opt
    return dict(rho=fq(rho), t=fq(t), optimum=fq(opt),
                primal_probabilities=[fq(x) for x in p],
                dual_total_gate_exit=[fq(F(x)) for x in dual],
                dual_site_slacks=[fq(x) for x in slacks],
                enumerated_minimum=fq(enum), enumerated_primal=[fq(x) for x in best],
                enumerated_feasible_active_sets=vertex_count,
                exact_gap="0/1")


def immutable_json(path, obj):
    with path.open("x") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)


def main():
    plan = [
        (1, 16, [F(3, 16), F(1, 4), F(3, 8)]),
        (2, 64, [F(3, 16), F(1, 5), F(1, 4), F(1, 3), F(3, 8), F(15, 32), F(49, 100)]),
        (3, 256, [F(3, 16), F(1, 5), F(7, 32), F(1, 4), F(5, 16), F(1, 3),
                  F(3, 8), F(7, 16), F(15, 32), F(49, 100), F(511, 1024)]),
    ]
    tags = [OUT / f"{PREFIX}_round{j}.json" for j, _, _ in plan]
    result_path = OUT / f"{PREFIX}_results.json"
    assert all(not p.exists() for p in tags + [result_path]), "No overwriting frozen outputs"
    round_receipts = []
    tic = time.time()
    for (j, den, ratios), path in zip(plan, tags):
        records = []
        for rho in ratios:
            # Uniform rational grid, kink, and absorption-threshold checks.
            nodes = set(F(k, den) for k in range(den + 1))
            nodes.update([rho, max(F(0), rho - F(1, den)), min(F(1), rho + F(1, den)),
                          EPS / 4, EPS / 8])
            records.extend(check_row(rho, t) for t in sorted(nodes))
        data = dict(status="PASS_EXACT_FINITE_GATE_EXIT_LP", round=j, grid_denominator=den,
                    rho_values=[fq(x) for x in ratios], rows=len(records), sites=SITES,
                    row_records=records, scope="Abstract finite probability rows and exact LP only; no actual FIRST or spatial measure.")
        immutable_json(path, data)
        receipt = dict(round=j, rows=len(records), file=path.name,
                       sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        round_receipts.append(receipt)
        print(json.dumps(receipt), flush=True)

    # Fixed-threshold and optimized ratio-dependent threshold certificates.
    thresholds = []
    for c in [F(0), F(1, 8), F(1, 4), F(1, 2), F(3, 4)]:
        for kappa in [F(3, 4), F(7, 8), F(1)]:
            u = F(4)
            rho = kappa / u
            cutoff = min(F(1), rho / (1 - c))
            rec = check_row(rho, cutoff)
            fee = u * cutoff
            assert fee == min(F(4), kappa / (1 - c))
            assert fee > EPS
            thresholds.append(dict(relative_coupling=fq(c), q_over_lambda=fq(kappa),
                                   u_over_lambda="4/1", cutoff=fq(cutoff),
                                   minimal_deleted_coefficient=fq(fee), extremizer=rec))
    additive = []
    for coefficient in [F(0), F(1, 2), F(1), F(2), F(8), F(1024)]:
        u, kappa = F(4), F(1)
        rho = kappa / u
        t = F(1) if coefficient < 1 else rho
        overlap = max(F(0), t - rho)
        required = u * (t - coefficient * overlap)
        optimum = 4 - 3 * coefficient if coefficient < 1 else F(1)
        assert required == optimum and required > EPS
        additive.append(dict(coupling_coefficient=fq(coefficient), t=fq(t), rho=fq(rho),
                             overlap=fq(overlap), sharp_absorption_coefficient=fq(optimum)))
    # Strict actual exit lower bound does not rescue a uniform bound: approach
    # Ee > 1-rho with exact positive slack while preserving J=0, t nearly rho.
    strict = []
    for den in [16, 64, 256]:
        rho = F(1, 4)
        slack = rho / den
        t = rho - slack
        probs = [F(0), t, 1 - rho + slack, F(0)]
        assert sum(probs) == 1 and min(probs) >= 0
        assert probs[2] > 1 - rho and probs[3] == 0
        strict.append(dict(slack=fq(slack), gate_mass=fq(t), exit_mean=fq(probs[2]),
                           overlap="0/1", traffic_over_lambda=fq(4 * t),
                           probabilities=[fq(p) for p in probs]))
    pointwise = []
    for n, den in [(512, 16), (4096, 64), (32768, 256)]:
        sigma = F(den + 1, den * n)
        bound = n * sigma / (4 + n * sigma)
        assert sigma > F(1, n) and bound > F(1, 5)
        pointwise.append(dict(n=n, sigma=fq(sigma), exit_lower_bound=fq(bound),
                              strict_margin_above_one_fifth=fq(bound - F(1, 5))))
    immutable_json(result_path, dict(status="PASS_THREE_ROUNDS_EXACT_GATE_OVERLAP_AND_THRESHOLD_GUARDS",
        total_rows=sum(x['rows'] for x in round_receipts), rounds=round_receipts,
        absorption_ceiling=fq(EPS), fixed_t_cutoff_at_ceiling=fq(EPS / 4),
        remaining_budget_after_existing_49_over_8192=fq(EPS / 2),
        fixed_t_cutoff_at_remaining_budget=fq(EPS / 8),
        optimal_adaptive_threshold_guards=thresholds, optimal_additive_guards=additive,
        strict_exit_approximation_rows=strict,
        pointwise_actual_exit_algebra_guards=pointwise,
        information_layer="LP thresholds use cap-average information only. Actual residual has the additional pointwise e>1/5, so LP zero-overlap extremizers are not feasible actual rows.",
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        elapsed_seconds=time.time() - tic,
        exact_arithmetic="Fraction throughout all comparisons; time is diagnostic floating only.",
        scope="Only general gate/exit row relaxation; no actual FIRST counterexample, new spatial fee, or theorem closure."))
    print("PASS three rounds", sum(x['rows'] for x in round_receipts), flush=True)


if __name__ == "__main__":
    main()
