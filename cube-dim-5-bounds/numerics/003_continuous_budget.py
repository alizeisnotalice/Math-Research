#!/usr/bin/env python3
"""Finite Chebyshev/tube budget scan; this is a planning computation only."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp


mp.mp.dps = 80
P = mp.mpf(35509) / (mp.mpf(764411904) * mp.pi ** 5)
CF = mp.mpf(14455) / (mp.mpf(31104) * mp.pi ** 5)
SECONDS_PER_ARB_NODE = mp.mpf('2.6')
G_PROXY = mp.mpf('0.000146450542661776629693871130288973764')
B_PROXY = mp.mpf('0.00003791184534255565')


def tube_error(epsilon: mp.mpf, subdivisions: int, degree: int, rho: mp.mpf) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    h = mp.mpf(1) / (2 * subdivisions)
    height = h * (rho - 1 / rho) / 2
    assert height < 6 * epsilon
    m_bound = P / (4 * (epsilon - height / 6) ** 4)
    e1 = 4 * m_bound * rho ** (-degree) / (rho - 1)
    lam = 2 * degree + 1
    return e1 * (1 + lam + lam ** 2 + lam ** 3), height, m_bound


def short_l54_error(epsilon: mp.mpf) -> mp.mpf:
    return mp.power(16, -mp.mpf(4) / 5) * mp.mpf(5) / 6 * mp.power(CF, mp.mpf(1) / 5) * mp.power(16 / mp.pi, mp.mpf(4) / 5) * mp.power(epsilon, mp.mpf(6) / 5)


def fmt(value: mp.mpf) -> str:
    return mp.nstr(value, 24)


def main() -> None:
    epsilons = [mp.mpf(1) / 2 ** exponent for exponent in (5, 6, 7, 8, 9, 10)]
    subdivisions = (4, 8, 16, 32)
    degrees = (8, 12, 16, 24, 32)
    rows = []
    # A fixed finite rho grid per requested parameter tuple.  The final entry is
    # strictly inside the tube boundary; log spacing resolves rho close to one.
    for epsilon in epsilons:
        for m in subdivisions:
            h = mp.mpf(1) / (2 * m)
            rho_max = (12 * epsilon / h + mp.sqrt((12 * epsilon / h) ** 2 + 4)) / 2
            rho_grid = [mp.e ** (mp.log(rho_max) * k / 256) for k in range(1, 256)]
            for degree in degrees:
                candidates = [(tube_error(epsilon, m, degree, rho)[0], rho, *tube_error(epsilon, m, degree, rho)[1:]) for rho in rho_grid]
                error, rho, height, m_bound = min(candidates, key=lambda item: item[0])
                nodes = m ** 4 * (degree + 1) ** 4
                rows.append({"epsilon": fmt(epsilon), "m": m, "d": degree, "rho": fmt(rho), "tube_height": fmt(height), "E4": fmt(error), "E4_over_g_proxy": fmt(error / G_PROXY), "M": fmt(m_bound), "short_L5_4_error": fmt(short_l54_error(epsilon)), "nodes": nodes, "serial_arb_seconds_at_2_6_per_node": fmt(SECONDS_PER_ARB_NODE * nodes), "serial_arb_days_at_2_6_per_node": fmt(SECONDS_PER_ARB_NODE * nodes / 86400)})
    # Pareto set: no row simultaneously improves E4 and serial node cost.
    pareto = []
    for row in rows:
        err = mp.mpf(row["E4"]); work = mp.mpf(row["serial_arb_seconds_at_2_6_per_node"])
        if not any(mp.mpf(other["E4"]) <= err and mp.mpf(other["serial_arb_seconds_at_2_6_per_node"]) <= work and (mp.mpf(other["E4"]) < err or mp.mpf(other["serial_arb_seconds_at_2_6_per_node"]) < work) for other in rows):
            pareto.append(row)
    feasible = [row for row in rows if mp.mpf(row["E4_over_g_proxy"]) <= mp.mpf('0.01')]
    short_proxy = []
    for epsilon in epsilons:
        short_error = short_l54_error(epsilon)
        i_proxy = (B_PROXY ** (mp.mpf(4) / 5) + short_error) ** (mp.mpf(5) / 4)
        c_proxy = 5 ** (mp.mpf(5) / 4) * i_proxy / G_PROXY ** (mp.mpf(5) / 4)
        short_proxy.append({"epsilon": fmt(epsilon), "short_L5_4_error": fmt(short_error), "Minkowski_I_proxy": fmt(i_proxy), "cost_proxy": fmt(c_proxy)})
    output = {"status": "planning_budget_only", "peak_constant": "35509/(764411904*pi^5)", "conditional_Cf": "14455/(31104*pi^5)", "short_error_formula": "16^(-4/5)*(5/6)*Cf^(1/5)*(16/pi)^(4/5)*epsilon^(6/5)", "assumptions": {"seconds_per_arb_node": "2.6", "g_proxy_from_002G_generic_sample": fmt(G_PROXY), "B_proxy_from_002D_deterministic_gauss_sum": fmt(B_PROXY), "warning": "Neither proxy is a certified global g or I_H bound; Cf is still under review."}, "short_time_minkowski_proxy": short_proxy, "rows": rows, "pareto_by_E4_and_serial_cost": sorted(pareto, key=lambda row: mp.mpf(row["serial_arb_seconds_at_2_6_per_node"])), "rows_with_E4_at_most_1pct_g_proxy": sorted(feasible, key=lambda row: mp.mpf(row["serial_arb_seconds_at_2_6_per_node"])), "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path("003_continuous_budget.json").write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"row_count": len(rows), "pareto": output["pareto_by_E4_and_serial_cost"], "one_percent_candidates": output["rows_with_E4_at_most_1pct_g_proxy"]}, indent=2))


if __name__ == "__main__":
    main()
