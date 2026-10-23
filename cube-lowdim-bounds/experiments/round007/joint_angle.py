"""Round 7: exact angular formula check for Lin's single-bump estimates.

Uses closed formulas extracted from BaseDensity.lean and BumpGenerator.lean.
It performs 401 angle evaluations (no quadrature) at the fixed r0 certificate radius.
The output is numerical evidence only; equality at the diagonal is established separately
from the source formulas.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
R = 4746319 / 2500000
EPSILON = 4981349 / 6250000
R0 = 86425249 / 100000000
N = 401


def tail_f(c: float) -> float:
    if c == 0.0:
        return 1 / (2 * R)
    return -(R / c**2) * math.log1p(-c / R) - 1 / c


def jfun(s: float) -> float:
    return 0.5 - 2 / s - 4 / s**2 * math.log1p(-s / 2)


def ffun(t: float) -> float:
    return (2 - (1 - t) * math.log(abs(1 - t))
            - (1 + t) * math.log(1 + t) + 2 * math.log(t))


def tfun(r: float) -> float:
    return 2 * (2 - r) * (1 - math.log(2) + math.log(r)) - 2 * (1 - r) * math.log1p(-r)


def tail_generator(r: float, d: float) -> float:
    return 2 * (2 * tail_f(r) + tail_f(d) + tail_f(-d)) / (R - 1)


def bump_generator(r: float, d: float) -> float:
    u, v = (r + d) / 2, (r - d) / 2
    return 4 * ((1 - v) * ffun(u / (1 - v)) + (1 - u) * ffun(v / (1 - u)))


def main() -> None:
    tail_floor = 4 / (R * (R - 1)) * jfun(2 * R0 / R)
    bump_ceiling = 4 * tfun(R0)
    scalar_net = tail_floor - EPSILON * bump_ceiling
    rows = []
    for k in range(N):
        d = R0 * k / N  # avoids the axis, includes the diagonal exactly
        tail = tail_generator(R0, d)
        bump = bump_generator(R0, d)
        rows.append({
            "d_over_r": d / R0,
            "tail_surplus": tail - tail_floor,
            "bump_deficit": bump_ceiling - bump,
            "joint_surplus": tail - EPSILON * bump - scalar_net,
        })
    out = {
        "R": R,
        "epsilon": EPSILON,
        "r0": R0,
        "samples": N,
        "tail_floor": tail_floor,
        "bump_ceiling": bump_ceiling,
        "scalar_net_at_diagonal": scalar_net,
        "minimum_tail_surplus": min(rows, key=lambda x: x["tail_surplus"]),
        "minimum_bump_deficit": min(rows, key=lambda x: x["bump_deficit"]),
        "minimum_joint_surplus": min(rows, key=lambda x: x["joint_surplus"]),
        "representatives": [rows[k] for k in (0, 100, 200, 300, 400)],
        "interpretation": "numerical check of closed formulas; source formulas prove both bounds attain equality at d=0",
    }
    (ROOT / "joint_angle_result.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
