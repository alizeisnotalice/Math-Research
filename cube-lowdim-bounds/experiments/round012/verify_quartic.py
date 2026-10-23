#!/usr/bin/env python3
"""Independently verify saved quartic leaves; no adaptive search or writes.

Only the previously reviewed round008 Fraction arithmetic/log/tail/b2 kernel
is imported. The quartic certifier itself is never imported or called.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent


def require(condition, label):
    if not condition:
        raise ValueError(label)


def core_path(explicit):
    if explicit is not None:
        return explicit.resolve()
    options = [
        ROOT.parent / "round008" / "certify.py",
        ROOT.parent.parent / "Math-Research-lowdim-agent" / "cube-lowdim-bounds"
        / "experiments" / "round008" / "certify.py",
    ]
    for path in options:
        if path.is_file():
            return path.resolve()
    raise FileNotFoundError("Pass --core <round008/certify.py>.")


def verify(path):
    data = json.loads((ROOT / "certificate.json").read_text())
    require(sha256(path.read_bytes()).hexdigest() == data["core_sha256"], "core hash")
    require(sha256((ROOT / "certify_quartic.py").read_bytes()).hexdigest()
            == data["certifier_sha256"], "certifier hash")
    for name, recorded_hash in data["analytic_dependencies"].items():
        require(sha256((ROOT / name).read_bytes()).hexdigest() == recorded_hash,
                "analytic dependency hash: " + name)
    require(set(data["analytic_dependencies"])
            == {"quartic_angular.md", "small_radius.json", "transfer_extension.md"},
            "required analytic dependencies")
    require(data["limits"]["tail_terms"] == 28 and data["limits"]["log_terms"] == 72,
            "kernel truncation metadata")
    attempts = data["attempts"]
    require(1 <= len(attempts) <= 2, "attempt count")
    successful = [a for a in attempts
                  if a["status"] == "certified_under_quartic_angular_and_transfer"]
    require(len(successful) == 1, "one successful saved certificate")
    saved = successful[0]
    parameters = saved["parameters"]
    require(parameters == {"R": "379/200", "e2": "3/4", "e3": "0", "e4": "21/100"},
            "intended rational candidate")
    R, e2, e4 = map(Q, (parameters["R"], parameters["e2"], parameters["e4"]))
    require(R > 1 and e2 >= 0 and e4 >= 0, "parameter domain")
    spec = importlib.util.spec_from_file_location("quartic_verifier_fraction_core", path)
    core = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = core
    spec.loader.exec_module(core)
    core.R = R
    I = core.I

    def quartic_diagonal(r):
        s = 1-r
        a = 1-r/2
        return 8*(a*r.square()+Q(22,3)*(a*a*a)
                  +2*(1+s*s*s)*core.log_i(r/2)
                  -2*s.square()*core.xlogx(s))

    def curvature_cap(r):
        if r.hi <= Q(2,3):
            rational = -Q(5,2)*r-3+18/r-4/r.square()
        elif r.lo >= Q(2,3):
            rational = (165-101*r)/8+1/(2*r.square())
        else:
            raise ValueError("saved leaf straddles 2/3")
        return 12*core.log_i(r/2)+rational

    def saved_leaf_lower(r):
        diagonal = (2*(2*core.tail_series(r)+I.point(1)/R)/(R-1)
                    -e2*core.b2_diag(r)-e4*quartic_diagonal(r))
        s = 1-r
        k = (e2*64*s.square()/(2-r).square()/(2-r)
             +I.point(1)/(R**3*(R-1))-4*e4*curvature_cap(r))
        return (diagonal+min(k.lo,Q(0))*r.hi*r.hi).lo

    leaves = saved["accepted_leaves"]
    require(bool(leaves), "nonempty leaves")
    require(Q(leaves[0]["lo"]) == Q(1,10) and Q(leaves[-1]["hi"]) == 1,
            "cover endpoints")
    verified_lowers = []
    for index, leaf in enumerate(leaves):
        lo, hi, recorded = map(Q, (leaf["lo"], leaf["hi"], leaf["lower"]))
        require(Q(1,10) <= lo < hi <= 1, "positive leaf width/domain")
        if index:
            require(Q(leaves[index-1]["hi"]) == lo, "exact adjacency")
        actual = saved_leaf_lower(I(lo, hi))
        require(actual >= 0, "recomputed nonnegative leaf")
        require(recorded == actual, "saved rational leaf lower exactly reproduced")
        verified_lowers.append(actual)
    inner = saved["inner"]
    require(inner["accepted_boxes"] == len(leaves), "saved leaf count")
    require(inner["splits"]+2 == len(leaves), "binary leaf/split count")
    require(inner["unresolved_boxes"] == 0 and not inner["unresolved_intervals"],
            "no unresolved boxes")
    require(inner["status"] == "certified" and inner["cover_verified"],
            "saved coverage metadata")
    minimum = min(verified_lowers)
    require(minimum == Q(inner["minimum_accepted_lower"]), "exact minimum lower")
    text = "\n".join(f"{leaf['lo']}/{leaf['hi']}" for leaf in leaves)
    require(sha256(text.encode()).hexdigest() == inner["accepted_leaf_sha256"],
            "saved leaf hash")

    # Endpoint-strip signs use only exact outward log enclosures.
    h4 = 8*(Q(553,150)+Q(729,1000)*core.log_point(Q(1,9))
            +core.log_point(Q(1,11)))
    quadratic = 3-5*core.log_point(2)
    require(h4.hi < 0 and quadratic.hi < 0, "small-radius strict signs")
    small = saved["small_radius"]
    supplied = json.loads((ROOT / "small_radius.json").read_text())
    require(h4.pair() == supplied["h4_uniform_upper_interval"]
            == small["h4_upper_interval"], "small-radius exact intervals")
    require(quadratic.hi == Q(small["quadratic_T_half_upper"]),
            "quadratic small-radius bound")
    require(small["h4_strictly_negative"]
            and small["quadratic_T_half_strictly_negative"]
            and small["matches_supplied_small_radius_interval"], "small metadata")

    outer = 1-(e2+e4)*R*(R-1)*(2-2*core.log_point(2))
    dominance = R-(e2+e4)*(R-1)/4
    previous = Q(6938554331,1790000000)
    cost = R*R/(R-1)-e2/6-e4/15
    checks = saved["outer_and_cost"]
    require(outer.lo > 0 and dominance > 0 and cost < previous,
            "outer/dominance/improvement strict signs")
    require(outer.pair() == checks["outer_slack"], "outer exact interval")
    require(dominance == Q(checks["kernel_dominance_margin"]), "dominance fraction")
    require(cost == Q(checks["cost"]) and previous == Q(checks["previous_upper"]),
            "exact cost and incumbent")
    require(previous-cost == Q(checks["improvement"]), "exact improvement")
    require(checks["strict"], "strict metadata")
    require(data["status"] == saved["status"], "top-level status")
    return {
        "status": "all_saved_rational_leaves_verified",
        "method": "no adaptive search; independent leaf formula; reused reviewed Fraction kernel",
        "leaf_count": len(leaves),
        "exact_cover": ["1/10", "1"],
        "all_recorded_lowers_exactly_reproduced": True,
        "minimum_lower_display": float(minimum),
        "small_h4_upper_display": float(h4.hi),
        "outer_slack_lower_display": float(outer.lo),
        "kernel_dominance": str(dominance),
        "cost": str(cost),
        "improvement": str(previous-cost),
        "source_and_analytic_hashes_match": True,
        "scope": "scalar certificate verified; quartic angular and transfer proofs remain explicit dependencies",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--core", type=Path)
    arguments = parser.parse_args()
    print(json.dumps(verify(core_path(arguments.core)), indent=2))
