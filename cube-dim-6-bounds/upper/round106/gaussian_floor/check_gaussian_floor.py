#!/usr/bin/env python3
"""One bounded exact/Arb verification; run only after primary static approval.

No subprocesses, sampling, old kernel producer, or frame-subset enumeration.
The SIG_DFL 177 second alarm stays armed through JSON receipt and exit.
"""
import signal
signal.signal(signal.SIGALRM, signal.SIG_DFL)
signal.alarm(177)

import hashlib
import itertools
import json
import math
from fractions import Fraction as Q
from pathlib import Path
import sys
import time

START = time.monotonic()
HERE = Path(__file__).resolve().parent
SOURCE = Path(__file__).resolve()
INPUT = HERE / "inputs.json"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    import sympy as sp
    from flint import arb, ctx
    ctx.prec = 192
    spec = json.loads(INPUT.read_text())
    assert spec["deterministic_orbit_example"]["uniform_y_subintervals"] == 64

    # Direct substitution u=Qx/(4t^2): t^-6 dt = -16 Qx^-5/2 u^3/2 du.
    C = sp.simplify(16 * sp.gamma(sp.Rational(5, 2)) / (4 * sp.pi)**3)
    assert sp.simplify(C - sp.gamma(sp.Rational(5, 2)) / (4 * sp.pi**3)) == 0
    assert sp.simplify(C - 3 / (16 * sp.pi**sp.Rational(5, 2))) == 0
    moment = sp.Rational(1, 32) * 2 / sp.sqrt(sp.pi)
    cap = sp.simplify(C * moment / 6**sp.Rational(5, 2))
    expected_cap = 3 / (256 * sp.pi**3 * 6**sp.Rational(5, 2))
    assert sp.simplify(cap - expected_cap) == 0
    assert Q(1, 165888) < Q(1, 65000)
    # pi>3 and sqrt(6)>2 imply cap < 1/165888.
    assert Q(3, 256 * 27 * 36 * 2) == Q(1, 165888)

    # Verify the actual body32 tight-frame identity, in exact integers.
    V = [(1,) + tail for tail in itertools.product((-1, 1), repeat=5)]
    gram_sum = [[sum(v[i]*v[j] for v in V) for j in range(6)] for i in range(6)]
    assert gram_sum == [[32 if i == j else 0 for j in range(6)] for i in range(6)]
    v = sp.ones(6, 1)
    Sigma = sp.eye(6) + 100 * v * v.T
    inv = sp.eye(6) - sp.Rational(100, 601) * v * v.T
    assert Sigma * inv == sp.eye(6)
    assert Sigma.det() == 601
    assert Sigma.trace() == 606

    def ball(r):
        return arb(r.numerator) / arb(r.denominator)

    rho = Q(100, 601)
    target = arb(1) / 20
    certificate_rows = []
    min_ball = None
    min_index = None
    for k in range(1, 6):
        for part in range(64):
            left, right = Q(part, 64), Q(part+1, 64)
            lower = arb(0)
            for ell in range(k+1):
                j = k - 2*ell
                for sign in (-1, 1):
                    # Qx(y)=k+y^2-rho*(j+sign*y)^2 is convex.
                    # Its exact maximum on this closed interval is an endpoint.
                    endpoints = [Q(k)+y*y-rho*(Q(j)+sign*y)**2 for y in (left, right)]
                    q_upper = max(endpoints)
                    assert q_upper > 0
                    qb = ball(q_upper)
                    term = arb(math.comb(k, ell)) / (qb*qb*qb.sqrt())
                    lower += term
            lower /= 2**(k+1)
            assert lower > target, (k, part, str(lower))
            row = {"k":k, "part":part, "interval":[str(left),str(right)],
                   "branch_lower_ball":str(lower), "proved_gt_1_over_20":True}
            certificate_rows.append(row)
            # Float is display/selection only, never used in acceptance.
            if min_ball is None or float(lower.mid()) < float(min_ball.mid()):
                min_ball, min_index = lower, [k, part]

    # h_sym > 1/(20 sqrt(601)) > 6^-5/2 / sqrt(101).
    # The second comparison is exactly equivalent after positive squaring.
    assert 400*601 < 6**5*101
    pi = arb.pi()
    cap_ball = arb(3)/(256*pi**3 * arb(6)**2 * arb(6).sqrt())
    deterministic_lower = arb(1)/(20*arb(601).sqrt())
    deterministic_original_cap = arb(1)/(arb(6)**2*arb(6).sqrt()*arb(101).sqrt())
    assert deterministic_lower > deterministic_original_cap
    return {
        "status":"PASS",
        "proof_status":"supporting exact/Arb calculations; mathematical reductions in theorem.md need primary review",
        "constant_C":str(C),
        "wholecube_samplewise_floor_ceiling":str(expected_cap),
        "wholecube_samplewise_floor_ceiling_ball":str(cap_ball),
        "cap_lt_1_over_165888_lt_target":True,
        "actual_frame_tight_identity":"sum vv^T=32 I",
        "deterministic_orbit_example":{
            "Sigma":"I+100 vv^T", "determinant":601,"T":101,
            "checked_closed_subintervals":len(certificate_rows),
            "every_branch_lower_gt":"1/20",
            "display_selected_smallest_branch_ball":str(min_ball),
            "display_selected_smallest_branch_index":min_index,
            "h_sym_strict_lower_ball":str(deterministic_lower),
            "original_pointwise_min_ceiling_ball":str(deterministic_original_cap),
            "strict_gap_squared_integer_comparison":[400*601,6**5*101]
        },
        "certificate_rows":certificate_rows,
        "not_claimed":["random expectation certificate for orbit-floor", "g0 proof", "complete I", "positive S", "new C6 bound"]
    }


def main():
    code = 0
    try:
        result = run()
    except BaseException as exc:
        code = 1
        result = {"status":"FAIL", "exception_type":type(exc).__name__,"exception":str(exc)}
    result.update({"scope":"one normalization + one deterministic matrix orbit certificate",
                   "source_sha256":sha(SOURCE),"input_sha256":sha(INPUT),
                   "wall_seconds":time.monotonic()-START,
                   "child_processes_created":0,"alarm_seconds":177,
                   "alarm_handler":"SIG_DFL","alarm_cancelled":False})
    (HERE/"results.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({key:value for key,value in result.items() if key != "certificate_rows"}, indent=2))
    return code


if __name__ == "__main__":
    sys.exit(main())
