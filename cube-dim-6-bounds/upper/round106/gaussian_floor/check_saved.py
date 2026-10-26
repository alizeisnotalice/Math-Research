#!/usr/bin/env python3
"""UNRUN saved-output mincheck, requiring a separate primary release.

Checks provenance, exact grid coverage, and saved certified endpoint balls.
It does not rerun covariance computations, branch sums, or the original producer.
It is not an independent proof of the moment-packing reduction.
"""
import signal
signal.signal(signal.SIGALRM, signal.SIG_DFL)
signal.alarm(177)

import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
START = time.monotonic()
INPUT = HERE / "mincheck_inputs.json"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    from flint import arb, ctx
    ctx.prec = 192
    spec = json.loads(INPUT.read_text())
    for name, expected in spec["saved_sha256"].items():
        assert Path(name).name == name
        assert sha(HERE/name) == expected, name
    released = json.loads((HERE/"primary_release.json").read_text())
    result = json.loads((HERE/"results.json").read_text())
    receipt = json.loads((HERE/"run_receipt.json").read_text())
    assert released["authorized"] is True and released["max_runs"] == 1
    assert result["source_sha256"] == released["source_sha256"] == sha(HERE/"check_gaussian_floor.py")
    assert result["input_sha256"] == released["input_sha256"] == sha(HERE/"inputs.json")
    assert result["status"] == receipt["status"] == "PASS"
    assert receipt["exit_code"] == 0 and receipt["run_count"] == 1
    assert result["child_processes_created"] == receipt["children_created_by_program"] == 0
    rows = result["certificate_rows"]
    assert len(rows) == 320
    seen = set()
    for row in rows:
        k, part = row["k"], row["part"]
        assert isinstance(k, int) and 1 <= k <= 5
        assert isinstance(part, int) and 0 <= part < 64
        assert (k, part) not in seen
        seen.add((k, part))
        assert list(map(Q, row["interval"])) == [Q(part,64),Q(part+1,64)]
        assert row["proved_gt_1_over_20"] is True
        assert arb(row["branch_lower_ball"]) > arb(1)/20
    assert seen == {(k,part) for k in range(1,6) for part in range(64)}
    cap = arb(result["wholecube_samplewise_floor_ceiling_ball"])
    assert 0 < cap < arb(1)/165888 < arb(1)/65000
    example = result["deterministic_orbit_example"]
    assert example["checked_closed_subintervals"] == 320
    assert example["strict_gap_squared_integer_comparison"] == [400*601,6**5*101]
    assert 400*601 < 6**5*101
    assert arb(example["h_sym_strict_lower_ball"]) > arb(example["original_pointwise_min_ceiling_ball"])
    assert "random expectation certificate for orbit-floor" in result["not_claimed"]
    assert "g0 proof" in result["not_claimed"] and "new C6 bound" in result["not_claimed"]
    return {"status":"PASS_saved_support_mincheck", "unique_closed_intervals":320,
            "recomputed_covariance_or_branch_sums":False,
            "independent_mathematical_proof":False,
            "saved_outputs_bound_to_released_source_and_input":True,
            "new_C6_claim":False}


def main():
    code = 0
    try:
        result = run()
    except BaseException as exc:
        code = 1
        result = {"status":"FAIL", "exception_type":type(exc).__name__,"exception":str(exc)}
    result.update({"input_sha256":sha(INPUT),"source_sha256":sha(Path(__file__).resolve()),
                   "wall_seconds":time.monotonic()-START,"child_processes_created":0,
                   "alarm_seconds":177,"alarm_handler":"SIG_DFL","alarm_cancelled":False,
                   "alarm_policy":"retained finite short-feedback choice, not a user runtime ceiling"})
    (HERE/"mincheck_results.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
    return code


if __name__ == "__main__":
    sys.exit(main())
