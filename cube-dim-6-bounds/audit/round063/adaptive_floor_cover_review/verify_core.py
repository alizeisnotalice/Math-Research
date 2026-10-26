#!/usr/bin/env python3
"""Seven independent pure-state scenarios, with no native or author-test execution."""
import signal
signal.alarm(180)
import sys
sys.dont_write_bytecode = True
from fractions import Fraction as Q
from pathlib import Path
from itertools import product
from tempfile import TemporaryDirectory
from datetime import datetime, timezone
import copy
import hashlib
import importlib.util
import json
import os
import time

ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rejects(call):
    try:
        call()
    except ValueError as error:
        return str(error)
    raise AssertionError("Expected rejection did not occur")


def main():
    tic = time.monotonic()
    start = datetime.now(timezone.utc).isoformat()
    assert not (ROOT/"results.json").exists(), "Refuse to overwrite the completed audit"
    pins = json.loads((ROOT/"input_manifest.json").read_text())
    for path,digest in pins["source_sha256"].items():
        assert sha(path) == digest, path
    spec = importlib.util.spec_from_file_location("f063_reviewed_core", pins["author_core"])
    C = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(C)
    pin = pins["author_input_manifest_sha256"]
    fresh = lambda:C.new_state(pin, test_only=True)
    checks = []

    def fixture_jet(attempt, value):
        expected = attempt["input"]
        return {**{key:copy.deepcopy(expected[key]) for key in ("model","model_definition","window","center")},
                "F_interval":[str(value),str(value)],"face_gradient_intervals":[["0","0"] for _ in range(5)],
                "source_proof":"F063_SYNTHETIC_CONSTANT_FUNCTION_ONLY_NO_BODY32_OBSERVATION",
                "provenance":{"kind":"C016_NATIVE_SAME_WINDOW_POSITIVE_BAND","precision_bits":512,
                              "full_group_terms":967680,"computed_sign_pairs":483840,"synthetic_fixture":True,
                              "native_sha256":"1"*64,"receipt_sha256":"2"*64,"gradient_gate_sha256":"3"*64}}

    def fixture_observation(attempt, value=Q(1), *, pass_box=False):
        jet = fixture_jet(attempt,value)
        observation={"input_sha256":attempt["input_sha256"],"receipt_exit_code":0,"stderr_empty":True,"jet":jet}
        if value > C.TARGET:
            expected=attempt["input"]
            center=list(map(Q,expected["center"]))
            bounds=C.box(expected["face_box"])
            a=Q(expected["window"]["a"])
            rows=[]
            for bits in product((0,1),repeat=5):
                delta=[bounds[i][bit]-center[i+1] for i,bit in enumerate(bits)]
                quadratic=sum(di*di for di in delta)/32
                exponent=-quadratic/a**2
                # Unresolved: 0 <= value*exp(exponent) <= value, since exponent<=0.
                # Passing root: value=2^100 and exponent=-40; exp(1)<4 proves
                # value*exp(-40)>2^20>1. No exponential evaluator is called.
                if pass_box:
                    assert value == Q(2)**100 and exponent == -40
                lower="1" if pass_box else "0"
                rows.append({"corner_bits":list(bits),"displacement":list(map(str,delta)),
                             "gradient_dot_lower":"0","minimizing_center_value":str(value),
                             "quadratic":str(quadratic),"exponent":str(exponent),
                             "certificate_lower":lower,"certificate_upper":str(value)})
            certificate={"model":C.MODEL,"window":copy.deepcopy(expected["window"]),"quadratic":"frame32",
                         "target":"1/65000","conditional_on_jet_source":jet["source_proof"],"corners":rows,
                         "tangent_bound_lower":"1" if pass_box else "0","tangent_bound_upper":str(value),
                         "status":"PASS_CONDITIONAL_ON_MATCHING_CENTER_JET" if pass_box else "UNRESOLVED"}
            observation.update(box_certificate=certificate,box_certificate_sha256=C.digest(certificate),
                               box_checker_sha256=C.BOX_CHECKER_SHA256)
        return observation

    # 1. Independent asymmetric rational geometry, no sampled grid.
    raw=[["1/8","7/8"],["1/4","3/4"],["0","1/2"],["1/3","2/3"],["0","1/4"]]
    tightened=C.tighten(raw)
    expected=[["1/3","7/8"],["1/3","3/4"],["1/3","1/2"],["1/3","1/2"],["0","1/4"]]
    assert tightened["box"] == expected and C.tighten(expected) == tightened
    center=list(map(Q,C.midpoint(expected)[1:]))
    assert all(center[i]>=center[i+1] for i in range(4))
    split=C.split_box(expected)
    assert split["axis"]==0 and split["seam"]=="29/48"
    seam=[Q(29,48),Q(1,2),Q(1,3),Q(1,3),Q(0)]
    for child in [split["lower"],split["upper"]]:
        assert not child["empty"]
        assert all(lo<=x<=hi for (lo,hi),x in zip(C.box(child["box"]),seam))
    checks.append({"configuration":1,"status":"PASS_PURE_GEOMETRY_CONFIGURATION",
                   "raw_box":raw,"tightened_box":expected,"seam_witness":list(map(str,seam)),
                   "general_proof":"proof_review.md"})

    # 2. Correctly hash-linked public appends accept a smaller unresolved root.
    initial=fresh()
    smaller=copy.deepcopy(initial)
    smaller["nodes"]["r"]=C.make_node("r",None,0,[["0","1/2"] for _ in range(5)],0)
    accepted=C.validate_state(smaller)
    assert accepted["root_partition_preserved"]
    with TemporaryDirectory(prefix="root_",dir=ROOT) as directory:
        before=C.append_snapshot(directory,initial,None)
        after=C.append_snapshot(directory,smaller,before)
        recovered=C.recover(directory)
        assert recovered["sequence"]==1 and recovered["sha256"]==after
        assert recovered["state"]["nodes"]["r"]["box"]==[["0","1/2"] for _ in range(5)]
    witness=["1"]*5
    assert not all(lo<=Q(x)<=hi for (lo,hi),x in zip(C.box(smaller["nodes"]["r"]["box"]),witness))
    checks.append({"configuration":2,"status":"BLOCKER_ACCEPTS_ROOT_COVERAGE_LOSS",
                   "public_api_path":["append_snapshot(initial)","append_snapshot(smaller_root,current_hash)","recover"],
                   "previous_snapshot_sha256":before,"accepted_snapshot_sha256":after,
                   "missed_ordered_face_witness":witness,"validator_claim":accepted,
                   "raw_snapshot_bytes_modified":False,"hashes_forged":False})

    # 3. A stale state object is accepted with the current predecessor hash.
    initial=fresh()
    inflight,original=C.reserve(initial)
    with TemporaryDirectory(prefix="orphan_",dir=ROOT) as directory:
        h0=C.append_snapshot(directory,initial,None)
        h1=C.append_snapshot(directory,inflight,h0)
        h2=C.append_snapshot(directory,initial,h1)
        recovered=C.recover(directory)
        resumed,reissued=C.reserve(recovered["state"])
        assert reissued["token"]==original["token"] and reissued["input_sha256"]==original["input_sha256"]
        assert recovered["state"]["inflight"] is None
        assert recovered["state"]["batch"]["native_reservations"]==0
    checks.append({"configuration":3,"status":"BLOCKER_ACCEPTS_ORPHAN_ROLLBACK_AND_REISSUE",
                   "public_api_path":["reserve","append_snapshot(inflight)","append_snapshot(stale_initial,current_hash)","recover","reserve"],
                   "inflight_snapshot_sha256":h1,"accepted_rollback_sha256":h2,
                   "original_token":original["token"],"reissued_token":reissued["token"],
                   "same_exact_input_sha256":original["input_sha256"],
                   "reservation_count_before":1,"reservation_count_after_rollback":0,
                   "reservation_count_after_reissue":resumed["batch"]["native_reservations"],
                   "recovery_action":recovered["recovery_action"],"native_launches":0,
                   "raw_snapshot_bytes_modified":False,"hashes_forged":False})

    # 4. Ordinary API transitions correctly enforce the retained first-batch cap.
    state=fresh()
    started=[]
    for _ in range(8):
        state,attempt=C.reserve(state)
        started.append(state["inflight"]["node_id"])
        state=C.ingest(state,attempt["token"],fixture_observation(attempt))
    assert len(set(started))==8 and state["batch"]["native_reservations"]==8
    ninth_error=rejects(lambda:C.reserve(state))
    assert "budget exhausted" in ninth_error
    checks.append({"configuration":4,"status":"PASS_RETAINED_HISTORY_EIGHT_NODE_EIGHT_ATTEMPT_CAP",
                   "started_nodes":started,"reservations":8,"ninth_reservation_error":ninth_error,
                   "synthetic_only":True,"actual_native_calls":0,"global_g_certificate":state["global_g_certificate"]})

    # 5. Reusing a previous cutoff jet does not pass the actual adapter.
    state,old_attempt=C.reserve(fresh())
    state=C.ingest(state,old_attempt["token"],fixture_observation(old_attempt,C.TARGET/2))
    assert state["nodes"]["r"]["cutoff_index"]==1
    state,new_attempt=C.reserve(state)
    wrong=fixture_observation(old_attempt,C.TARGET/2)
    wrong["input_sha256"]=new_attempt["input_sha256"]
    error=rejects(lambda:C.ingest(state,new_attempt["token"],wrong))
    assert "window mismatch" in error
    assert len(state["nodes"]["r"]["attempts"])==2
    checks.append({"configuration":5,"status":"PASS_SAME_WINDOW_GUARD",
                   "old_a":old_attempt["input"]["window"]["a"],"new_a":new_attempt["input"]["window"]["a"],
                   "reuse_rejection":error,"prior_attempt_retained":True,"synthetic_only":True})

    # 6. Passed leaves are protected by both scheduling and snapshot append.
    initial=fresh()
    state,attempt=C.reserve(initial)
    passed=C.ingest(state,attempt["token"],fixture_observation(attempt,Q(2)**100,pass_box=True))
    assert C.validate_state(passed)["all_nonempty_leaves_passed"] and not passed["global_g_certificate"]
    dispatch_error=rejects(lambda:C.reserve(passed))
    with TemporaryDirectory(prefix="passed_",dir=ROOT) as directory:
        h0=C.append_snapshot(directory,initial,None)
        h1=C.append_snapshot(directory,passed,h0)
        changed=copy.deepcopy(passed)
        changed["nodes"]["r"]["proof"]["token"]="changed"
        mutation_error=rejects(lambda:C.append_snapshot(directory,changed,h1))
    checks.append({"configuration":6,"status":"PASS_PASSED_LEAF_NO_RECOMPUTE",
                   "dispatch_rejection":dispatch_error,"snapshot_change_rejection":mutation_error,
                   "synthetic_only":True,"all_nonempty_leaves_passed_is_not_a_global_g_certificate":True})

    # 7. The public append/recover path accepts a test-to-production gate change.
    initial=fresh()
    changed=copy.deepcopy(initial)
    changed["environment"]="MAIN_RELEASED"
    assert changed["release"] is None
    with TemporaryDirectory(prefix="release_",dir=ROOT) as directory:
        h0=C.append_snapshot(directory,initial,None)
        h1=C.append_snapshot(directory,changed,h0)
        recovered=C.recover(directory)
        reserved,attempt=C.reserve(recovered["state"])
    checks.append({"configuration":7,"status":"BLOCKER_ACCEPTS_MISSING_RELEASE_AFTER_ENVIRONMENT_CHANGE",
                   "public_api_path":["append_snapshot(TEST_ONLY)","append_snapshot(environment_MAIN_RELEASED_release_None,current_hash)","recover","reserve"],
                   "accepted_snapshot_sha256":h1,"release":reserved["release"],
                   "reserved_token":attempt["token"],"actual_native_calls":0,
                   "external_provenance_forged":False,"raw_snapshot_bytes_modified":False,"hashes_forged":False})

    assert len(checks)==7
    for path,digest in pins["source_sha256"].items():
        assert sha(path)==digest,path
    output={"status":"REJECT_IMPLEMENTATION_FOR_RELEASE_PURE_GEOMETRY_RETAINED",
            "actual_configurations":7,"checks":checks,"concrete_blocking_scenarios":[2,3,7],
            "native_calls":0,"kernel_calls":0,"author_formal_program_runs":0,"old_formal_replays":0,
            "matlab_calls":0,"subagents":0,"global_g_certificate":False,"new_C6_bound":False,
            "author_files_unchanged":True,"pid":os.getpid(),"python":sys.version,
            "start_utc":start,"end_utc":datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds":time.monotonic()-tic,"hard_process_limit_seconds":180,
            "source_sha256":{name:sha(ROOT/name) for name in ["registration.md","input_manifest.json","verify_core.py"]}}
    (ROOT/"results.json").write_text(json.dumps(output,indent=2)+"\n")
    print(json.dumps({key:output[key] for key in ["status","actual_configurations","concrete_blocking_scenarios","elapsed_seconds","native_calls"]},indent=2))


if __name__=="__main__":
    main()
