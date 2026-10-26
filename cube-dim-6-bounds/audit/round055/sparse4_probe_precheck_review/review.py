"""One registered independent structural batch; never imports a producer module."""
import sys
sys.dont_write_bytecode = True
import datetime, gzip, hashlib, json, os, re, time
from pathlib import Path
from fractions import Fraction as Q
from independent_structure import (
    Budget, canonical_line, check_frame, check_pairs, check_group, check_scaled_fan,
)
P = Path(__file__).resolve().parent
def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read_json(path):
    return json.loads(Path(path).read_text())
def verify_pins(pins):
    for path, digest in pins.items():
        assert sha(path) == digest, ("source hash changed", path)
def parse_integer(token):
    assert re.fullmatch(r"0|-?[1-9][0-9]*", token), ("noncanonical integer", token)
    value = int(token)
    assert -(1 << 63) < value < (1 << 63)
    return value

def main():
    assert sys.argv[1:] == ["--run"]
    assert not (P / "execution_receipt.json").exists(), "one attempt; never restart"
    cfg = read_json(P / "input_manifest.json")
    assert cfg["status"] == "STABLE_PRODUCER_INPUT_PINNED_FOR_ONE_STRUCTURAL_REVIEW"
    verify_pins(cfg["sha256"])
    for f, h in cfg["auditor_source_sha256"].items():
        assert sha(P / f) == h
    receipt = {
        "status": "RUNNING", "pid": os.getpid(),
        "actual_started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "input_manifest_sha256": sha(P / "input_manifest.json"),
        "production_native_calls": 0, "author_writer_runs": 0,
        "old_formal_repeats": 0, "matrix_operation_cap": 30000,
    }
    (P / "execution_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    start = time.monotonic()
    budget = Budget()
    try:
        frame = read_json(cfg["frame"])
        fan = json.loads(gzip.decompress(Path(cfg["fan"]).read_bytes()))
        producer = Path(cfg["producer"])
        register = read_json(producer / "registration.json")
        assert register["point_numerators"] == [50, 50, 50, 45, 27, 0]
        assert register["denominator"] == 50
        assert register["T_exact"] == "32/825" and register["beta_exact"] == "1/80"
        directions, columns, frame_result = check_frame(frame, budget)
        pair_result = check_pairs(
            directions, register["point_numerators"], 50,
            read_json(producer / "S2_certificate.json"), budget,
        )
        representatives, group_result = check_group(register["point_numerators"], budget)
        rays, cones, scale_result = check_scaled_fan(fan, columns, budget)

        raw = (producer / "native_input.txt").read_text().split()
        assert raw[0] == "C055_SPARSE4_BETA1_80_X050_T32_825"
        values = list(map(parse_integer, raw[1:]))
        assert values[:7] == [478, 13229, 1920, 50, 24, 32, 825]
        offset = 7
        for ray in rays:
            assert values[offset:offset + 7] == [ray["Q"]] + ray["R"]
            offset += 7
        for cone in cones:
            assert values[offset:offset + 7] == [cone["scaled_determinant"]] + cone["ray_ids"]
            offset += 7
        images = [tuple(values[offset + 6 * i:offset + 6 * (i + 1)]) for i in range(1920)]
        offset += 6 * 1920
        assert offset == len(values)
        assert len(set(images)) == 1920
        assert {canonical_line(x) for x in images} == set(representatives)
        assert all(x < tuple(-a for a in x) for x in images)
        assert sum(x[0] == 0 for x in images) == 320
        assert all(x == 0 or abs(x) in (27, 45, 50) for row in images for x in row)

        # Reconstruct the exact canonical serialization used by the adapter,
        # without importing or running it.
        packed = {
            "rays": [{"R": ray["R"], "Q": ray["Q"],
                      "scale": Q(fan["rays"][i]["q_true_152"]).denominator}
                     for i, ray in enumerate(rays)],
            "cones": [[cone["scaled_determinant"], cone["ray_ids"]] for cone in cones],
            "point_orbit_representatives": images,
            "orbit_weight": 24, "point_denominator": 50,
        }
        canonical_digest = hashlib.sha256(
            json.dumps(packed, separators=(",", ":"), sort_keys=True).encode()
        ).hexdigest()
        cert = read_json(producer / "input_certificate.json")
        assert canonical_digest == cert["canonical_data_sha256"]
        assert cert["input_sha256"] == sha(producer / "native_input.txt")
        assert len(images) * len(cones) == 25399680
        assert len(images) * len(cones) * 24 == 609592320
        assert 2 * len(images) * len(cones) == 50799360

        time_exact = Q(32, 825)
        assert time_exact / 4096 == Q(1, 105600)
        assert time_exact / 4000 == Q(1, 103125)
        # Coefficients of pi powers: (pi/2)/(64*pi^6)=1/(128*pi^5).
        assert Q(1, 2) / 64 == Q(1, 128)
        normalization = {
            "H_normalizer": "64*pi^6",
            "core_G0_nonzero_value": "(pi/2)*sign(b)",
            "core_G0_normalizer": "64*pi^6",
            "equivalent_sign_Gamma_normalizer": "128*pi^5",
            "zero_node_value": "0",
            "constant_functional_all_higher_normalized_jets": "0",
            "derivatives_of_spatial_G_claimed": False,
            "T_exact": str(time_exact),
            "fixed_time_64_necessary_threshold": str(time_exact / 4096),
            "fixed_time_62_5_necessary_threshold": str(time_exact / 4000),
        }

        build = read_json(producer / "build_record.json")
        assert sha(producer / "sparse4_probe.cpp") == build["source_sha256"]
        assert sha(producer / "scratch/sparse4_probe") == build["binary_sha256"]
        for path, h in build["runtime_dylib_sha256"].items():
            assert sha(path) == h
        assert sha(build["headers_archive"]) == build["headers_sha256"]
        source = (producer / "sparse4_probe.cpp").read_text()
        launcher = (producer / "launch_once.py").read_text()
        required_source = [
            "static constexpr slong PREC=512",
            "arb_zero(g[j].a)",
            "if(arb_contains_zero(d.a))",
            "terms!=25399680ULL||terms*24!=609592320ULL",
            "arb_mul_2exp_si(scale.a,scale.a,1);arb_pow_ui(scale.a,scale.a,6,PREC)",
            "arb_sub(K.a,G.a,H.a,PREC)",
            "if(next(in)!=478||next(in)!=13229||next(in)!=1920)",
            "if(arb_contains_zero(d.a))throw",
        ]
        for fragment in required_source:
            assert fragment in source, fragment
        for fragment in [
            "MAIN_RELEASED_ONE_NATIVE_AFTER_B055", "B055_review_accepted_by_main",
            "PASS_STRICT_OFF_TRUE_S2", "child.wait(timeout=1800)",
            "Refuse second native attempt", "F050_AND_F052_FINAL_PASS_PINNED",
            "global_mT_certified=False", "new_Cn_bound=False",
        ]:
            assert fragment in launcher, fragment
        # Source inspection supplies semantics; these assertions only bind
        # reviewed invariants to the exact source, and are not an independent
        # proof of arbitrary C++ behavior.
        driver_result = {
            "source_review_scope": "manual complete source review plus frozen-string binding",
            "native_input_all_records_matched": True,
            "input_header": raw[0], "canonical_input_digest": canonical_digest,
            "paired_tuple_traversals": 25399680,
            "G_logical_terms": 609592320, "H_logical_terms": 609592320,
            "G_DD5_evaluations": 25399680, "H_DD5_evaluations": 25399680,
            "joint_field_DD5_evaluations": 50799360,
            "precision_bits": 512, "launcher_native_timeout_seconds": 1800,
            "source_sha256": build["source_sha256"],
            "binary_sha256": build["binary_sha256"],
            "native_executed_by_auditor": False,
        }
        verify_pins(cfg["sha256"])
        budget.guard()
        result = {
            "status": "PASS_PRECHECK_STRUCTURE_AND_SOURCE_ONLY",
            "configurations": 6, "frame": frame_result, "true_S2": pair_result,
            "group": group_result, "full_fan_scaling": scale_result,
            "normalization": normalization, "driver": driver_result,
            "matrix_operations": budget.matrix_operations,
            "production_G_H_p_native": 0, "author_writer_runs": 0,
            "old_formal_repeats": 0, "new_C6_bound": False,
            "production_launch_authorized_by_this_audit": False,
        }
        (P / "results.json").write_text(json.dumps(result, indent=2) + "\n")
        receipt.update(status=result["status"], exit_code=0,
                       matrix_operations=budget.matrix_operations, configurations=6,
                       actual_seconds=time.monotonic() - start,
                       results_sha256=sha(P / "results.json"))
    except BaseException as error:
        receipt.update(status="STOP_AUDIT_FAILURE_NO_RETRY", exit_code=2,
                       error=repr(error), matrix_operations=budget.matrix_operations,
                       actual_seconds=time.monotonic() - start)
        raise
    finally:
        receipt["completed_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        receipt["all_handles_closed"] = True
        (P / "execution_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({k: v for k, v in receipt.items() if k != "input_manifest_sha256"}))

if __name__ == "__main__":
    main()
