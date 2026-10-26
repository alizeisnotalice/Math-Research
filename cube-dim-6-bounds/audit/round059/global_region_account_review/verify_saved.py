#!/usr/bin/env python3
"""Four independent saved-only A059 checks; no root finder or kernel code."""
import signal
signal.alarm(180)
from pathlib import Path
from fractions import Fraction as Q
from decimal import Decimal, localcontext
from datetime import datetime, timezone
import hashlib
import json
import os
import sys
import time

ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dyadic_upper(obj):
    exponent = int(obj.get("common_binary_exponent", obj.get("upper_binary_exponent")))
    return Q(int(obj["upper_mantissa"])) * Q(2) ** exponent


def display(value):
    with localcontext() as context:
        context.prec = 22
        return str(Decimal(value.numerator)/Decimal(value.denominator))


def main():
    tic = time.monotonic()
    start = datetime.now(timezone.utc).isoformat()
    if (ROOT/"results.json").exists():
        raise RuntimeError("Refuse to overwrite a completed independent review")
    pins = json.loads((ROOT/"input_manifest.json").read_text())
    for path,digest in pins["source_sha256"].items():
        assert sha(path) == digest, path
    paths = {key:Path(value) for key,value in pins["roles"].items()}
    read = lambda key:json.loads(paths[key].read_text())
    author = read("author_results")
    checks = []

    # 1. Saved endpoint transcription and source manifest chain only.
    cap, middle, low = read("cap_results"), read("middle_results"), read("low_results")
    low_manifest, low_audit = read("low_manifest"), read("low_audit_manifest")
    for entry in read("author_source_index")["sources"]:
        assert sha(entry["source"]) == entry["sha256"]
    assert low_manifest["whitelist_sha256"]["results.json"] == sha(paths["low_results"])
    assert low_manifest["whitelist_sha256"]["theorem.md"] == sha(paths["low_theorem"])
    assert low_audit["status"] == "PASS_FINAL_FROZEN"
    assert low_audit["source_sha256_unchanged"][str(paths["low_manifest"])] == sha(paths["low_manifest"])
    assert low_audit["all_positive_time"] and low_manifest["all_positive_time"]
    assert low_audit["mean_B_range"] == ["0","1/4"]
    bc,be,bd = map(Q,[author["paid_cap_B"],author["paid_low_B"],author["middle_only_short_B"]])
    cap_upper = dyadic_upper(cap["configurations"][5]["complete_cap_norm_upper"])
    assert cap["status"] == "COMPLETE_FULL_TIME_CAP_ONLY"
    assert cap_upper < bc == Q(2118205441,10**15)
    assert Q(low["union"]["simple_rational_union_upper"]) == be == Q(13404268290817,10**18)
    assert dyadic_upper(low["union"]["new_full_time_union_norm_upper"]) < be
    assert Q(low_audit["independent_union_interval"][1]) < be
    assert Q(low_audit["strict_rational_union_upper"]) == be
    assert middle["status"] == "COMPLETE_ALL31_UNSIGNED_MIDDLES_SHORT_UPPER"
    assert middle["shared_constants"]["rho"] == "1/64"
    assert middle["shared_constants"]["epsilon"] == "1/1000"
    assert Q(middle["union"]["simple_rational_union_upper"]) == bd == Q(4689838525803,10**18)
    assert dyadic_upper(middle["union"]["short_union_norm_upper"]) < bd
    checks.append({"configuration":1,"status":"PASS_SAVED_SCOPE_AND_ENDPOINT_CHAIN",
                   "cap_full_G_B":str(bc),"low_full_G_B":str(be),"middle_only_short_B":str(bd),
                   "cap_source_upper_strict_margin":str(bc-cap_upper),
                   "old_budget_evaluations":0})

    # 2. Powers of saved endpoints, with no new root enclosure construction.
    roots = {item["label"]:item for item in author["root_certificates"]}
    assert len(roots) == 5
    root_checks = []
    for label,item in roots.items():
        lower,upper,value = map(Q,[item["lower"],item["upper"],item["radicand"]])
        degree = item["degree"]
        assert isinstance(degree,int) and degree > 0
        assert 0 <= lower < upper and lower**degree <= value < upper**degree
        root_checks.append({"label":label,"degree":degree,"lower_power_le_input":True,
                            "input_lt_upper_power":True,"strict_lower_power":lower**degree < value})
    cr,er,pr,tr,rr = [roots[key] for key in ["cap_B_to_I","low_B_to_I","paid_disjoint_regions_B",
                                            "conditional_target_I","sufficient_remaining_B_allowance"]]
    paid_i = Q(author["paid_I_conservative_upper"])
    room = Q(author["sufficient_unpaid_I_allowance"])
    assert Q(cr["radicand"]) == bc**6 and cr["degree"] == 5
    assert Q(er["radicand"]) == be**6 and er["degree"] == 5
    assert paid_i == Q(cr["upper"])+Q(er["upper"])
    assert Q(pr["radicand"]) == paid_i**5 and pr["degree"] == 6
    assert Q(author["paid_disjoint_B_upper"]) == Q(pr["upper"])
    assert Q(tr["radicand"]) == Q(author["conditional_Bcrit"])**6 and tr["degree"] == 5
    assert room == Q(tr["lower"])-paid_i
    assert Q(rr["radicand"]) == room**5 and rr["degree"] == 6
    assert Q(author["sufficient_unpaid_B_allowance_lower"]) == Q(rr["lower"])
    checks.append({"configuration":2,"status":"PASS_FIVE_SAVED_INTEGER_POWER_CERTIFICATES",
                   "root_certificates":root_checks,"serialized_dependency_links":True,"root_search_calls":0})

    # 3. Consume the saved witness rows; do not enumerate new spatial cases.
    rows = author["boolean_partition_cases"]
    assert len(rows) == 16
    witnessed = set()
    for row in rows:
        bits = tuple(row[key] for key in ["C","E","D","H"])
        assert all(type(bit) is bool for bit in bits)
        witnessed.add(bits)
        first = next((i for i,bit in enumerate(bits) if bit),4)
        assert row["part"] == first
    assert len(witnessed) == 16  # all possible four-bit rows are now covered.
    assert author["counts"] == [1,5,10,10,5] and sum(author["counts"]) == 31
    checks.append({"configuration":3,"status":"PASS_SAVED_PRIORITY_WITNESSES",
                   "distinct_saved_boolean_rows":16,"coverage_is_proved_in":"proof_review.md",
                   "new_spatial_samples":0,"pair_enumerations":0})

    # 4. Exact interpretation and monotone sufficient-allowance relations.
    crit,g = Q(author["conditional_Bcrit"]),Q(author["conditional_global_g"])
    paid_b,room_b = Q(author["paid_disjoint_B_upper"]),Q(author["sufficient_unpaid_B_allowance_lower"])
    assert crit == Q(4,48750) and g == Q(1,65000)
    assert (6*crit)**6 == 64**5*g**6
    assert room > 0 and room_b > 0 and room_b**6 <= room**5
    assert Q(tr["lower"])**5 < crit**6
    assert paid_b**6 >= paid_i**5
    assert author["checker_sha256"] == sha(paths["author_checker"])
    for key in ["global_I_certified","global_g_certified","new_C6_bound"]:
        assert author[key] is False
    assert author["kernel_or_native_calls"] == 0 and author["source_production_replays"] == 0
    assert author["unpaid_terms"] == ["middle long-time norm L_D","high-endpoint integral I_H",
                                       "far-region integral I_F","global Green floor g"]
    checks.append({"configuration":4,"status":"PASS_CONDITIONAL_ALLOWANCE_INTERPRETATION",
                   "paid_I_upper":str(paid_i),"paid_B_upper":str(paid_b),
                   "conditional_Bcrit":str(crit),"sufficient_unpaid_I_allowance":str(room),
                   "sufficient_unpaid_B_allowance_lower":str(room_b),
                   "display_only":{"paid_B":display(paid_b),"Bcrit":display(crit),
                                   "remaining_B_allowance":display(room_b)},
                   "non_strict_remaining_I_test_also_sufficient_for_saved_strict_lower_target":True,
                   "unpaid_integral_measured":False})
    for path,digest in pins["source_sha256"].items():
        assert sha(path) == digest, path
    assert len(checks) == 4
    output={"status":"PASS_INDEPENDENT_REGION_ACCOUNT_AND_PARTIAL_BUDGET_ONLY",
            "actual_saved_only_configurations":4,"checks":checks,
            "source_hashes_unchanged":True,"kernel_calls":0,"old_budget_evaluations":0,
            "original_059_checker_runs":0,"root_search_calls":0,"global_baseline_replays":0,
            "matlab_calls":0,"subagents":0,"global_I_proved":False,"global_g_proved":False,"new_C6_bound":False,
            "pid":os.getpid(),"python":sys.version,"start_utc":start,
            "end_utc":datetime.now(timezone.utc).isoformat(),"elapsed_seconds":time.monotonic()-tic,
            "hard_process_limit_seconds":180,
            "source_sha256":{name:sha(ROOT/name) for name in ["registration.md","input_manifest.json","proof_review.md","verify_saved.py"]}}
    (ROOT/"results.json").write_text(json.dumps(output,indent=2)+"\n")
    print(json.dumps({"status":output["status"],"configurations":4,"elapsed_seconds":output["elapsed_seconds"],
                      "displays_only":checks[3]["display_only"],"kernel_calls":0},indent=2))


if __name__ == "__main__":
    main()
