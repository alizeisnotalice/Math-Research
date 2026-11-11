#!/usr/bin/env python3
"""Stream completed r093 MC/Gaussian rows into a small auditable summary.

This script reads only the explicitly named /tmp outputs and r093 source files;
it does not rerun simulations and never copies the raw JSONL into the repository.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

ROOT = Path("/Volumes/data/CodexProjects/T")
R093 = ROOT / "work/r093"
OUT = R093 / "extra93/summary93.json"
WEIGHT_RUNS = [
    {
        "name": "strong_weight_probe93",
        "rows": Path("/tmp/strong_weight_probe93/rows.jsonl"),
        "summary": Path("/tmp/strong_weight_probe93/summary.json"),
        "resources": Path("/tmp/strong_weight_probe93_resources.json"),
        "source": Path("/tmp/strong_weight_probe93.py"),
    },
    {
        "name": "strong_weight_probe93b",
        "rows": Path("/tmp/strong_weight_probe93b/rows.jsonl"),
        "summary": Path("/tmp/strong_weight_probe93b/summary.json"),
        "resources": Path("/tmp/strong_weight_probe93b_resources.json"),
        "source": Path("/tmp/strong_weight_probe93b.py"),
    },
]
GAUSSIAN = {
    "name": "strong_fullcpu93",
    "rows": Path("/tmp/strong_fullcpu93/rows.jsonl"),
    "summary": Path("/tmp/strong_fullcpu93/summary.json"),
    "resources": Path("/tmp/strong_fullcpu93_resources.json"),
    "source": Path("/tmp/strong_fullcpu93.py"),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path: Path):
    return json.loads(path.read_text())


def read_jsonl(path: Path):
    with path.open() as f:
        for raw in f:
            if raw.strip():
                yield json.loads(raw), hashlib.sha256(raw.encode()).hexdigest()


def mu_root() -> float:
    lo, hi = 1.0, 2.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if mid - 1 - math.exp(-mid) < 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


MU = mu_root()
SIGMA = math.sqrt(MU)


def critical_target(n: int):
    root = math.sqrt(n)
    lo, hi = 0.0, max(2.0, math.sqrt(math.log(n)) + 1.0)
    for _ in range(100):
        mid = (lo + hi) / 2
        q = math.erfc(mid / math.sqrt(2.0)) / 2.0
        rhs = SIGMA * mid / root
        if q > rhs:
            lo = mid
        else:
            hi = mid
    z = (lo + hi) / 2
    tau = SIGMA * z / root
    b = math.log1p(tau) - math.log1p(-tau)
    return {"z": z, "tau": tau, "b": b}


def main():
    if any(not p.exists() for run in WEIGHT_RUNS for p in run.values() if isinstance(p, Path)):
        missing = [str(p) for run in WEIGHT_RUNS for p in run.values()
                   if isinstance(p, Path) and not p.exists()]
        raise FileNotFoundError(missing)
    if any(not p.exists() for p in GAUSSIAN.values() if isinstance(p, Path)):
        missing = [str(p) for p in GAUSSIAN.values()
                   if isinstance(p, Path) and not p.exists()]
        raise FileNotFoundError(missing)

    critical_rows = []
    batch_summaries = []
    source_hashes = {}
    input_hashes = {}
    per_n = {}
    all_rows_count = 0
    seen_specs = set()
    seen_seeds = set()
    total_draw_pairs = 0
    total_values_across_both_laws = 0
    failure_union = 0.0
    weighted_elapsed_s = 0.0
    resources_list = []
    threshold_stats = {
        name: {"count": 0, "min_mean": None, "min_lower": None, "min_upper": None}
        for name in ("q_0.02", "q_0.05", "q_0.1", "critical")
    }

    for run in WEIGHT_RUNS:
        summary = read_json(run["summary"])
        resources = read_json(run["resources"])
        if not summary.get("pass_batch") or summary.get("actual_workers") != 12:
            raise ValueError(f"non-passing or non-12-worker batch: {run['name']}")
        if resources.get("exit_code") != 0 or resources.get("stop_reason") is not None:
            raise ValueError(f"resource report not clean: {run['name']}")
        if resources.get("reserve_GiB", 0) < 1:
            raise ValueError(f"reserve below 1 GiB: {run['name']}")
        batch_summaries.append({"name": run["name"], **summary})
        source_hashes[str(run["source"])] = sha256(run["source"])
        for key in ("rows", "summary", "resources"):
            input_hashes[str(run[key])] = sha256(run[key])
        resources_list.append({"name": run["name"], **resources})
        weighted_elapsed_s += float(summary["seconds"])

        run_rows = 0
        for row, line_hash in read_jsonl(run["rows"]):
            run_rows += 1
            all_rows_count += 1
            n = int(row["N"])
            m = int(row["m"])
            w = float(row["w"])
            c = float(row["c"])
            rep = int(row["replicate"])
            spec = (n, m, w, c, rep)
            if spec in seen_specs:
                raise ValueError(f"duplicate parameter/replicate job across batches: {spec}")
            seen_specs.add(spec)
            expected_seed = int.from_bytes(hashlib.sha256(repr(spec).encode()).digest()[:8], "big")
            if int(row["seed"]) != expected_seed:
                raise ValueError(f"seed does not match source derivation for {spec}")
            if int(row["seed"]) in seen_seeds:
                raise ValueError("seed collision between distinct Monte Carlo jobs")
            seen_seeds.add(int(row["seed"]))
            deficit = (1.0 - w) * m
            if deficit <= 0:
                raise ValueError("weight deficit must be positive")
            if float(row["minimum_samplewise_gap"]) < -1e-10:
                raise ValueError("common-random-number dominance check failed")
            target = critical_target(n)
            expected_eta = math.sqrt(math.log(4e8) / (2 * int(row["samples"]))) + 2 / int(row["samples"])
            if abs(float(row["CDF_error_with_interpolation_allowance"]) - expected_eta) > 1e-15:
                raise ValueError("unexpected DKW/interpolation allowance")
            if abs(float(row["per_job_failure_bound"]) - 1e-8) > 1e-20:
                raise ValueError("unexpected per-job DKW union bound")
            near = [r for r in row["rows"]
                    if abs(float(r["CDF_probability"]) - target["tau"]) <= 5e-12]
            if len(near) != 1:
                raise ValueError(f"critical probability not uniquely present: N={n}, found={len(near)}")
            if len(row["rows"]) != 4:
                raise ValueError(f"expected 4 distinct threshold rows, found {len(row['rows'])}")
            qrow = near[0]
            raw_low, raw_high = map(float, qrow["DKW_gap_interval"])
            low, high = raw_low / deficit, raw_high / deficit
            value = float(qrow["gap_per_weight_deficit"])
            record = {
                "run": run["name"], "N": n, "m": m, "w": w, "c": c,
                "replicate": rep, "seed": int(row["seed"]),
                "critical_CDF_probability": float(qrow["CDF_probability"]),
                "b": target["b"], "mean_gap_per_weight_deficit": value,
                "normalized_DKW_lower": low, "normalized_DKW_upper": high,
                "upper_minus_b": high - target["b"], "row_sha256": line_hash,
            }
            critical_rows.append(record)
            per_n.setdefault(str(n), []).append(record)
            for threshold_row in row["rows"]:
                q = float(threshold_row["CDF_probability"])
                if abs(q - target["tau"]) <= 5e-12:
                    category = "critical"
                else:
                    fixed = min((0.02, 0.05, 0.1), key=lambda x: abs(q - x))
                    if abs(q - fixed) > 1e-14:
                        raise ValueError(f"unrecognized tail probability {q}")
                    category = f"q_{fixed:g}"
                raw_lo, raw_hi = map(float, threshold_row["DKW_gap_interval"])
                tr = {
                    "run": run["name"], "N": n, "m": m, "w": w, "c": c,
                    "replicate": rep, "seed": int(row["seed"]),
                    "CDF_probability": q,
                    "mean_gap_per_weight_deficit": float(threshold_row["gap_per_weight_deficit"]),
                    "normalized_DKW_lower": raw_lo / deficit,
                    "normalized_DKW_upper": raw_hi / deficit,
                    "row_sha256": line_hash,
                }
                st = threshold_stats[category]
                st["count"] += 1
                if st["min_mean"] is None or tr["mean_gap_per_weight_deficit"] < st["min_mean"]["mean_gap_per_weight_deficit"]:
                    st["min_mean"] = tr
                if st["min_lower"] is None or tr["normalized_DKW_lower"] < st["min_lower"]["normalized_DKW_lower"]:
                    st["min_lower"] = tr
                if st["min_upper"] is None or tr["normalized_DKW_upper"] < st["min_upper"]["normalized_DKW_upper"]:
                    st["min_upper"] = tr
            total_draw_pairs += int(row["samples"])
            total_values_across_both_laws += 2 * int(row["samples"])
            failure_union += float(row["per_job_failure_bound"])
        if run_rows != int(summary["jobs"]):
            raise ValueError(f"row count mismatch for {run['name']}: {run_rows} != {summary['jobs']}")

    gaussian_summary = read_json(GAUSSIAN["summary"])
    gaussian_resources = read_json(GAUSSIAN["resources"])
    if not gaussian_summary.get("pass_batch") or gaussian_summary.get("actual_workers") != 12:
        raise ValueError("Gaussian batch did not pass with 12 workers")
    if gaussian_resources.get("exit_code") != 0 or gaussian_resources.get("stop_reason") is not None:
        raise ValueError("Gaussian resource report not clean")
    if gaussian_resources.get("reserve_GiB", 0) < 1:
        raise ValueError("Gaussian reserve below 1 GiB")
    source_hashes[str(GAUSSIAN["source"])] = sha256(GAUSSIAN["source"])
    for key in ("rows", "summary", "resources"):
        input_hashes[str(GAUSSIAN[key])] = sha256(GAUSSIAN[key])
    resources_list.append({"name": GAUSSIAN["name"], **gaussian_resources})

    gaussian_cases = 0
    gaussian_k = set()
    gaussian_bs = set()
    rank_over_b = []
    max_precision_rel_diff = 0.0
    max_root_residual_ratio = 0.0
    precision_pairs_consistent = 0
    for row, _line_hash in read_jsonl(GAUSSIAN["rows"]):
        gaussian_cases += 1
        gaussian_k.add(int(row["k"]))
        gaussian_bs.add(str(row["B"]))
        prs = sorted(row["precision_rows"], key=lambda x: int(x["working_precision"]))
        if len(prs) != 2 or [int(x["working_precision"]) for x in prs] != [80, 160]:
            raise ValueError("expected 80- and 160-digit rows")
        lo, hi = prs
        for pr in prs:
            target = float(pr["root_log_residual_target"])
            for key in ("z_log_residual", "zB_log_residual"):
                ratio = float(pr[key]) / target if target else math.inf
                max_root_residual_ratio = max(max_root_residual_ratio, ratio)
                if ratio > 1:
                    raise ValueError(f"Gaussian root residual exceeds target: {key}")
        for key in ("rank_shift", "rank_shift_over_B", "z", "zB"):
            a, b = float(lo[key]), float(hi[key])
            denom = max(abs(a), abs(b))
            rel = abs(a - b) / denom if denom else 0.0
            max_precision_rel_diff = max(max_precision_rel_diff, rel)
        if float(lo["rank_shift"]) <= 0 or float(hi["rank_shift"]) <= 0:
            raise ValueError("expected positive Gaussian rank shift")
        if float(lo["rank_shift_over_B"]) <= 0 or float(hi["rank_shift_over_B"]) <= 0:
            raise ValueError("expected positive normalized rank shift")
        precision_pairs_consistent += 1
        rank_over_b.extend((float(lo["rank_shift_over_B"]), float(hi["rank_shift_over_B"])))
    if gaussian_cases != int(gaussian_summary["cases"]):
        raise ValueError("Gaussian row count mismatch")

    critical_rows.sort(key=lambda r: r["upper_minus_b"])
    standout = critical_rows[:10]
    all_upper_below = [r for r in critical_rows if r["normalized_DKW_upper"] < r["b"]]
    min_mean = min(critical_rows, key=lambda r: r["mean_gap_per_weight_deficit"])
    min_lower = min(critical_rows, key=lambda r: r["normalized_DKW_lower"])
    min_upper = min(critical_rows, key=lambda r: r["normalized_DKW_upper"])
    global_stats = {
        "critical_rows": len(critical_rows),
        "case_replicates": all_rows_count,
        "distinct_parameter_replicate_specs": len(seen_specs),
        "distinct_seeds": len(seen_seeds),
        "duplicate_spec_count": all_rows_count - len(seen_specs),
        "per_job_sample_draw_pairs": 2_000_000,
        "total_sample_draw_pairs": total_draw_pairs,
        "total_values_in_two_coupled_laws": total_values_across_both_laws,
        "per_job_two_law_DKW_failure_bound": 1e-8,
        "family_union_bound_over_weighted_jobs": failure_union,
        "minimum_observed_mean": min_mean,
        "minimum_normalized_DKW_lower_endpoint": min_lower,
        "minimum_normalized_DKW_upper_endpoint": min_upper,
        "count_normalized_DKW_upper_below_b": len(all_upper_below),
        "Dkw_upper_below_b_by_N": {
            n: sum(r["normalized_DKW_upper"] < r["b"] for r in rows)
            for n, rows in sorted(per_n.items(), key=lambda kv: int(kv[0]))
        },
        "critical_rows_by_N": {
            n: {
                "count": len(rows),
                "min_mean": min(r["mean_gap_per_weight_deficit"] for r in rows),
                "min_normalized_lower": min(r["normalized_DKW_lower"] for r in rows),
                "min_normalized_upper": min(r["normalized_DKW_upper"] for r in rows),
                "min_b": min(r["b"] for r in rows),
                "max_b": max(r["b"] for r in rows),
                "count_normalized_DKW_upper_below_b": sum(r["normalized_DKW_upper"] < r["b"] for r in rows),
            }
            for n, rows in sorted(per_n.items(), key=lambda kv: int(kv[0]))
        },
        "ten_closest_DKW_upper_to_b_rows_by_smallest_upper_minus_b": standout,
    }

    # The exact side-vector increment law for c=-log(r): Z_r=E+c B,
    # E~Exp(1), B~Bernoulli(r/(1+r)) independent. For p=1/(1+exp(c)),
    # cB has the same law. The probe uses r_i=e^-c for every coordinate and
    # two weight levels: 1 on N-m entries and w on m entries.
    model_review = {
        "corresponds_to_fixed_side_vector_model": True,
        "law_identity": "Z_r = Exp(1) + c*Bernoulli(p), c=-log(r), p=r/(1+r)=1/(1+exp(c)); summands independent",
        "coupling_identity": "full=base+part and weighted=base+w*part, so full-weighted=(1-w)*part>=0 samplewise",
        "sampled_subfamily": "common side ratio r=e^-c for all coordinates; weights 1 on N-m coordinates and w on m coordinates; not a search over arbitrary heterogeneous (r_i,w_i)",
        "finite_grids": {
            "N": [4, 8, 16, 32, 64, 128],
            "m": "unique values among 1,N//4,N//2",
            "w": [0.25, 0.5, 0.75, 0.9, 0.95, 0.98, 0.995],
            "c": "j/8 for j=1..48, plus mu",
            "replicates": 4,
            "note": "two replicates omit selected customary w,c combinations according to source code; the wider batch contains targeted repeated cases",
        },
        "critical_probability": "tau_N=sigma*z/sqrt(N), Q(z)=sigma*z/sqrt(N); b=log((1+tau_N)/(1-tau_N))",
    }

    result = {
        "scope": "Finite Monte Carlo diagnostics for a two-level weighted, common-side-ratio subfamily of the fixed-side-vector model; not a proof, uniform derivative bound, physical halo certificate, or evidence of a finite onset dimension.",
        "model_review": model_review,
        "critical_quantile_gap_summary": global_stats,
        "all_tail_threshold_sensitivity_summary": {
            "metric_name": "empirical quantile-gap finite-difference per unit total weight deficit; descriptive sensitivity ratio, not a proven derivative bound",
            "categories": threshold_stats,
            "total_threshold_rows": sum(x["count"] for x in threshold_stats.values()),
            "minimum_over_all_threshold_categories": {
                "mean": min((x["min_mean"] for x in threshold_stats.values()), key=lambda r: r["mean_gap_per_weight_deficit"]),
                "normalized_DKW_lower": min((x["min_lower"] for x in threshold_stats.values()), key=lambda r: r["normalized_DKW_lower"]),
                "normalized_DKW_upper": min((x["min_upper"] for x in threshold_stats.values()), key=lambda r: r["normalized_DKW_upper"]),
            },
        },
        "finite_probability_review": {
            "status": "valid finite-family probability statement, conditional on iid sampling within each job and the stated DKW inequality",
            "per_law_DKW_tail": "2*exp(-2*SAMPLES*eta^2)",
            "CDF_error_with_interpolation_allowance": math.sqrt(math.log(4e8) / (2 * 2_000_000)) + 2 / 2_000_000,
            "two_law_union": "per_job_failure_bound=1e-8; common-random-number correlation between laws does not invalidate the union bound",
            "uniformity": "DKW is uniform over the CDF, so all four CDF probabilities in a job are covered simultaneously",
            "quantile_interpolation": "eta includes the source's 2/SAMPLES allowance for NumPy linear-quantile rank interpolation; endpoint probabilities are clamped to [0,1]",
            "across_jobs": "simple union bound over 680+22644 weighted jobs; no independence across jobs is assumed",
            "family_failure_bound_over_weighted_jobs": failure_union,
            "scope_limit": "These are finite-sample simultaneous intervals for the selected simulations only, not a confidence statement about all parameter values or all side/weight vectors.",
        },
        "gaussian_diagnostics": {
            "cases": gaussian_cases,
            "actual_workers": gaussian_summary["actual_workers"],
            "k_min": min(gaussian_k), "k_max": max(gaussian_k), "k_count": len(gaussian_k),
            "sampled_B_values": sorted(gaussian_bs, key=lambda x: float(x.split("/")[0]) / float(x.split("/")[1]) if "/" in x else float(x)),
            "sampled_B_min": min(float(x.split("/")[0]) / float(x.split("/")[1]) if "/" in x else float(x) for x in gaussian_bs),
            "sampled_B_max": max(float(x.split("/")[0]) / float(x.split("/")[1]) if "/" in x else float(x) for x in gaussian_bs),
            "rank_shift_over_B_min": min(rank_over_b),
            "rank_shift_over_B_max": max(rank_over_b),
            "precision_pairs_80_160_checked": precision_pairs_consistent,
            "max_relative_difference_80_vs_160_over_z_zB_rank_shift_and_rank_shift_over_B": max_precision_rel_diff,
            "max_root_residual_over_precision_target": max_root_residual_ratio,
            "scope": "Gaussian rank-window/root diagnostic only; no physical N0 implication.",
        },
        "batches": batch_summaries + [{"name": GAUSSIAN["name"], **gaussian_summary}],
        "resources": resources_list,
        "resource_interpretation": "Guard resource reports are sampled macOS memory-availability/RSS estimates, not hard OS reservations. Reported seconds are batch elapsed wall time; per-process CPU time was not recorded.",
        "total_batch_elapsed_wall_seconds": sum(float(x["seconds"]) for x in batch_summaries + [{"seconds": gaussian_summary["seconds"]}]),
        "source_sha256": {
            **source_hashes,
            str(R093 / "core93.md"): sha256(R093 / "core93.md"),
            str(R093 / "audit93.md"): sha256(R093 / "audit93.md"),
            str(R093 / "experiment93.py"): sha256(R093 / "experiment93.py"),
            str(ROOT / "Math-Research-strong-maximal/strong-maximal-endpoint/code/experiment92.py"): sha256(ROOT / "Math-Research-strong-maximal/strong-maximal-endpoint/code/experiment92.py"),
            str(Path(__file__).resolve()): sha256(Path(__file__).resolve()),
        },
        "input_sha256": input_hashes,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "output": str(OUT),
        "critical_rows": len(critical_rows),
        "min_mean": min_mean,
        "min_normalized_DKW_lower": min_lower,
        "min_normalized_DKW_upper": min_upper,
        "count_DKW_upper_below_b": len(all_upper_below),
        "family_union_bound": failure_union,
        "gaussian_cases": gaussian_cases,
        "rank_shift_over_B_min": min(rank_over_b),
        "rank_shift_over_B_max": max(rank_over_b),
    }, indent=2))


if __name__ == "__main__":
    main()
