#!/usr/bin/env python3
"""Recorded third calibration: no new face grid and no QMC batch restart."""
from __future__ import annotations

import hashlib
import itertools
import json
import time
from math import gamma, pi, sin
from pathlib import Path

import numpy as np
from scipy.stats import levy_stable, qmc


ALPHA, SAMPLES = 1.5, 32768
LOG2_48 = np.linspace(-12.0, 10.0, 48)
SUB24 = np.rint(np.linspace(0, 47, 24)).astype(int)
SEEDS = (2026102801, 2026102802)
SOFT_SECONDS = 115.0
GRID_STEP, GRID_MAX = 0.01, 128.0
SOURCE_DIR = Path(__file__).parents[1] / "qmc_face_i_second"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Density:
    def __init__(self):
        self.grid = np.arange(0.0, GRID_MAX + GRID_STEP / 2, GRID_STEP)
        self.values = levy_stable.pdf(self.grid, ALPHA, 0.0)
        self.tail = ALPHA * gamma(ALPHA) * sin(pi * ALPHA / 2) / pi
        self.zero_formula = gamma(1.0 / ALPHA) / (pi * ALPHA)

    def evaluate(self, values):
        absolute = np.abs(values)
        within = absolute <= GRID_MAX
        result = np.empty_like(absolute)
        result[within] = np.interp(absolute[within], self.grid, self.values)
        result[~within] = self.tail * absolute[~within] ** (-1.0 - ALPHA)
        return result, int((~within).sum())


def cms(uniforms):
    angle = pi * (uniforms[:, :11] - 0.5)
    exponential = -np.log(uniforms[:, 11:])
    return (np.sin(ALPHA * angle) / np.cos(angle) ** (1.0 / ALPHA)
            * (np.cos((1.0 - ALPHA) * angle) / exponential) ** ((1.0 - ALPHA) / ALPHA))


def atomic(path, value):
    raw = json.dumps(value, indent=2, sort_keys=True).encode() + b"\n"
    temporary = path.with_suffix(".tmp")
    temporary.write_bytes(raw)
    temporary.replace(path)
    return hashlib.sha256(raw).hexdigest()


def select_points():
    files = [SOURCE_DIR / f"face_batch_{number}_results.json" for number in (1, 2)]
    rows = []
    for filename in files:
        batch = json.loads(filename.read_text())
        for point in batch["points"]:
            rows.append({"source_file": filename.name, "source_sha256": sha(filename), **point})
    values = np.asarray([row["truncated_green_proxy"] for row in rows])
    median = float(np.median(values))
    selections = {
        "old_sample_min": rows[int(np.argmin(values))],
        "old_sample_max": rows[int(np.argmax(values))],
        "old_sample_median": rows[int(np.argmin(np.abs(values - median)))],
    }
    if len({(item["source_file"], item["index"]) for item in selections.values()}) != 3:
        raise AssertionError("selection unexpectedly duplicated")
    return selections, median


def run_one(label, source, seed, density):
    target = Path(__file__).with_name(f"{label}_conditional_seed_{seed}.json")
    if target.exists():
        old = json.loads(target.read_text())
        if old.get("status") != "complete_noncertified" or old.get("source_selection") != source:
            raise RuntimeError(f"refusing to overwrite recorded or partial result: {target}")
        return old, True
    start = time.monotonic()
    directions = {(1,) + signs for signs in itertools.product((-1, 1), repeat=4)}
    base = np.ones((5, 5)) - 2.0 * np.eye(5)
    base[:, 0] *= -1
    inverse = np.linalg.inv(base)
    base_columns = {tuple(int(x) for x in base[:, j]) for j in range(5)}
    other = tuple(sorted(directions - base_columns))
    assert len(other) == 11 and round(abs(np.linalg.det(base))) == 48
    sample = qmc.Sobol(d=22, scramble=True, seed=seed).random_base2(15)
    y = cms(sample) @ np.asarray(other, dtype=float)
    point = np.asarray([1.0] + source["y"], dtype=float)
    times, log_times = 2.0 ** LOG2_48, LOG2_48 * np.log(2.0)
    raw = {"status": "checkpoint_in_progress", "label": label, "seed": seed,
           "source_selection": source, "samples": SAMPLES, "alpha": "3/2",
           "log2_times_48": [float(x) for x in LOG2_48], "subsample_24_indices": [int(x) for x in SUB24],
           "density_values": [], "tail_density_evaluations": 0,
           "density_zero_scipy": float(density.values[0]), "density_zero_gamma_formula": density.zero_formula,
           "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    for node, time_value in enumerate(times):
        coordinates = (point / time_value ** (1.0 / ALPHA) - y) @ inverse.T
        values, tails = density.evaluate(coordinates)
        weights = values.prod(axis=1) * time_value ** (-5.0 / ALPHA) / 48.0
        raw["density_values"].append(float(weights.mean()))
        raw["tail_density_evaluations"] += tails
        raw["elapsed_seconds"] = time.monotonic() - start
        atomic(target, raw)
        if raw["elapsed_seconds"] > SOFT_SECONDS:
            raw["status"] = "checkpoint_soft_limit"
            atomic(target, raw)
            return raw, False
    density_values = np.asarray(raw["density_values"])
    raw["green_48_proxy"] = float(np.trapezoid(times * density_values, log_times))
    raw["green_24_subsample_proxy"] = float(np.trapezoid(times[SUB24] * density_values[SUB24], log_times[SUB24]))
    raw["status"] = "complete_noncertified"
    raw["elapsed_seconds"] = time.monotonic() - start
    atomic(target, raw)
    return raw, False


def main():
    selected, old_median = select_points()
    density = Density()
    all_results, resumed = {}, []
    for label, source in selected.items():
        for seed in SEEDS:
            outcome, was_recorded = run_one(label, source, seed, density)
            all_results.setdefault(label, []).append(outcome)
            if was_recorded:
                resumed.append(f"{label}:{seed}")
            if outcome["status"] != "complete_noncertified":
                print(json.dumps({"status": "stopped_with_checkpoint", "label": label, "seed": seed,
                                  "file": f"{label}_conditional_seed_{seed}.json"}, indent=2))
                return
    overview = {label: {"old_12_proxy": source["truncated_green_proxy"],
                        "green_48_by_block": [item["green_48_proxy"] for item in all_results[label]],
                        "green_24_by_block": [item["green_24_subsample_proxy"] for item in all_results[label]]}
                for label, source in selected.items()}
    print(json.dumps({"status": "all_six_calibrations_complete", "old_median_proxy": old_median,
                      "overview": overview, "already_recorded_skipped": resumed}, indent=2))


if __name__ == "__main__":
    main()
