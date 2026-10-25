#!/usr/bin/env python3
"""Bounded non-certified alpha=3/2 face-I QMC magnitude screen.

Each invocation creates exactly one registered 64-point batch.  A checkpoint
JSON is overwritten only by progressive completion of that same invocation;
the program refuses to restart a pre-existing batch result.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import time
from fractions import Fraction
from math import gamma, pi, sin
from pathlib import Path

import numpy as np
from scipy.stats import levy_stable, qmc


ALPHA = 1.5
Q = 10.0 / 7.0
FACE_POINTS = 64
SAMPLES = 2048
LOG2_TIMES = np.linspace(-8.0, 8.0, 12)
GRID_STEP, GRID_MAX = 0.01, 128.0
SOFT_SECONDS = 115.0
BATCH_PARAMETERS = {
    1: {"face_seed": 2026102601, "conditional_seed": 2026102701},
    2: {"face_seed": 2026102602, "conditional_seed": 2026102702},
}


def check_budget(start):
    if time.monotonic() - start > SOFT_SECONDS:
        raise TimeoutError("120-second soft limit")


class DensityTable:
    def __init__(self):
        self.grid = np.arange(0.0, GRID_MAX + GRID_STEP / 2, GRID_STEP)
        self.pdf = levy_stable.pdf(self.grid, ALPHA, 0.0)
        self.tail = ALPHA * gamma(ALPHA) * sin(pi * ALPHA / 2.0) / pi
        self.zero_formula = gamma(1.0 / ALPHA) / (pi * ALPHA)

    def __call__(self, values):
        x = np.abs(values)
        inside = x <= GRID_MAX
        output = np.empty_like(x)
        output[inside] = np.interp(x[inside], self.grid, self.pdf)
        output[~inside] = self.tail * x[~inside] ** (-1.0 - ALPHA)
        return output, int((~inside).sum())


def cms(uniforms):
    angle = pi * (uniforms[:, :11] - 0.5)
    exponential = -np.log(uniforms[:, 11:])
    return (np.sin(ALPHA * angle) / np.cos(angle) ** (1.0 / ALPHA)
            * (np.cos((1.0 - ALPHA) * angle) / exponential) ** ((1.0 - ALPHA) / ALPHA))


def atomic_write(path, payload):
    raw = json.dumps(payload, indent=2, sort_keys=True).encode() + b"\n"
    temporary = path.with_suffix(".tmp")
    temporary.write_bytes(raw)
    temporary.replace(path)
    return hashlib.sha256(raw).hexdigest()


def aggregate(rows):
    if not rows:
        return None
    green = np.asarray([row["truncated_green_proxy"] for row in rows])
    green_blocks = np.asarray([row["four_conditional_block_green_proxies"] for row in rows])
    i_proxy = float(np.mean(green ** Q))
    g_sample = float(np.min(green))
    cost = float((10.0 / 3.0) ** Q * i_proxy / g_sample ** Q)
    block_i = np.mean(green_blocks ** Q, axis=0)
    block_g = np.min(green_blocks, axis=0)
    block_cost = (10.0 / 3.0) ** Q * block_i / block_g ** Q
    return {
        "point_count": len(rows), "I_proxy": i_proxy, "sampled_minimum_G_proxy": g_sample,
        "conditional_cost_proxy": cost,
        "four_conditional_block_I_proxies": [float(x) for x in block_i],
        "four_conditional_block_minimum_G_proxies": [float(x) for x in block_g],
        "four_conditional_block_cost_proxies": [float(x) for x in block_cost],
        "conditional_block_cost_std": float(np.std(block_cost, ddof=1)),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch", type=int, choices=(1, 2), required=True)
    args = parser.parse_args()
    config = BATCH_PARAMETERS[args.batch]
    output_path = Path(__file__).with_name(f"face_batch_{args.batch}_results.json")
    if output_path.exists():
        raise FileExistsError(f"refusing to restart existing recorded batch: {output_path}")
    started = time.monotonic()
    # The base columns are body directions after the first-column sign correction.
    base = np.ones((5, 5)) - 2.0 * np.eye(5)
    base[:, 0] *= -1.0
    inverse = np.linalg.inv(base)
    base_columns = {tuple(int(x) for x in base[:, j]) for j in range(5)}
    directions = {(1,) + signs for signs in itertools.product((-1, 1), repeat=4)}
    other = tuple(sorted(directions - base_columns))
    assert round(abs(np.linalg.det(base))) == 48 and len(other) == 11
    density = DensityTable()
    face_engine = qmc.Sobol(d=4, scramble=True, seed=config["face_seed"])
    faces = face_engine.random_base2(6)  # exactly 64 new points
    conditional_engine = qmc.Sobol(d=22, scramble=True, seed=config["conditional_seed"])
    stable = cms(conditional_engine.random_base2(11))
    y = stable @ np.asarray(other, dtype=float)
    times, log_times = 2.0 ** LOG2_TIMES, np.log(2.0 ** LOG2_TIMES)
    payload = {
        "status": "in_progress_noncertified_alpha_3_over_2_face_I_screen",
        "batch": args.batch, "seeds": config, "alpha": "3/2", "q": "10/7",
        "face_points_requested": FACE_POINTS, "samples_per_point": SAMPLES,
        "time_nodes_log2": [float(x) for x in LOG2_TIMES],
        "formula": "p_t=t^(-5/alpha)/48 E product s_alpha(B^-1(t^(-1/alpha)x-Y))",
        "base": "J_5-2I_5, first column sign flipped, abs(det)=48",
        "points": [], "tail_density_evaluations": 0,
        "density": {"grid_step": GRID_STEP, "grid_max_abs": GRID_MAX,
                    "zero_gamma_formula": density.zero_formula, "zero_scipy": float(density.pdf[0]),
                    "outside_grid": "leading stable density tail; noncertified"},
        "limitations": ["finite log-time range", "SciPy interpolation", "leading-tail substitution",
                        "scrambled Sobol and block variation are noncertified", "sampled min is not global g"],
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    stopped = False
    for index, face in enumerate(faces, 1):
        try:
            densities, block_densities = [], []
            point = np.concatenate(([1.0], face))
            for time_value in times:
                coordinates = (point / time_value ** (1.0 / ALPHA) - y) @ inverse.T
                factors, tails = density(coordinates)
                payload["tail_density_evaluations"] += tails
                weights = factors.prod(axis=1) * time_value ** (-5.0 / ALPHA) / 48.0
                densities.append(float(weights.mean()))
                block_densities.append([float(item.mean()) for item in np.array_split(weights, 4)])
            green = float(np.trapezoid(times * np.asarray(densities), log_times))
            by_block = [float(np.trapezoid(times * np.asarray(item), log_times))
                        for item in np.asarray(block_densities).T]
            payload["points"].append({"index": index, "y": [float(x) for x in face],
                                      "truncated_green_proxy": green,
                                      "four_conditional_block_green_proxies": by_block})
            payload["aggregate"] = aggregate(payload["points"])
            payload["elapsed_seconds"] = time.monotonic() - started
            atomic_write(output_path, payload)
            check_budget(started)
        except TimeoutError:
            stopped = True
            break
    payload["status"] = "checkpoint_soft_limit" if stopped else "complete_noncertified"
    payload["elapsed_seconds"] = time.monotonic() - started
    payload["aggregate"] = aggregate(payload["points"])
    result_hash = atomic_write(output_path, payload)
    print(json.dumps({"status": payload["status"], "points_completed": len(payload["points"]),
                      "aggregate": payload["aggregate"], "elapsed_seconds": payload["elapsed_seconds"]}, indent=2))
    print("sha256", result_hash)


if __name__ == "__main__":
    main()
