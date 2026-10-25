#!/usr/bin/env python3
"""Freeze four scrambled-Sobol, rational-quantized face samples."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.stats import qmc

HERE = Path(__file__).resolve().parent
D = 2**30
SEEDS = (6005001, 6005002, 6005003, 6005004)
PER_REPLICATE = 128


def write_points(path: Path, rows: list[tuple[int, ...]]) -> None:
    path.write_text(f"{D} {len(rows)}\n" + "\n".join(" ".join(map(str, r)) for r in rows) + "\n")


def main() -> None:
    rows: list[tuple[int, ...]] = []
    labels=[]
    for replicate, seed in enumerate(SEEDS):
        values = qmc.Sobol(d=5, scramble=True, seed=seed).random_base2(m=7)
        integers = np.floor(values * D).astype(np.int64)
        assert np.all(integers > 0) and np.all(integers < D)
        for index, row in enumerate(integers):
            rows.append((D, *(int(v) for v in row)))
            labels.append({"global_index":len(rows)-1,"replicate":replicate,"index_in_replicate":index})
    full=HERE/"sobol_4x128.txt"; write_points(full,rows)
    # Interleave the first two points of every replicate for the initial sign check.
    pilot_indices=[0,1,128,129,256,257,384,385]
    pilot=HERE/"sobol_pilot8.txt"; write_points(pilot,[rows[i] for i in pilot_indices])
    metadata={
      "status":"frozen_scrambled_sobol_rational_face_samples",
      "dimension":5,"face":"x1=1, y in (0,1)^5","denominator":D,
      "scipy_version":scipy.__version__,"sobol_scramble":True,"seeds":list(SEEDS),
      "points_per_replicate":PER_REPLICATE,"total_points":len(rows),"pilot_original_indices":pilot_indices,
      "full_sha256":hashlib.sha256(full.read_bytes()).hexdigest(),"pilot_sha256":hashlib.sha256(pilot.read_bytes()).hexdigest(),
      "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"labels":labels,
    }
    (HERE/"input_metadata.json").write_text(json.dumps(metadata,indent=2)+"\n")
    print(json.dumps({k:v for k,v in metadata.items() if k != "labels"},indent=2))


if __name__ == "__main__": main()
