#!/usr/bin/env python3
"""Convert the frozen exact fan to primitive integer rays and fixed face points."""
from __future__ import annotations

import hashlib
import json
import math
import random
import argparse
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
FAN = HERE.parent.parent / "round002" / "fan" / "exact_cells_results.json"
DEN = 257
SEED = 6003003


def primitive_ray(ray, rate):
    scale = math.lcm(*(item.denominator for item in ray), rate.denominator)
    values = [int(scale * item) for item in ray]
    q = int(scale * rate)
    divisor = math.gcd(math.gcd(*(abs(v) for v in values)), abs(q))
    return [v // divisor for v in values], q // divisor


def points(denominator, seed):
    rng = random.Random(seed)
    permutations = []
    for _ in range(5):
        p = list(range(32))
        rng.shuffle(p)
        permutations.append(p[:32])
    result = []
    for row in range(32):
        coordinate = [denominator]
        for column in range(5):
            lo = permutations[column][row] * denominator // 32
            hi = (permutations[column][row] + 1) * denominator // 32
            coordinate.append(rng.randrange(max(1, lo + 1), min(denominator, hi)))
        # All coordinates are in (0,1), and require at least three absolute
        # values including 1; this excludes the known S2 structure.
        assert len(set(coordinate)) >= 3
        result.append(coordinate)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--denominator", type=int, default=DEN)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--points-output", default="face_points.txt")
    args = parser.parse_args()
    frozen = json.loads(FAN.read_text())
    with (HERE / "fan_integer_input.txt").open("w") as output:
        output.write(f"{len(frozen['cells'])}\n")
        for cell in frozen["cells"]:
            rays = [tuple(Fraction(x) for x in ray) for ray in cell["rays_a0_equals_one"]]
            rates = [Fraction(x) for x in cell["rates"]]
            converted = [primitive_ray(ray, rate) for ray, rate in zip(rays, rates, strict=True)]
            # det of the independently scaled rays is the numerator needed by
            # alpha=det/prod(q); it is integral for primitive inputs.
            from sympy import Matrix
            determinant = abs(int(Matrix.hstack(*(Matrix(v) for v, _ in converted)).det()))
            output.write(f"{determinant}\n")
            for ray, rate in converted:
                output.write(" ".join(str(value) for value in (rate, *ray)) + "\n")
    fixed_points = points(args.denominator, args.seed)
    with (HERE / args.points_output).open("w") as output:
        output.write(f"{args.denominator} {len(fixed_points)}\n")
        for point in fixed_points:
            output.write(" ".join(map(str, point)) + "\n")
    metadata = {
        "scope": "exact integer fan input and deterministic stratified face points",
        "fan_sha256": hashlib.sha256(FAN.read_bytes()).hexdigest(),
        "denominator": args.denominator,
        "seed": args.seed,
        "points": 32,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (HERE / "input_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
