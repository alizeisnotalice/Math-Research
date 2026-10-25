#!/usr/bin/env python3
"""Read-only source audit; one fresh, subdivided Arb integration at 384 bits.

python check_stable015.py --upper-dir PATH --numerics-dir PATH
Requires python-flint. It does not execute or import the supplied research code.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import platform
import time
from fractions import Fraction
from pathlib import Path
import flint
from flint import acb, arb, ctx


def determinant(matrix):
    """Fraction-free Bareiss elimination, independently of source permutations."""
    a = [list(row) for row in matrix]
    sign, previous = 1, 1
    for k in range(len(a) - 1):
        if a[k][k] == 0:
            r = next(i for i in range(k + 1, len(a)) if a[i][k])
            a[k], a[r] = a[r], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, len(a)):
            for j in range(k + 1, len(a)):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator // previous
        for i in range(k + 1, len(a)):
            a[i][k] = 0
        previous = pivot
    return sign * a[-1][-1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--upper-dir", type=Path, required=True)
    parser.add_argument("--numerics-dir", type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    inputs = {
        "upper_md": args.upper_dir / "round009_alpha_shorttime.md",
        "upper_py": args.upper_dir / "round009_alpha_shorttime.py",
        "upper_json": args.upper_dir / "round009_alpha_shorttime_results.json",
        "kernel_py": args.numerics_dir / "010_arb_stable_density.py",
        "kernel_json": args.numerics_dir / "010_arb_stable_density.json",
    }
    hashes = {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in inputs.items()}
    upper = json.loads(inputs["upper_json"].read_text())
    kernel = json.loads(inputs["kernel_json"].read_text())
    assert hashes["kernel_py"] == kernel["source_sha256"]
    assert kernel["status"] == "complete"
    ctx.prec = 384
    alpha = arb(3) / 2
    c = (alpha + 1).gamma() * (arb.pi() * alpha / 2).sin() / arb.pi()
    tail_constant = 2 * c / alpha + 2 * c / (2 - alpha)
    a = alpha / (alpha + 1)
    K = tail_constant / (2 * (1 - a) * a ** alpha)
    A = (1 / alpha).gamma() / (arb.pi() * alpha)
    Cf = (arb(16) / 48) * K * arb(16) ** (1 + alpha) * A ** 4
    B = arb(5) / 6 * Cf ** (arb(3) / 10) * K ** (arb(7) / 10)
    constants = {"c": c, "tail_constant": tail_constant, "K": K, "A": A, "Cf": Cf, "B": B}
    for name, value in constants.items():
        assert abs(value - arb(upper[name + "_decimal"])) < arb("1e-38")
    assert abs(K - arb(25) / 6 * (arb(10) / (3 * arb.pi())).sqrt()) < arb("1e-110")
    assert Fraction(-5, 3) * Fraction(3, 10) + Fraction(7, 10) == Fraction(1, 5)
    directions = tuple((1,) + e for e in itertools.product((-1, 1), repeat=4))
    v = (1, 1, 1, 1, 1)
    w = ((1, -1, -1, -1, 1), (1, -1, -1, 1, -1),
         (1, -1, 1, -1, -1), (1, 1, -1, -1, -1))
    quotients = tuple(tuple(z[j] - z[0] for j in range(1, 5)) for z in w)
    assert determinant((v,) + w) == determinant(quotients) == -48
    reflection_count = 0
    for chosen in directions:
        transform = lambda z: tuple(z[i] * chosen[i] for i in range(5))
        assert transform(chosen) == v
        assert set(map(transform, directions)) == set(directions)
        reflection_count += 1
    U = arb(4)
    tail = 2 / (3 * arb.pi() * U) * (-U ** 3).exp()
    rows = []
    assert [row["x"] for row in kernel["rows"]] == [0, 1, 2, 8, 32, 128]
    for source in kernel["rows"]:
        point = source["x"]
        x = arb(point)
        # u exp(-u^3) cos(xu^2) is entire, so analytic flag needs no branch check.
        f = lambda u, analytic: u * (-u ** 3).exp() * (x * u ** 2).cos()
        pieces = [acb.integral(f, k, k + 1) for k in range(4)]
        truncated = (2 / acb.pi() * sum(pieces, acb(0))).real
        assert truncated.is_finite()
        certified = truncated + arb(0, tail.upper())
        source_truncated = arb(source["truncated_integral_arb_ball"])
        source_full = arb(source["s1_arb_ball_with_tail"])
        assert truncated in source_truncated
        assert certified in source_full
        rows.append({"x": point, "independent_truncated_ball": str(truncated),
                     "independent_full_ball": str(certified),
                     "source_truncated_contains_new_ball": True,
                     "source_full_contains_new_ball": True})
    assert A in arb(kernel["rows"][0]["s1_arb_ball_with_tail"])
    payload = {
        "status": "passed", "input_sha256": hashes,
        "audit_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "python": platform.python_version(), "python_flint": flint.__version__,
        "precision_bits": ctx.prec, "source_code_executed": False,
        "kernel_independent_runs": 1, "interval_split": [0, 1, 2, 3, 4],
        "constant_balls": {name: str(value) for name, value in constants.items()},
        "exact_full_determinant": -48, "exact_quotient_determinant": -48,
        "body_reflections_verified": reflection_count,
        "time_exponent_exact": "1/5", "integrated_time_exponent_exact": "6/5",
        "tail_ball": str(tail), "kernel_rows": rows,
        "independence_limit": "Fresh script, 384-bit subdivided integration; same python-flint/Arb backend and same analytic u-substitution as source, not an independent interval implementation.",
        "scope": "Auxiliary constants, analytic short-time interface, six finite 1D kernel values only; no global g/I/C5 result.",
        "elapsed_seconds": time.monotonic() - started,
    }
    destination = Path(__file__).with_name("check_results.json")
    destination.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": payload["status"], "elapsed_seconds": payload["elapsed_seconds"],
                      "input_sha256": hashes, "result": str(destination)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
