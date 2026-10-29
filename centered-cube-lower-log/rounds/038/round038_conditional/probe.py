"""Deterministic rational true-box conditional-moment checks; no Monte Carlo."""
from pathlib import Path
from fractions import Fraction as F
import itertools
import math
import hashlib
import json
import sys
import time
import numpy as np

HERE = Path(__file__).resolve().parent


def quantiles(n, r):
    cumulative = list(itertools.accumulate(math.comb(n, s) for s in range(n+1)))
    upper = [next(s for s, c in enumerate(cumulative) if r*c >= b*2**n)
             for b in range(1, r)] + [n]
    lower = [0] + [v+1 for v in upper[:-1]]
    prior = [sum((F(math.comb(n, s), 2**n) for s in range(a, z+1)), F(0))
             for a, z in zip(lower, upper)]
    assert sum(prior) == 1
    return lower, upper, prior


def box_mass(x, atom, length, delta):
    clip = lambda t: max(F(0), min(F(1), t))
    return sum((clip((x-atom-s+length/2+delta/2)/delta)
                - clip((x-atom-s-length/2+delta/2)/delta)
                for s in [-1, 0, 1]), F(0))


def polynomial(z, lower, upper):
    # z[i][b] is (minus, plus); one factor per coordinate in each monomial.
    total = F(0)
    for b, (lo, hi) in enumerate(zip(lower, upper)):
        coeff = [F(1)]
        for coordinate in z:
            minus, plus = coordinate[b]
            nxt = [F(0)]*(len(coeff)+1)
            for s, c in enumerate(coeff):
                nxt[s] += c*minus
                nxt[s+1] += c*plus
            coeff = nxt
        if lo <= hi:
            total += sum(coeff[lo:hi+1], F(0))
    return total


def run_case(n, r, length):
    started = time.time()
    k = 4
    delta = F(1, 64*k*r)
    lower, upper, prior = quantiles(n, r)
    atoms = [[(F(j)+(F(2*b+1, 2*r) if j<k//2 else 1-F(2*b+1, 2*r)))/k
              for j in range(k)] for b in range(r)]
    points = []
    transition_rows = []
    kernels = []
    for i in range(n):
        b, cell, side = i % r, i % k, (1 if i % 2 == 0 else -1)
        center = (atoms[b][cell]+side*length/2) % 1
        pair = [(center-delta/8) % 1, (center+delta/8) % 1]
        points.append(pair)
        transition_rows.append({"coordinate": i, "table": b, "cell": cell,
                                "side": side, "points": list(map(str, pair)),
                                "bin_width": str(delta/4),
                                "selected_box_masses": [str(box_mass(x, atoms[b][cell], length, delta))
                                                        for x in pair]})
        coordinate = []
        for x in pair:
            coordinate.append([[box_mass(x, a, length, delta)/(k*length)
                                for a in table] for table in atoms])
        kernels.append(coordinate)
    cells = list(itertools.product(range(k), repeat=n))
    labels = [next(b for b, cut in enumerate(upper) if sum(j<k//2 for j in c)<=cut)
              for c in cells]
    histogram = [labels.count(b) for b in range(r)]
    assert [F(h, k**n) for h in histogram] == prior
    values = []
    max_float_relative_error = 0.
    for hidden in itertools.product([0, 1], repeat=n):
        weights = [kernels[i][hidden[i]] for i in range(n)]
        direct = sum((math.prod(weights[i][b][c[i]] for i in range(n))
                      for c, b in zip(cells, labels)), F(0))
        z = [[(sum(w[b][k//2:], F(0)), sum(w[b][:k//2], F(0)))
              for b in range(r)] for w in weights]
        assert direct == polynomial(z, lower, upper)
        # Independent NumPy direct evaluation of the same real-box query.
        wfloat = np.array(weights, dtype=float)
        cellarray = np.asarray(cells)
        labelarray = np.asarray(labels)
        f = np.ones(len(cells))
        for i in range(n):
            f *= wfloat[i, labelarray, cellarray[:, i]]
        approx = float(f.sum())
        max_float_relative_error = max(max_float_relative_error,
            abs(approx-float(direct))/max(1, abs(float(direct))))
        values.append(direct)
    zmean, relative_max = [], []
    kernel_moments = []
    for i in range(n):
        coordinate = []
        rels = []
        rows = []
        for b in range(r):
            bm = []
            for sign, js in [("minus", range(k//2, k)), ("plus", range(k//2))]:
                pair = [sum((kernels[i][h][b][j] for j in js), F(0)) for h in [0, 1]]
                mean = sum(pair)/2
                variance = sum((v-mean)**2 for v in pair)/2
                assert mean > 0
                rel = variance/mean**2
                rels.append(rel)
                bm.append(mean)
                rows.append({"b": b, "class": sign, "mean": str(mean),
                             "variance": str(variance), "relative_variance": str(rel)})
            coordinate.append(bm)
        zmean.append(coordinate)
        relative_max.append(max(rels))
        kernel_moments.append(rows)
    mean = sum(values)/len(values)
    var = sum((v-mean)**2 for v in values)/len(values)
    predicted = polynomial(zmean, lower, upper)
    assert mean == predicted
    cap = math.prod(1+v for v in relative_max)-1
    relative_variance = var/mean**2
    assert relative_variance <= cap
    return {"n": n, "K": k, "R": r, "L": str(length), "delta": str(delta),
            "segments": [{"b": b, "lower": a, "upper": z, "empty": a>z,
                          "prior": str(p)} for b, (a,z,p) in enumerate(zip(lower,upper,prior))],
            "source_mass": "1", "hidden_queries": len(values), "source_cells": len(cells),
            "transitions": transition_rows, "kernel_moments": kernel_moments,
            "conditional_mean": str(mean), "mean_substitution_exact": mean == predicted,
            "conditional_variance": str(var), "relative_variance": str(relative_variance),
            "coordinate_relative_variance_max": list(map(str, relative_max)),
            "relative_variance_cap": str(cap), "cap_minus_actual": str(cap-relative_variance),
            "cap_valid_exact": relative_variance<=cap,
            "relative_variance_float": float(relative_variance), "cap_float": float(cap),
            "max_float_relative_error": max_float_relative_error,
            "seconds": time.time()-started}


def main():
    output = HERE/"RESULTS.json"
    if output.exists():
        raise RuntimeError("Frozen deterministic run exists; refusing replacement.")
    start = time.time()
    cases = []
    for n in [3, 5]:
        for r in [2, 4]:
            for length in [F(10,9), F(7,6), F(11,9)]:
                result = run_case(n, r, length)
                cases.append(result)
                print(json.dumps({k: result[k] for k in ["n","R","L","relative_variance_float",
                                      "cap_float","max_float_relative_error","seconds"]}), flush=True)
    out = {"scope": "exact rational finite true-box conditional moments; no MC",
           "source_SHA256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           "plan_SHA256": hashlib.sha256((HERE/"PLAN.md").read_bytes()).hexdigest(),
           "python": sys.executable, "python_version": sys.version, "numpy": np.__version__,
           "cases": cases, "seconds": time.time()-start}
    output.write_text(json.dumps(out, separators=(",", ":")))
    print("COMPLETE", out["seconds"], flush=True)


if __name__ == "__main__":
    main()
