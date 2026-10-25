#!/usr/bin/env python3
"""Finite exact/Arb audit for the alpha=3/2 stratum and tube interface.

The analytic inequalities are proved in round010_alpha_strata.md.  This script
checks only finite determinant certificates and explicit mathematical constants.
It does not evaluate Green functions or certify target g/I bounds.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
import hashlib
import json
from pathlib import Path
import time

from flint import arb, ctx

ctx.prec = 256
started = time.monotonic()


def det_integer(columns):
    """Fraction-free Bareiss elimination, using Python integers only."""
    a = [list(row) for row in zip(*columns)]
    n, sign, previous = len(a), 1, 1
    for k in range(n-1):
        if a[k][k] == 0:
            pivot = next((j for j in range(k+1, n) if a[j][k]), None)
            if pivot is None:
                return 0
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        current = a[k][k]
        for i in range(k+1, n):
            for j in range(k+1, n):
                numerator = a[i][j]*current-a[i][k]*a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator//previous
            a[i][k] = 0
        previous = current
    return sign*a[-1][-1]


V = [(1,)+p for p in product((-1, 1), repeat=4)]
U = [tuple(-2*int(j in S) for j in range(4))
     for size in range(1, 5) for S in combinations(range(4), size)]
quotient_best = [0]*15
quotient_witness = [None]*15
quotient_histogram = {}
for indexes in combinations(range(15), 4):
    d = abs(det_integer([U[j] for j in indexes]))
    quotient_histogram[d] = quotient_histogram.get(d, 0)+1
    for j in indexes:
        if d > quotient_best[j]:
            quotient_best[j], quotient_witness[j] = d, list(indexes)
assert sum(quotient_histogram.values()) == 1365
assert max(quotient_histogram) == 48
assert all(quotient_best)
reciprocal_sums = [sum((Fraction(1, quotient_best[j])
                         for j, u in enumerate(U) if u[k]), Fraction())
                   for k in range(4)]
assert len(set(reciprocal_sums)) == 1

pair_best = {pair: 0 for pair in combinations(range(16), 2)}
pair_witness = {pair: None for pair in pair_best}
body_histogram = {}
for indexes in combinations(range(16), 5):
    d = abs(det_integer([V[j] for j in indexes]))
    body_histogram[d] = body_histogram.get(d, 0)+1
    for pair in combinations(indexes, 2):
        if d > pair_best[pair]:
            pair_best[pair], pair_witness[pair] = d, list(indexes)
assert sum(body_histogram.values()) == 4368
assert len(pair_best) == 120 and min(pair_best.values()) > 0
assert max(body_histogram) == 48


def ball(x):
    return arb(x.numerator)/arb(x.denominator) if isinstance(x, Fraction) else arb(x)


def textball(x):
    return x.str(55)


pi = arb.pi()
alpha, q = arb(3)/2, arb(10)/7
K = (arb(2)/pi).sqrt()/arb(3)**(arb(1)/4)
Kold = arb(25)/6*(arb(10)/(3*pi)).sqrt()
A = (arb(2)/3).gamma()*arb(2)/(3*pi)
quotient_cauchy_volume = Fraction(14455, 191102976)
P = (ball(quotient_cauchy_volume)*(arb(11)/3).gamma()
     *arb(15)**(arb(4)/3)/(2*pi)**4)
Pdet = A**4/48
Cf = K*arb(16)**(arb(7)/2)*P
B = arb(5)/6*K*(arb(16)**(arb(7)/2)*P)**(arb(3)/10)
B_old = (arb(5)/6*Kold
         *(arb(16)**(arb(7)/2)*Pdet)**(arb(3)/10))
C = K*arb(16)**(arb(5)/2)*A**3*ball(reciprocal_sums[0])
R = arb(5)/2*P**(arb(3)/5)*C**(arb(2)/5)
D = K*arb(16)**(arb(5)/2)*R
corner_norm = (arb(14)/9)**(arb(7)/10)
full_face_I_bound = (D*(15+corner_norm))**q
checks = {"K < 607/1000": K < arb(607)/1000,
          "P < 1/125000": P < arb(1)/125000,
          "B < 27/100": B < arb(27)/100,
          "D < 3": D < arb(3)}
assert all(checks.values()), checks
cutoffs = []
for denominator in (4096, 16384, 65536):
    epsilon = arb(1)/denominator
    tube_norm = D*(corner_norm*epsilon**(arb(9)/5)
                   +15*epsilon**(arb(14)/5))
    cutoffs.append({"epsilon_or_delta": f"1/{denominator}",
                    "short_time_Lq_norm_bound": textball(B*epsilon**(arb(6)/5)),
                    "corner_tube_Lq_norm_bound": textball(tube_norm),
                    "corner_tube_integral_G_to_q_bound": textball(tube_norm**q)})

payload = {
    "model": "alpha=3/2, 16 independent body directions in R^5",
    "scope": "Finite geometry and explicit constants only; no g/I-target or C5 certificate.",
    "runtime": {"arb_precision_bits": ctx.prec, "elapsed_seconds": time.monotonic()-started},
    "quotient_geometry": {
        "number_of_four_sets": 1365,
        "absolute_determinant_histogram": quotient_histogram,
        "directions": U,
        "best_determinants": quotient_best,
        "best_basis_indexes_zero_based": quotient_witness,
        "coordinate_reciprocal_sums": [str(s) for s in reciprocal_sums]},
    "body_geometry": {
        "number_of_five_sets": 4368,
        "absolute_determinant_histogram": body_histogram,
        "number_of_pairs_completed": 120,
        "minimum_best_pair_determinant": min(pair_best.values()),
        "maximum_best_pair_determinant": max(pair_best.values()),
        "pair_witnesses": [{"pair_indexes_zero_based": list(pair),
                            "best_determinant": pair_best[pair],
                            "basis_indexes_zero_based": pair_witness[pair]}
                           for pair in pair_best]},
    "constants": {
        "K_contour": textball(K), "K_round009": textball(Kold),
        "A_one_dimensional_peak": textball(A),
        "P_quotient_peak_Jensen_volume": textball(P),
        "P_quotient_peak_old_det48": textball(Pdet),
        "P_new_over_old": textball(P/Pdet),
        "Cf_face_peak": textball(Cf),
        "B_short_time": textball(B), "B_round009": textball(B_old),
        "B_new_over_old": textball(B/B_old),
        "C_quotient_nonzero_short_time": textball(C),
        "R_weighted_quotient_potential": textball(R),
        "D_full_G_face_envelope": textball(D),
        "whole_positive_face_I_very_coarse_bound": textball(full_face_I_bound)},
    "rational_ball_checks": checks,
    "finite_cutoff_budgets": cutoffs,
    "files_sha256": {
        Path(__file__).name: hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "round010_alpha_strata_registration.json": hashlib.sha256(
            Path(__file__).with_name("round010_alpha_strata_registration.json").read_bytes()).hexdigest(),
        "round002_quotient_peak_proof.md": hashlib.sha256(
            Path(__file__).with_name("round002_quotient_peak_proof.md").read_bytes()).hexdigest()}
}
raw = (json.dumps(payload, sort_keys=True, indent=2)+"\n").encode()
destination = Path(__file__).with_name("round010_alpha_strata_results.json")
destination.write_bytes(raw)
print(json.dumps({"output": str(destination),
                  "sha256": hashlib.sha256(raw).hexdigest(),
                  "constants": payload["constants"],
                  "cutoffs": cutoffs,
                  "checks": checks,
                  "elapsed_seconds": payload["runtime"]["elapsed_seconds"]}, indent=2))
