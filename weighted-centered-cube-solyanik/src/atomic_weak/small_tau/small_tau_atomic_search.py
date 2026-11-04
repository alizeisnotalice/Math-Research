#!/usr/bin/env python3
"""Signed exploratory search for small-tau, exponentially separated atoms.

Positions use r_j=exp(-j*d/tau), masses use m_j=q_j*r_j.  Thus the search
allows ratios exp(d/tau), while the weighted input cost of a scale is the
well-conditioned quantity q_j*exp(-j*d).  Decimal stores the physical tiny
positions and masses with arbitrary decimal exponents; no IEEE underflow
represents a candidate as zero.  Results are numerical candidates only.
"""
import argparse
import hashlib
import json
import math
import os
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from decimal import Decimal, getcontext
from pathlib import Path

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True

import numpy as np
import scipy
from scipy.optimize import differential_evolution

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / "signed_checkpoints"
DECIMAL_PRECISION = 90
getcontext().prec = DECIMAL_PRECISION
D = Decimal


def source_hash():
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def signature(task):
    payload = {"schema": "round26-small-tau-v1", "source_sha256": source_hash(),
               "task": task, "python": sys.version, "platform": platform.platform(),
               "numpy": np.__version__, "scipy": scipy.__version__,
               "decimal_precision": DECIMAL_PRECISION,
               "model": {"positions": "exp(-j*d/tau)", "masses": "q_j*r_j",
                         "d_bounds": [0.1, 8.0], "log_q_bounds": [-12.0, 5.0],
                         "threshold": 1.0}}
    text = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return {"sha256": hashlib.sha256(text.encode()).hexdigest(), "payload": payload}


def parse_tau(text):
    numerator, denominator = text.split("/")
    return D(numerator) / D(denominator)


def tau_power(positive, tau):
    assert positive > 0
    return (tau * positive.ln()).exp()


def merge(intervals):
    answer = []
    for left, right in sorted(intervals, key=lambda row: row[0]):
        if left >= right:
            continue
        if answer and left <= answer[-1][1]:
            answer[-1] = (answer[-1][0], max(answer[-1][1], right))
        else:
            answer.append((left, right))
    return answer


def block_halo(atoms):
    """Exact atomic strict halo modulo endpoints, using contiguous blocks."""
    atoms = sorted(atoms, key=lambda row: row[0])
    intervals = []
    for left in range(len(atoms)):
        for right in range(left, len(atoms)):
            # Do not form a global prefix sum: at tau=1/256 a physical
            # micro-scale mass can be 10^-1000 beside a 10^-5 outer mass.
            # Prefix subtraction then loses the micro-scale block entirely.
            mass = sum((atoms[index][1] for index in range(left, right + 1)), D(0))
            a = atoms[right][0] - mass / 2
            b = atoms[left][0] + mass / 2
            if a < b:
                intervals.append((a, b))
    return merge(intervals)


def weighted_mass(intervals, tau):
    tau = D(tau)
    total = D(0)
    for left, right in intervals:
        if right <= 0:
            total += (tau_power(-left, tau) - tau_power(-right, tau)) / tau
        elif left >= 0:
            total += (tau_power(right, tau) - tau_power(left, tau)) / tau
        else:
            total += (tau_power(-left, tau) + tau_power(right, tau)) / tau
    return total


def score(atoms, tau):
    tau = D(tau)
    denominator = sum((mass * tau_power(abs(point), tau) / abs(point)
                       for point, mass in atoms), D(0))
    assert denominator > 0
    return weighted_mass(block_halo(atoms), tau) / denominator


def atoms_from_vector(family, vector, tau, k):
    tau = D(tau)
    d = D(str(vector[0]))
    q = [D(str(z)).exp() for z in vector[1:]]
    radii = [(-D(j) * d / tau).exp() for j in range(k)]
    if family == "symmetric_chain":
        assert len(q) == k
        positive = [(radii[j], q[j] * radii[j]) for j in range(k)]
        return sorted([(-x, mass) for x, mass in positive] + positive)
    if family == "asymmetric_comb":
        assert len(q) == k + 1
        return [(-D(1), q[0])] + [(radii[j], q[j + 1] * radii[j]) for j in range(k)]
    raise ValueError(family)


def model_cost_check(family, vector, tau, k):
    """Independent algebraic form for the positive denominator, used as a control."""
    d = D(str(vector[0])); q = [D(str(z)).exp() for z in vector[1:]]
    if family == "symmetric_chain":
        terms = [q[j] * (-D(j) * d).exp() for j in range(k)]
        return 2 * sum(terms, D(0))
    terms = [q[j + 1] * (-D(j) * d).exp() for j in range(k)]
    return q[0] + sum(terms, D(0))


def verify_cost(family, atoms, vector, tau, k):
    tau = D(tau)
    direct = sum((mass * tau_power(abs(point), tau) / abs(point)
                  for point, mass in atoms), D(0))
    expected = model_cost_check(family, vector, tau, k)
    assert abs(direct - expected) <= D("1e-70") * max(D(1), abs(expected))


def load_or_quarantine(path, expected):
    if not path.exists():
        return None
    try:
        record = json.loads(path.read_text())
    except json.JSONDecodeError:
        record = None
    if record and record.get("run_signature", {}).get("sha256") == expected["sha256"]:
        return record
    stale = path.parent / "stale"; stale.mkdir(exist_ok=True)
    target = stale / (path.stem + ".mismatch.json")
    index = 1
    while target.exists():
        target = stale / (path.stem + f".mismatch-{index}.json"); index += 1
    path.replace(target)
    return None


def job(task):
    family, tau_text, k, seed, maxiter, popsize = task
    tau = parse_tau(tau_text)
    name = f"{family}_tau{tau_text.replace('/','_')}_K{k:02d}_seed{seed}"
    path = CACHE / f"{name}.json"
    expected = signature({"family": family, "tau": tau_text, "K": k, "seed": seed,
                          "maxiter": maxiter, "popsize": popsize})
    cached = load_or_quarantine(path, expected)
    if cached is not None:
        return cached
    dimension = 1 + k + (1 if family == "asymmetric_comb" else 0)
    bounds = [(0.1, 8.0)] + [(-12.0, 5.0)] * (dimension - 1)
    rng = np.random.default_rng(seed)
    population = np.array([[rng.uniform(a, b) for a, b in bounds]
                           for _ in range(popsize * dimension)])
    # Small q values at added scales embed a stable old-scale input.
    if family == "symmetric_chain":
        pair = [0.5, math.log(4 / 3)] + [-10.0] * (k - 1)
    else:
        pair = [0.5, math.log(4 / 3), math.log(4 / 3)] + [-10.0] * (k - 1)
    population[0] = np.array(pair)
    calls = 0
    def objective(vector):
        nonlocal calls
        calls += 1
        atoms = atoms_from_vector(family, vector, tau, k)
        verify_cost(family, atoms, vector, tau, k)
        value = score(atoms, tau)
        return -float(tau * value)
    started = time.time()
    result = differential_evolution(objective, bounds, init=population, seed=seed,
                                    maxiter=maxiter, popsize=popsize, tol=0.0, atol=0.0,
                                    polish=False, updating="immediate")
    atoms = atoms_from_vector(family, result.x, tau, k)
    verify_cost(family, atoms, result.x, tau, k)
    value = score(atoms, tau)
    halo = block_halo(atoms)
    record = {"status": "completed", "family": family, "tau": tau_text, "K": k,
              "atom_count": len(atoms), "seed": seed, "maxiter": maxiter,
              "popsize": popsize, "evaluations": calls, "seconds": time.time()-started,
              "tau_times_score": str(tau * value), "score": str(value),
              "d": str(D(str(result.x[0]))), "log_q": [float(z) for z in result.x[1:]],
              "atoms": [[str(x),str(m)] for x,m in atoms],
              "halo": [[str(a),str(b)] for a,b in halo],
              "message": str(result.message), "run_signature": expected}
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(record, indent=2)); temporary.replace(path)
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--maxiter", type=int, default=300)
    parser.add_argument("--popsize", type=int, default=10)
    parser.add_argument("--seeds", type=int, default=1)
    parser.add_argument("--ks", type=int, nargs="+", default=[3, 5])
    args = parser.parse_args(); CACHE.mkdir(exist_ok=True)
    taus = ["1/64", "1/128", "1/256"]
    tasks=[]
    for tau_index,tau in enumerate(taus):
        for family in ("symmetric_chain", "asymmetric_comb"):
            for k in args.ks:
                for repetition in range(args.seeds):
                    seed=620000+10000*tau_index+1000*(family=="asymmetric_comb")+31*k+repetition
                    tasks.append((family,tau,k,seed,args.maxiter,args.popsize))
    print(json.dumps({"event":"start","tasks":len(tasks),"workers":args.workers,
                      "taus":taus,"ks":args.ks,"maxiter":args.maxiter,"popsize":args.popsize,
                      "decimal_precision":DECIMAL_PRECISION},sort_keys=True),flush=True)
    rows=[]
    with ProcessPoolExecutor(max_workers=min(args.workers,len(tasks))) as executor:
        futures={executor.submit(job,task):task for task in tasks}
        for future in as_completed(futures):
            row=future.result(); rows.append(row)
            print(json.dumps({"event":"done","family":row["family"],"tau":row["tau"],
                              "K":row["K"],"tauR":row["tau_times_score"],
                              "evaluations":row["evaluations"],"seconds":row["seconds"]},sort_keys=True),flush=True)
    rows.sort(key=lambda r:(r['tau'],r['family'],r['K'],r['seed']))
    (ROOT/'small_tau_summary.json').write_text(json.dumps(rows,indent=2))
    print(json.dumps({"event":"complete","records":len(rows),
                      "evaluations":sum(r['evaluations'] for r in rows)},sort_keys=True))


if __name__=="__main__": main()
