#!/usr/bin/env python3
"""Search a non-reflection-symmetric multi-scale comb.

The configuration is -delta_(-1) plus atoms at exp(-jd), j=0,...,K-1.
It is not covered by the previous symmetric-chain experiment, and for K>=8
has more atoms than the archived arbitrary-position DE search.  Scores are
exploratory double-precision atomic lower bounds, never certificates.
"""
import argparse, json, math, os, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"
import numpy as np
from scipy.optimize import differential_evolution
from multiscale_search import block_halo, score

ROOT = Path(__file__).resolve().parent
CHECKPOINTS = ROOT / "asymmetric_checkpoints"


def comb(vector, k):
    d = vector[0]
    masses = [2.0 ** value for value in vector[1:]]
    # The first mass belongs to -1, then positive positions 1,e^-d,... .
    return [(-1.0, masses[0])] + [(math.exp(-d*j), masses[j+1])
                                      for j in range(k)]


def job(task):
    tau, k, seed, maxiter, popsize = task
    name = f"tau{int(100*tau):02d}_K{k:02d}_seed{seed}"
    output = CHECKPOINTS / f"{name}.json"
    if output.exists():
        return json.loads(output.read_text())
    started = time.time()
    rng = np.random.default_rng(seed)
    bounds = [(0.04, 1.0)] + [(-24.0, 5.0)] * (k + 1)
    dimension = len(bounds)
    initial_population = np.array([[rng.uniform(a, b) for a, b in bounds]
                                   for _ in range(popsize*dimension)])
    # Embed the certified two-atom geometry, then a biased decaying comb.
    initial_population[0] = np.array([0.3, math.log2(4/3), math.log2(4/3)]
                                     + [-22.0]*(k-1))
    initial_population[1] = np.array([0.20, math.log2(4/3), math.log2(4/3)]
                                     + [math.log2(0.5)-0.7*j for j in range(k-1)])
    calls = 0
    def objective(vector):
        nonlocal calls
        calls += 1
        value = score(comb(vector, k), tau)
        return -value if math.isfinite(value) else 1e100
    initial = max(-objective(row) for row in initial_population)
    result = differential_evolution(objective, bounds, init=initial_population,
                                    seed=seed, maxiter=maxiter, popsize=popsize,
                                    tol=0.0, atol=0.0, polish=False,
                                    updating="immediate")
    atoms = comb(result.x, k)
    record = {"status": "completed", "tau": tau, "K": k, "seed": seed,
              "maxiter": maxiter, "popsize": popsize, "evaluations": calls,
              "seconds": time.time()-started, "initial_score": initial,
              "final_score": score(atoms,tau), "d": float(result.x[0]),
              "log2_masses": [float(y) for y in result.x[1:]],
              "atoms": [[float(x),float(m)] for x,m in atoms],
              "halo": [[float(a),float(b)] for a,b in block_halo(atoms)],
              "message": str(result.message)}
    temporary=output.with_suffix(".tmp")
    temporary.write_text(json.dumps(record,indent=2,allow_nan=False)); temporary.replace(output)
    return record


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--workers",type=int,default=4)
    parser.add_argument("--maxiter",type=int,default=1000)
    parser.add_argument("--popsize",type=int,default=12)
    parser.add_argument("--seeds",type=int,default=2)
    parser.add_argument("--ks",type=int,nargs="+",default=[8,12,16])
    args=parser.parse_args(); CHECKPOINTS.mkdir(exist_ok=True)
    tasks=[(tau,k,350000+int(100*tau)*1000+37*k+rep,args.maxiter,args.popsize)
           for tau in (0.25,0.5) for k in args.ks for rep in range(args.seeds)]
    print(json.dumps({"event":"start","family":"asymmetric_comb","tasks":len(tasks),
                      "workers":args.workers,"ks":args.ks,"maxiter":args.maxiter,
                      "popsize":args.popsize},sort_keys=True),flush=True)
    rows=[]
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures={executor.submit(job,task):task for task in tasks}
        for future in as_completed(futures):
            row=future.result(); rows.append(row)
            print(json.dumps({"event":"done","tau":row["tau"],"K":row["K"],
                              "seed":row["seed"],"initial":row["initial_score"],
                              "final":row["final_score"],"evals":row["evaluations"],
                              "seconds":row["seconds"]},sort_keys=True),flush=True)
    rows.sort(key=lambda row:(row["tau"],row["K"],row["seed"]))
    (ROOT/"asymmetric_summary.json").write_text(json.dumps(rows,indent=2))
    print(json.dumps({"event":"complete","rows":len(rows)},sort_keys=True))

if __name__=="__main__": main()
