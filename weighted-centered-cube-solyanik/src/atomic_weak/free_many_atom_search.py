#!/usr/bin/env python3
"""Free-position 12/16-atom follow-up to the archived <=8-atom search.

One atom is fixed at +1 by common dilation.  All other locations in [-1,1]
and all positive masses vary.  Unlike the archived search, positions as close
as 1e-12 to zero are admitted, and K>8.  The singular weighted input cost is
included literally.  This remains an exploratory search, not a proof.
"""
import argparse, json, math, os, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True
import numpy as np
from scipy.optimize import differential_evolution
from multiscale_search import block_halo, score
from run_provenance import load_matching_or_quarantine, signature, write_record

ROOT=Path(__file__).resolve().parent
CACHE_ROOT=ROOT/'signed_cache'/'free_many'


def unpack(vector,n):
    return [(float(x),2.0**float(y))
            for x,y in zip(list(vector[:n-1])+[1.0],vector[n-1:])]


def geometric_seed(n):
    # Pair at +/-1 plus symmetric logarithmically spaced atoms of vanishing mass.
    positions=[-1.0]
    needed=n-2
    for j in range(needed):
        k=1+j//2
        positions.append((1 if j%2==0 else -1)*math.exp(-0.45*k))
    # Mass coordinates follow the n-1 free positions and then the +1 anchor.
    masses=[math.log2(4/3)]+[-20.0]*needed+[math.log2(4/3)]
    return np.array(positions+masses)


def job(task):
    tau,n,seed,maxiter,popsize=task
    output=CACHE_ROOT/f'tau{int(tau*100):02d}_N{n:02d}_seed{seed}.json'
    expected = signature(Path(__file__),
                         {"tau": tau, "N": n, "seed": seed,
                          "maxiter": maxiter, "popsize": popsize},
                         {"family": "free_position_atoms", "free_position_bounds": [-1.0, 1.0],
                          "fixed_position": 1.0, "log2_mass_bounds": [-26.0, 5.0],
                          "nearzero_rejection": 1e-12},
                         [ROOT / "multiscale_search.py", ROOT / "run_provenance.py"])
    cached=load_matching_or_quarantine(output, expected)
    if cached is not None: return cached
    started=time.time(); rng=np.random.default_rng(seed)
    bounds=[(-1.0,1.0)]*(n-1)+[(-26.0,5.0)]*n
    dimension=len(bounds)
    population=np.array([[rng.uniform(a,b) for a,b in bounds]
                         for _ in range(popsize*dimension)])
    pair=np.array([-1.0]+[0.35]*(n-2)+[math.log2(4/3)]+[-24.0]*(n-2)+[math.log2(4/3)])
    population[0]=pair
    population[1]=geometric_seed(n)
    calls=0; nearzero=0
    def objective(vector):
        nonlocal calls,nearzero
        calls+=1; atoms=unpack(vector,n)
        if min(abs(x) for x,m in atoms)<1e-12:
            nearzero+=1; return 1e100
        value=score(atoms,tau)
        return -value if math.isfinite(value) else 1e100
    initial=max(-objective(row) for row in population)
    result=differential_evolution(objective,bounds,init=population,seed=seed,
                                  maxiter=maxiter,popsize=popsize,tol=0,atol=0,
                                  polish=False,updating='immediate')
    atoms=unpack(result.x,n)
    record={'status':'completed','tau':tau,'N':n,'seed':seed,'maxiter':maxiter,
            'popsize':popsize,'evaluations':calls,'nearzero_rejections':nearzero,
            'seconds':time.time()-started,'initial_score':initial,
            'final_score':score(atoms,tau),'atoms':[[float(x),float(m)] for x,m in atoms],
            'halo':[[float(a),float(b)] for a,b in block_halo(atoms)],
            'message':str(result.message)}
    record = write_record(output, record, expected)
    return record


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=15)
    parser.add_argument('--maxiter',type=int,default=1200);parser.add_argument('--popsize',type=int,default=14)
    parser.add_argument('--seeds',type=int,default=2);parser.add_argument('--ns',type=int,nargs='+',default=[12,16])
    args=parser.parse_args();CACHE_ROOT.mkdir(parents=True, exist_ok=True)
    tasks=[(tau,n,450000+int(tau*100)*1000+41*n+rep,args.maxiter,args.popsize)
           for tau in (.25,.5) for n in args.ns for rep in range(args.seeds)]
    print(json.dumps({'event':'start','family':'free_many','tasks':len(tasks),'workers':args.workers,
                      'ns':args.ns,'maxiter':args.maxiter,'popsize':args.popsize},sort_keys=True),flush=True)
    rows=[]
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures={executor.submit(job,task):task for task in tasks}
        for future in as_completed(futures):
            row=future.result();rows.append(row)
            print(json.dumps({'event':'done','tau':row['tau'],'N':row['N'],'seed':row['seed'],
                              'initial':row['initial_score'],'final':row['final_score'],
                              'evals':row['evaluations'],'nearzero':row['nearzero_rejections'],
                              'seconds':row['seconds']},sort_keys=True),flush=True)
    rows.sort(key=lambda r:(r['tau'],r['N'],r['seed']))
    (CACHE_ROOT/'summary.json').write_text(json.dumps(rows,indent=2))
    print(json.dumps({'event':'complete','rows':len(rows)},sort_keys=True))

if __name__=='__main__':main()
