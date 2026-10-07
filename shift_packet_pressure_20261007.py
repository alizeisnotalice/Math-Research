#!/usr/bin/env python3
"""Fresh finite probes of the proved random-grid packet obstruction.

No FIRST gates. This is a product-uniform hardband proxy, not a test of
the full geom traffic. All means have simultaneous finite-sample bounds;
floating arithmetic error is not interval-certified.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib, json, math, random

HERE = Path(__file__).resolve().parent
OUT = HERE / 'shift_packet_pressure_20261007_results.json'
if OUT.exists():
    raise SystemExit('Refusing to overwrite completed data')
records = []
N = 16384
# Nine configurations, mean and two probability estimates per configuration.
epsilon = math.sqrt(math.log(2 * 27 / .05) / (2 * N))
for ri, n in enumerate((8, 32, 128)):
    for ci, c in enumerate((F(1,4), F(1,2), F(3,4))):
        seed = 202610071307 + 1000 * ri + ci
        rng = random.Random(seed)
        sumz = 0.0
        hits = [0, 0]
        containment_hits = 0
        for _ in range(N):
            z = 1.0
            intact = True
            for j in range(n):
                if rng.random() < float(c):
                    intact = False
                    u = rng.random()
                    z *= max(u, 1-u)
            sumz += z
            hits[0] += z > .5
            hits[1] += z > .25
            containment_hits += intact
        exact_mean = (1-c/4)**n
        exact_containment = (1-c)**n
        mean = sumz/N
        assert abs(mean-float(exact_mean)) <= epsilon
        alpha = float(c)/n
        dstar = alpha + math.sqrt(alpha*alpha+alpha)
        log_opt = 2*n*math.asinh(math.sqrt(alpha))
        direct = n*(math.log1p(dstar)-math.log1p(-alpha/dstar))
        assert abs(log_opt-direct) < 1e-11
        derivative_residual = dstar*dstar-2*alpha*dstar-alpha
        # Exact comparisons at rational scale choices only. Independent of MC.
        scale_rows = []
        for L in (1, math.ceil(math.sqrt(n)), n):
            delta = F(L,n)
            value = ((1+delta)/(1-c/L))**n
            scale_rows.append({'L':L, 'delta':str(delta),
                               'cost_exact':str(value),
                               'log_cost':math.log(value.numerator)-math.log(value.denominator)})
        records.append(dict(round=ri+1,n=n,c=str(c),seed=seed,samples=N,
            mean_exact=str(exact_mean),mean=mean,
            mean_joint_interval=[max(0.,mean-epsilon),min(1.,mean+epsilon)],
            containment_exact=str(exact_containment),containment_hits=containment_hits,
            containment_scope='Diagnostic count only; no extra coverage claim. Zero does not mean probability zero.',
            dominance=[dict(eta=str(eta),hits=h,estimate=h/N,
                            joint_interval=[max(0.,h/N-epsilon),min(1.,h/N+epsilon)],
                            analytic_upper_exact=str(min(F(1),exact_mean/eta)))
                       for eta,h in zip((F(1,2),F(1,4)),hits)],
            delta_star=dstar,log_optimal_containment_cost=log_opt,
            derivative_float_residual=derivative_residual,
            rational_scale_costs=scale_rows))
result = dict(status='complete',source_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    rounds=3,configurations=9,samples_per_config=N,joint_hoeffding_epsilon=epsilon,
    confidence_scope='At least .95 simultaneous over 27 registered means/probabilities in exact sampling model; does not include arithmetic error or continuum geom.',
    records=records)
OUT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'output':str(OUT),'configs':len(records),'epsilon':epsilon,
                  'rows':[{'n':r['n'],'c':r['c'],'mean':r['mean'],
                           'p_half':r['dominance'][0]['estimate'],
                           'log_best_cost':r['log_optimal_containment_cost']} for r in records]},ensure_ascii=False))
