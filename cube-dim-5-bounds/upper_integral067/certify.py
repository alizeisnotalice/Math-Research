#!/usr/bin/env python3
"""One registered macrobox: positive time bands, reused zero-order stencil."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import isqrt, prod
import argparse, gzip, hashlib, json, os, resource, sys, time
import flint
from flint import arb, ctx, fmpq

HERE = Path(__file__).resolve().parent
if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def atomic(p, value):
    t = p.with_suffix('.tmp')
    t.write_text(json.dumps(value, sort_keys=True, separators=(',', ':')) + '\n')
    os.replace(t, p)


def rq(x):
    x = F(str(x))
    return fmpq(x.numerator, x.denominator)


def ep(x, upper=False):
    return str((x.upper() if upper else x.lower()).fmpq())


def enc(x):
    return {'lower': ep(x), 'upper': ep(x, True), 'ball': str(x)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    began = time.monotonic()
    reg = json.loads((HERE/'registration.json').read_text())
    assert sha(HERE/'mechanism.md') == reg['mechanism_sha256']
    for v in reg['sources'].values():
        assert sha(v['path']) == v['sha256']
    ctx.prec = reg['precision_bits']
    assert ctx.prec == 384
    fingerprint = hashlib.sha256((sha(__file__) + sha(HERE/'registration.json')).encode()).hexdigest()
    src = {k: Path(v['path']) for k, v in reg['sources'].items()}
    low = json.loads(src['low063'].read_text())
    cache = json.loads(src['cache063'].read_text())
    assert cache['next_simplex'] == 3504 and cache['center'] == 'low'
    assert cache['fingerprint'] == low['fingerprint']
    c = list(map(F, low['selected_center']['point']))
    D = low['selected_center']['D']
    eps, h = F(reg['eps']), F(reg['halfwidth'])
    assert eps == F(1,2048) and h == F(1,512)
    assert all(0 < x-h < x+h < 1 for x in c[1:])
    table = json.loads(gzip.decompress(src['vertices017'].read_bytes()))
    nodes = [sum(int(x*D)*z for x,z in zip(c,v)) for v in table['vertices']]
    assert len(set(nodes)) == 280 and min(abs(k) for k in nodes) > 0
    a = 144*D*eps
    A, B = a.numerator, a.denominator
    weights = {}
    for key, w in cache['stencil'][0].items():
        assert int(key) != 0
        k2 = F(int(key)-A*A, B*B)
        assert k2.denominator == 1 and k2 >= 0
        k2 = int(k2)
        assert isqrt(k2)**2 == k2 and k2 in {k*k for k in nodes}
        weights[k2] = weights.get(k2,F(0))+F(w)
    assert len(weights) == 140 and sum(weights.values()) == 0
    # Verify the metric and the box maximum, without invoking symbolic solvers.
    M = [[2,1,1,1,1],[1,-1,-1,-1,2],[1,-1,-1,2,-1],
         [1,-1,2,-1,-1],[1,2,-1,-1,-1]]
    N = [[1,1,1,1,1],[1,-1,-1,-1,1],[1,-1,-1,1,-1],
         [1,-1,1,-1,-1],[1,1,-1,-1,-1]]
    assert all(sum(M[i][k]*N[j][k] for k in range(5)) == 6*(i==j)
               for i in range(5) for j in range(5))
    Q = [[sum(F(M[k][i+1]*M[k][j+1],36) for k in range(5))
          for j in range(4)] for i in range(4)]
    assert all(Q[i][j] == F(i==j,4)-F(1,36) for i in range(4) for j in range(4))
    metric = [sum(s[i]*Q[i][j]*s[j] for i in range(4) for j in range(4))
              for s in product([-1,1],repeat=4)]
    assert max(metric) == 1
    ts = list(map(F,reg['time_endpoints']))
    assert ts == [F(2**i,2048)*(1+F(j,8)) for i in range(11) for j in range(8)]+[F(1)]
    cp = HERE/'checkpoint.json'
    rows = [{'t': str(eps), 'H': low['center_H'], 'source': 'reused063'}]
    previous_segments = []
    if cp.exists():
        assert args.resume, 'Existing checkpoint; use --resume, never silently replay.'
        saved = json.loads(cp.read_text())
        assert saved['fingerprint'] == fingerprint
        rows, previous_segments = saved['rows'], saved['segments']
    receipt = {'pid':os.getpid(),'source_sha256':sha(__file__),
               'registration_sha256':sha(HERE/'registration.json'),
               'fingerprint':fingerprint,'python':sys.executable,'flint':flint.__version__,
               'new_H_evaluations':0,'new_jet_calls':0,'new_cone_evaluations':0}
    factor = arb(D)**4/(144*(2*arb.pi())**5)
    minimal = {'D':D,'B':B,'center':list(map(str,c)),
               'weights':{str(k):str(w) for k,w in sorted(weights.items())},
               'from_cache_sha256':sha(src['cache063'])}
    atomic(HERE/'zero_stencil.json',minimal)
    for t in ts[len(rows):]:
        value = factor*sum((arb(rq(w))*arb(rq(B*B*k2+(144*D*B*t)**2)).log()
                            for k2,w in weights.items()),arb(0))
        assert value.is_finite() and value > 0
        rows.append({'t':str(t),'H':enc(value),'source':'new_scalar_log_evaluation'})
        receipt['new_H_evaluations'] += 1
        atomic(cp,{'fingerprint':fingerprint,'rows':rows,
                   'segments':previous_segments+[{'pid':os.getpid(),
                   'elapsed_seconds':time.monotonic()-began,'completed_endpoints':len(rows)}]})
        if time.monotonic()-began >= reg['stop_check_seconds'] and len(rows)<len(ts):
            receipt['status']='checkpointed_exit75'
            atomic(HERE/'run_receipt.json',receipt)
            raise SystemExit(75)
    assert len(rows) == 89
    bands = []
    L = U = F(0)
    for i,row in enumerate(rows):
        t = F(row['t'])
        R = ( (1+arb(rq(h*h/(20*t*t)))).sqrt()
               +arb(rq(h))/(2*arb(5).sqrt()*arb(rq(t))) )**10
        RL, RU = F(ep(R)), F(ep(R,True))
        if i==0:
            assert RL<=3125<=RU
            # Exact special value avoids an unnecessary outward rounding loss.
            RL=RU=F(3125)
        current = row['H']
        if i+1<len(rows):
            nxt=rows[i+1]['H']
            ml=F(current['lower'])-F(nxt['upper'])
            mu=F(current['upper'])-F(nxt['lower'])
            stop=rows[i+1]['t']
        else:
            ml,mu=F(current['lower']),F(current['upper'])
            stop='infinity'
        assert 0 < ml <= mu and 1 <= RL <= RU
        lc,uc=ml/RU,mu*RU
        L += lc
        U += uc
        bands.append({'start':str(t),'stop':stop,'mass_lower':str(ml),'mass_upper':str(mu),
                      'R_lower':str(RL),'R_upper':str(RU),
                      'contribution_lower':str(lc),'contribution_upper':str(uc)})
    volume = (2*h)**4
    lower = arb(rq(volume))*arb(rq(L))**(arb(5)/4)
    upper = arb(rq(volume))*arb(rq(U))**(arb(5)/4)
    Hlo,Hhi=F(rows[0]['H']['lower']),F(rows[0]['H']['upper'])
    oldscan = next(r for r in low['scans'] if r['h']==str(h))
    comparison = {
        'h':str(h),'old_center':low['center_H'],
        'old_jet_variation':oldscan['variation'],
        'old_absolute_resolvent_remainder':oldscan['all_time_resolvent_remainder'],
        'old_relative_remainder_upper':oldscan['relative_remainder_upper'],
        'old_chord_integral_upper':oldscan['chord_integral_upper'],
        'new_U_over_center_lower':str(U/Hlo),
        'new_L_over_center_upper':str(L/Hhi),
        'old_chord_over_new_constant_upper':str(F(oldscan['chord_integral_upper'])/F(ep(upper,True)))
    }
    result={'fingerprint':fingerprint,'center':list(map(str,c)),
            'halfwidth':str(h),'box_lower':list(map(str,(x-h for x in c[1:]))),
            'box_upper':list(map(str,(x+h for x in c[1:]))),'volume':str(volume),
            'eps':str(eps),'whole_box_H_lower':str(L),'whole_box_H_upper':str(U),
            'power_integral_lower':ep(lower),'power_integral_upper':ep(upper,True),
            'power_integral_lower_expression':str(lower),'power_integral_upper_expression':str(upper),
            'meaningful_gate_passed':U<=F(reg['meaningful_ratio_gate'])*Hlo,
            'ideal_gate_passed':U<=F(reg['ideal_ratio_gate'])*Hlo,
            'finite_time_endpoint_count':len(rows),'positive_band_count':len(bands),
            'time_domain':'[1/2048,infinity)','bands':bands,'comparison063':comparison,
            'metric_corner_values':list(map(str,metric)),
            'new_global_I_or_g_or_C5':False}
    atomic(HERE/'results.json',result)
    receipt.update({'status':'complete','elapsed_seconds':time.monotonic()-began,
                    'peak_rss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                    'endpoint_checkpoint_sha256':sha(cp),'results_sha256':sha(HERE/'results.json'),
                    'zero_stencil_sha256':sha(HERE/'zero_stencil.json')})
    atomic(HERE/'run_receipt.json',receipt)
    summary={'status':'complete','seconds':receipt['elapsed_seconds'],
             'H_center':float(Hlo),'L':float(L),'U':float(U),'U_over_H':float(U/Hlo),
             'integral_upper':float(F(ep(upper,True))),
             'old_R_over_H':float(F(oldscan['relative_remainder_upper'])),
             'old_chord_over_new':float(F(comparison['old_chord_over_new_constant_upper'])),
             'meaningful':result['meaningful_gate_passed'],'ideal':result['ideal_gate_passed'],
             'new_H_evaluations':receipt['new_H_evaluations']}
    atomic(HERE/'summary.json',summary)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
