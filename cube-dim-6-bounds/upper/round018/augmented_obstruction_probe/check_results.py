#!/usr/bin/env python3
"""B018 saved-data checker: exact dyadics and true 38-column S2, zero kernels."""
import sys
sys.dont_write_bytecode = True
import hashlib
import itertools
import json
import math
import os
from decimal import Decimal, localcontext
from fractions import Fraction as Q
from pathlib import Path
from time import perf_counter

HERE = Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def interval(obj):
    assert set(obj) == {'lower_mantissa','upper_mantissa','common_binary_exponent'}
    lo, hi, ex = (int(obj[k]) for k in ('lower_mantissa','upper_mantissa','common_binary_exponent'))
    assert lo <= hi and abs(ex) <= 10000
    den, mul = 1 << max(0,-ex), 1 << max(0,ex)
    return Q(lo*mul,den), Q(hi*mul,den)

def dec(q):
    with localcontext() as ctx:
        ctx.prec = 65
        return str(Decimal(q.numerator)/Decimal(q.denominator))

def main():
    tic = perf_counter()
    manifest = json.loads((HERE/'input_manifest.json').read_text())
    receipt = json.loads((HERE/'run_receipt.json').read_text())
    raw = json.loads((HERE/'native.stdout.json').read_text())
    for row in manifest['inputs']:
        assert sha(row['path']) == row['sha256']
    assert sha(HERE/'input_manifest.json') == receipt['input_manifest_sha256']
    assert sha(HERE/'native.stdout.json') == receipt['stdout_sha256']
    assert sha(HERE/'native.stderr.txt') == receipt['stderr_sha256']
    assert (HERE/'native.stderr.txt').read_text() == ''
    assert receipt['status'] == 'SINGLE_NATIVE_EXIT0_STOPPED'
    assert receipt['native_attempts'] == 1 and receipt['returncode'] == 0
    assert receipt['elapsed_seconds'] < 120 and receipt['all_native_children_waited']
    assert receipt['started_epoch']-manifest['started_epoch']+receipt['elapsed_seconds'] < 1200
    assert manifest['T'] == '2/45' and manifest['beta'] == '1'
    assert manifest['inherits_base_kernel_value'] is False
    vs = [(1,*s) for s in itertools.product((-1,1),repeat=5)] + [tuple(int(i==j) for i in range(6)) for j in range(6)]
    assert [list(v) for v in vs] == manifest['model']['directions']
    frame_sha = hashlib.sha256(json.dumps(manifest['model'],sort_keys=True,separators=(',',':')).encode()).hexdigest()
    assert frame_sha == manifest['frame_sha256']
    assert raw['status'] == 'generic_nonzero_supported'
    assert raw['precision_bits'] == 512 and raw['full_group_terms'] == 967680
    assert raw['computed_sign_pairs'] == 483840 and raw['zero_node_terms'] == 0
    assert 0 <= raw['confluent_terms'] <= 967680 and raw['confluent_terms'] % 2 == 0
    ns, den = manifest['numerators'], manifest['denominator']
    x = [Q(n,den) for n in ns]
    assert x[0] == 1 and all(x[i] > x[i+1] > 0 for i in range(5))
    k, p = interval(raw['K']), interval(raw['p_T'])
    gradient = [interval(v) for v in raw['gradient']]
    assert len(gradient) == 6
    native_euler = interval(raw['Euler_x_grad_plus5K_plusTp'])
    assert k[0] > 0 and p[0] > 0 and native_euler[0] <= 0 <= native_euler[1]
    for pair in [k,p,*gradient]:
        assert pair[1]-pair[0] < max(Q(1),abs(sum(pair)/2))*Q(2)**-256
    T = Q(2,45)
    elo,ehi = 5*k[0]+T*p[0],5*k[1]+T*p[1]
    for xx,gg in zip(x,gradient):
        elo += xx*gg[0]
        ehi += xx*gg[1]
    assert elo <= 0 <= ehi
    m64,m625 = Q(1,92160),Q(1,90000)
    assert T/4096 == m64 and T/4000 == m625
    assert k[1] < m64 < m625
    denominator = 10**16
    lower = Q(k[0].numerator*denominator//k[0].denominator,denominator)
    upper = Q(-((-k[1].numerator*denominator)//k[1].denominator),denominator)
    assert lower < k[0] <= k[1] < upper < m64
    # Actual new S2: all 38 choose 2 pairs, including both kinds of axis pairs.
    norm2 = sum(y*y for y in x)
    candidates = []
    pair_types = {'body_body':0,'body_axis':0,'axis_axis':0}
    for i,v in enumerate(vs):
        a = sum(t*t for t in v)
        xv = sum(t*z for t,z in zip(v,x))
        for j in range(i):
            w = vs[j]
            b = sum(t*t for t in w)
            xw = sum(t*z for t,z in zip(w,x))
            c = sum(t*s for t,s in zip(v,w))
            determinant = a*b-c*c
            assert determinant > 0
            d2 = norm2-(b*xv*xv-2*c*xv*xw+a*xw*xw)/determinant
            candidates.append((d2,j,i))
            typ = 'body_body' if i < 32 else ('body_axis' if j < 32 else 'axis_axis')
            pair_types[typ] += 1
    assert len(candidates) == 703 and pair_types == {'body_body':496,'body_axis':192,'axis_axis':15}
    distance2,j,i = min(candidates)
    assert distance2 > 0
    with localcontext() as ctx:
        ctx.prec = 45
        distance_display = str((Decimal(distance2.numerator)/Decimal(distance2.denominator)).sqrt())
    q = (k[0]/T,k[1]/T)
    dq = ((T*p[0]-k[1])/T**2,(T*p[1]-k[0])/T**2)
    result = {'status':'PASS_STRICT_FIXED_BETA1_T2_OVER45_OBSTRUCTION',
        'pid':os.getpid(),'elapsed_seconds':perf_counter()-tic,'additional_kernel_calls':0,
        'native_calls_in_batch':1,'beta':'1','T':'2/45','frame_sha256':frame_sha,
        'numerators':ns,'denominator':den,'K_interval_exact':[str(v) for v in k],
        'K_midpoint_display':dec(sum(k)/2),'K_width_exact':str(k[1]-k[0]),'K_width_display':dec(k[1]-k[0]),
        'coarse_rational_K_enclosure':[str(lower),str(upper)],
        'm64_minus_K_upper_exact':str(m64-k[1]),'m64_minus_K_upper_display':dec(m64-k[1]),
        'm625_minus_K_upper_exact':str(m625-k[1]),'rejects_strict64':True,'rejects625':True,
        'point_cost_interval_exact':[str(T/(64*k[1])),str(T/(64*k[0]))],
        'point_cost_lower_display':dec(T/(64*k[1])),
        'q_interval_exact':[str(v) for v in q],'q_midpoint_display':dec(sum(q)/2),
        'dq_dT_interval_exact':[str(v) for v in dq],'dq_dT_midpoint_display':dec(sum(dq)/2),
        'dq_dT_strict_sign':'negative' if dq[1]<0 else ('positive' if dq[0]>0 else 'unresolved'),
        'gradient_intervals_exact':[[str(a),str(b)]for a,b in gradient],
        'p_interval_exact':[str(v)for v in p],
        'independent_Euler_interval_exact':[str(elo),str(ehi)],
        'true_augmented_S2':{'pair_count':703,'pair_types':pair_types,'distance_squared_exact':str(distance2),
                            'distance_display':distance_display,'closest_pair_indices':[j,i],
                            'strictly_off_true_S2':True},
        'native_confluent_terms':raw['confluent_terms'],'native_zero_terms':raw['zero_node_terms'],
        'new_C6_bound':False,'global_minimum_certified':False,
        'scope':'only fixed beta=1/T=2/45 is obstructed by this continuous face point',
        'inputs':[{'path':str(HERE/f),'sha256':sha(HERE/f)}for f in
                  ('registration.md','probe.py','input_manifest.json','run_receipt.json','native.stdout.json','native.stderr.txt','check_results.py')]}
    (HERE/'checked_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ('status','pid','elapsed_seconds','K_midpoint_display','K_width_display',
          'coarse_rational_K_enclosure','m64_minus_K_upper_display','point_cost_lower_display','q_midpoint_display',
          'dq_dT_midpoint_display','dq_dT_strict_sign','true_augmented_S2','native_confluent_terms')},indent=2))

if __name__ == '__main__':
    main()
