#!/usr/bin/env python3
"""Reconstruct only C017's eleven registered rational scalar checks.

No imports of producer code, numerical kernels, Sturm routines or new times.
"""
from pathlib import Path
from fractions import Fraction as Q
import hashlib, json, os, platform, time

HERE = Path(__file__).resolve().parent
TOP = Path('/Users/zhengzhihao/Desktop/T')
SRC = TOP/'Math-Research-cube6-curvature-astra/cube-dim-6-bounds/upper/round017/time_family_obstruction'
B15 = TOP/'Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round015/spatial_minimum_search'
F15 = HERE.parents[1]/'round015/spatial_obstruction_review/results.json'
FILES = {}

def read(p):
    data = p.read_bytes()
    FILES[str(p)] = hashlib.sha256(data).hexdigest()
    return json.loads(data)

def pin(p, expected):
    data = p.read_bytes()
    FILES[str(p)] = hashlib.sha256(data).hexdigest()
    assert FILES[str(p)] == expected, str(p)

def dyadic(row):
    e = int(row['common_binary_exponent'])
    factor = Q(1 << e) if e >= 0 else Q(1, 1 << -e)
    return tuple(int(row[k])*factor for k in ('lower_mantissa','upper_mantissa'))

def main():
    start = time.perf_counter()
    delivery = read(SRC/'delivery_receipt.json')
    for row in delivery['files']:
        pin(SRC/row['name'],row['sha256'])
        assert (SRC/row['name']).stat().st_size == row['bytes']
    meta = read(SRC/'input_manifest.json')
    cert = read(SRC/'scalar_certificate.json')
    for row in meta['inputs']:
        pin(Path(row['path']),row['sha256'])
    assert FILES[str(SRC/'input_manifest.json')] == cert['manifest_sha256']
    assert FILES[str(SRC/'check_scalar.py')] == cert['source_sha256']
    native = read(B15/'native_outputs/02.stdout.json')
    checkpoint = read(B15/'checkpoint.json')
    previous = read(F15)
    pin(F15,'b9aba4a6119f59b89481ac97463992aa1c6f0a20d87fda2d396c48d418d3666e')
    old = previous['results'][1]
    rec = checkpoint['native_results'][1]
    point = [100000000000000,99999998999999,99999997999988,
             90268521812170,53714088833376,157107163435]
    assert rec['numerators'] == old['numerators'] == point
    assert rec['denominator'] == old['denominator'] == 10**14
    assert rec['T'] == cert['T0'] == '6/125'
    assert rec['stdout_sha256'] == old['raw_sha256'] == FILES[str(B15/'native_outputs/02.stdout.json')]
    assert native['status'] == 'generic_nonzero_supported'
    assert native['full_group_terms'] == 967680
    L,U = dyadic(native['K'])
    l,u = dyadic(native['p_T'])
    assert [str(L),str(U)] == cert['K_interval'] == rec['check']['K_interval'] == old['K_exact_dyadic_bounds']
    assert [str(l),str(u)] == cert['p_interval'] == rec['check']['p_interval']
    # Previous independently audited geometry is reused, not recomputed here.
    assert old['S2_distance_squared_exact'] == '6020805515960314567115330493/40000000000000000000000000000'
    t0,a,b = Q(6,125),Q(223,5000),Q(251,5000)
    assert cert['time_interval'] == [str(a),str(b)] == meta['interval']
    c0 = Q(1,4100)
    # Algebraically separated A/T+B*T^32, unlike the producer's bracket formula.
    Aleft,Aright = U-l*t0/33,U-u*t0/33
    Bleft,Bright = l/(33*t0**32),u/(33*t0**32)
    left,right = Aleft/a+Bleft*a**32,Aright/b+Bright*b**32
    q0,dq = U/t0,(u-L/t0)/t0
    named = [
      ('positive_ordered_K_interval',0<L<=U,L),
      ('positive_ordered_p_interval',0<l<=u,l),
      ('registered_interval_straddles_T0',0<a<t0<b,min(a,t0-a,b-t0)),
      ('existing_q_upper_below_1_over_4096',q0<Q(1,4096),Q(1,4096)-q0),
      ('existing_p_upper_below_q_lower',u<L/t0,L/t0-u),
      ('both_A_coefficients_positive',Aright>0,Aright),
      ('left_envelope_below_1_over_4100',left<c0,c0-left),
      ('right_envelope_below_1_over_4100',right<c0,c0-right),
      ('left_envelope_below_1_over_4096',left<Q(1,4096),Q(1,4096)-left),
      ('right_envelope_below_1_over_4096',right<Q(1,4096),Q(1,4096)-right),
      ('center_q_below_selected_common_bound',q0<c0,c0-q0)]
    assert len(named) == cert['named_scalar_check_count'] == 11
    assert [n for n,_,_ in named] == [r['name'] for r in cert['named_scalar_checks']]
    assert all(ok and margin>0 for _,ok,margin in named)
    assert all(r['PASS'] for r in cert['named_scalar_checks'])
    assert str(left) == cert['left_envelope_q_exact']
    assert str(right) == cert['right_envelope_q_exact']
    assert str(q0) == cert['q0_upper_exact']
    assert str(dq) == cert['q_prime_T0_upper_exact'] == old['dq_dT_exact_bounds'][1]
    assert cert['common_q_upper'] == str(c0)
    assert cert['supported_mixture_model_cost_lower'] == str(1/(64*c0)) == '1025/16'
    for key in ('all_T_obstruction_certified','unrestricted_time_mixtures_obstructed',
                'global_C6_improvement_impossibility_claimed','beta1_frame_affected'):
        assert cert[key] is False
    out = {'status':'PASS_CONTINUOUS_INTERVAL_AND_J_SUPPORTED_PROBABILITY_MIXTURES',
           'pid':os.getpid(),'python':platform.python_version(),
           'named_scalar_checks':[{'name':n,'PASS':ok,'strict_margin_exact':str(m)} for n,ok,m in named],
           'named_scalar_check_count':11,'new_kernel_calls':0,'new_G_calls':0,
           'new_time_evaluations':0,'sturm_or_root_calls':0,'MC_calls':0,'MATLAB_calls':0,
           'T0':str(t0),'time_interval':[str(a),str(b)],'common_q_upper':str(c0),
           'left_envelope_q_exact':str(left),'right_envelope_q_exact':str(right),
           'q0_upper_exact':str(q0),'q_prime_T0_upper_exact':str(dq),
           'A_left_exact':str(Aleft),'A_right_exact':str(Aright),
           'model_cost_lower':'1025/16','all_T_obstruction':False,
           'unrestricted_time_mixtures_obstructed':False,'other_frames_obstructed':False,
           'new_C6_bound':False,'input_sha256':FILES,
           'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'elapsed_seconds':time.perf_counter()-start}
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('status','pid','elapsed_seconds','named_scalar_check_count',
                                      'new_kernel_calls','new_time_evaluations')}))

if __name__=='__main__':
    main()
