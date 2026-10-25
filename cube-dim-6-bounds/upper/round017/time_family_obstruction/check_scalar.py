#!/usr/bin/env python3
"""Twelve or fewer exact scalar checks; never evaluates any kernel."""
from pathlib import Path
from fractions import Fraction as Q
from time import perf_counter
import hashlib,json,os,signal

HERE=Path(__file__).resolve().parent
B015=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round015/spatial_minimum_search')
NATIVE=B015/'native_outputs/02.stdout.json'
T0=Q(6,125);LEFT=Q(223,5000);RIGHT=Q(251,5000)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(name,value):(HERE/name).write_text(json.dumps(value,indent=2)+'\n')
def interval(row):
    power=Q(2)**int(row['common_binary_exponent'])
    return tuple(Q(int(row[k]))*power for k in ['lower_mantissa','upper_mantissa'])

start=perf_counter()
sources=[HERE/'registration.md',HERE/'theorem.md',Path(__file__),NATIVE,
         B015/'checkpoint.json',B015/'checked_results.json',B015/'report.md']
manifest={'status':'REGISTERED_BEFORE_SCALAR_ARITHMETIC','new_kernel_calls':0,
          'max_named_scalar_checks':12,'max_seconds':1200,'T0':str(T0),
          'interval':[str(LEFT),str(RIGHT)],'preferred_c':'1/4100','fallback_c':'1/4096',
          'inputs':[{'path':str(p),'sha256':sha(p)}for p in sources]}
write('input_manifest.json',manifest)
signal.signal(signal.SIGALRM,lambda a,b:(_ for _ in ()).throw(TimeoutError('1200-second scalar cap')))
signal.alarm(1200)
try:
    native=json.loads(NATIVE.read_text());checkpoint=json.loads((B015/'checkpoint.json').read_text())
    record=checkpoint['native_results'][1]
    assert record['stdout_sha256']==sha(NATIVE)
    assert record['T']=='6/125'and record['numerators']==[100000000000000,99999998999999,99999997999988,90268521812170,53714088833376,157107163435]
    assert record['denominator']==100000000000000
    assert native['status']=='generic_nonzero_supported'and native['full_group_terms']==967680
    L,U=interval(native['K']);l,u=interval(native['p_T'])
    assert [str(L),str(U)]==record['check']['K_interval']
    assert [str(l),str(u)]==record['check']['p_interval']
    checks=[]
    def check(name,value):
        checks.append({'name':name,'PASS':bool(value)})
        assert len(checks)<=12
        return bool(value)
    assert check('positive_ordered_K_interval',0<L<=U)
    assert check('positive_ordered_p_interval',0<l<=u)
    assert check('registered_interval_straddles_T0',0<LEFT<T0<RIGHT)
    assert check('existing_q_upper_below_1_over_4096',U/T0<Q(1,4096))
    assert check('existing_p_upper_below_q_lower',u<L/T0)
    assert check('both_A_coefficients_positive',U-u*T0/33>0)
    qleft=(U-l*T0/33*(1-(LEFT/T0)**33))/LEFT
    qright=(U+u*T0/33*((RIGHT/T0)**33-1))/RIGHT
    preferred=check('left_envelope_below_1_over_4100',qleft<Q(1,4100))
    preferred=check('right_envelope_below_1_over_4100',qright<Q(1,4100))and preferred
    assert check('left_envelope_below_1_over_4096',qleft<Q(1,4096))
    assert check('right_envelope_below_1_over_4096',qright<Q(1,4096))
    c0=Q(1,4100)if preferred else Q(1,4096)
    assert check('center_q_below_selected_common_bound',U/T0<c0)
    out={'status':'PASS_CONTINUOUS_INTERVAL_AND_SUPPORTED_MIXTURES_ONLY','pid':os.getpid(),
         'named_scalar_checks':checks,'named_scalar_check_count':len(checks),'new_kernel_calls':0,
         'new_G_calls':0,'new_time_evaluations':0,'T0':str(T0),'time_interval':[str(LEFT),str(RIGHT)],
         'common_q_upper':str(c0),'K_interval':list(map(str,[L,U])),'p_interval':list(map(str,[l,u])),
         'left_envelope_q_exact':str(qleft),'right_envelope_q_exact':str(qright),
         'q0_upper_exact':str(U/T0),'q_prime_T0_upper_exact':str((u-L/T0)/T0),
         'display_only':{'q_left_upper':float(qleft),'q_right_upper':float(qright),'q0_upper':float(U/T0),
                         'p0_upper':float(u),'q_prime_T0_upper':float((u-L/T0)/T0)},
         'supported_mixture_model_cost_lower':str(1/(64*c0)),
         'all_T_obstruction_certified':False,'unrestricted_time_mixtures_obstructed':False,
         'global_C6_improvement_impossibility_claimed':False,'beta1_frame_affected':False,
         'source_sha256':sha(__file__),'manifest_sha256':sha(HERE/'input_manifest.json'),
         'elapsed_seconds':perf_counter()-start}
    write('scalar_certificate.json',out)
    print(json.dumps({k:v for k,v in out.items()if k not in ['K_interval','p_interval','left_envelope_q_exact','right_envelope_q_exact','q0_upper_exact','q_prime_T0_upper_exact']},indent=2))
finally:signal.alarm(0)
