"""One registered old point only. Does not start the independent exact full pass."""
from pathlib import Path
from fractions import Fraction as Q
import subprocess,time,json,hashlib,sys
from flint import arb,ctx,fmpq
sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent
assert json.loads((HERE/'anchor_validation.json').read_text())['status'].startswith('PASS_ALL6')
SOURCE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round003/direct_kernel')
OLD=SOURCE/'local_minimum/exact_local_reference_results.json'
old=json.loads(OLD.read_text());fan=SOURCE/'fan_integer_input.txt'
cmd=[str(HERE/'scratch/generic_green'),'full',str(fan),str(old['denominator'])]+list(map(str,old['point_numerators']))
tic=time.perf_counter();p=subprocess.run(cmd,capture_output=True,text=True,timeout=300);elapsed=time.perf_counter()-tic
assert p.returncode==0,p.stderr
native=json.loads(p.stdout);assert native['status']=='generic_nonzero_supported'and native['full_group_terms']==967680
(HERE/'native_result.json').write_text(json.dumps(native,indent=2)+'\n')

def exact_bounds(row):
    e=int(row['common_binary_exponent']);return[Q(int(row[k]))*Q(2)**e for k in ['lower_mantissa','upper_mantissa']]

ctx.prec=1024
g=(arb(fmpq(int(old['Q_over_pi5_numerator']),int(old['Q_over_pi5_denominator'])))/arb.pi()**5)
l,u=exact_bounds(native['G']);lo=arb(fmpq(l.numerator,l.denominator));hi=arb(fmpq(u.numerator,u.denominator))
assert lo<g.lower()and hi>g.upper(),'STOP: G exact certificate disagreement'
el,eu=exact_bounds(native['Euler_x_grad_plus5G']);assert el<=0<=eu
result={'status':'PASS_NATIVE_G_AND_EULER; GRADIENT_X2_INDEPENDENT_REFERENCE_PENDING',
    'full_evaluator_executions_so_far':1,'new_spatial_points':0,'native_returncode':p.returncode,
    'native_external_seconds':elapsed,'native_kernel_seconds':native['kernel_seconds'],
    'old_value_only_fmpq_seconds':old['elapsed_seconds'],
    'old_value_only_over_native_value_plus6gradient_time':old['elapsed_seconds']/elapsed,
    'native_G_contains_exact_old_Q_over_pi5':True,'Euler_interval_contains_zero':True,
    'old_exact_source_sha256':hashlib.sha256(OLD.read_bytes()).hexdigest(),'fan_sha256':hashlib.sha256(fan.read_bytes()).hexdigest(),
    'native_source_sha256':hashlib.sha256((HERE/'generic_green.cpp').read_bytes()).hexdigest(),
    'binary_sha256':hashlib.sha256((HERE/'scratch/generic_green').read_bytes()).hexdigest(),'command':cmd}
(HERE/'native_validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
