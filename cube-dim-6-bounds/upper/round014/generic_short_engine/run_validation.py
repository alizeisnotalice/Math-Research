from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,subprocess,time,sys
from flint import arb,ctx,fmpq
sys.set_int_max_str_digits(0);ctx.prec=1024
HERE=Path(__file__).resolve().parent
assert json.loads((HERE/'anchor_validation.json').read_text())['status']=='PASS_ALL6_PHI_DD_PARTIAL_AND_P_IDENTITY_CASES'
ROOT=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
GFILE=ROOT/'upper/round003/direct_kernel/local_minimum/exact_local_reference_results.json'
HFILE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-theory/cube-dim-6-bounds/upper/round007/g_time_moments/results.json')
DGFILE=HERE.parents[1]/'round013/generic_green_engine/gradient_exact_reference.json'
DHFILE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round010/short_minimum_diagnostic/gradient_results.json')
FAN=ROOT/'upper/round003/direct_kernel/fan_integer_input.txt'
g=json.loads(GFILE.read_text());h=json.loads(HFILE.read_text());dg=json.loads(DGFILE.read_text());dh=json.loads(DHFILE.read_text())
assert g['point_numerators']==h['center_numerators']==dg['point_numerators']==dh['center_numerators']
assert g['denominator']==h['denominator']==dg['denominator']==dh['denominator']
assert hashlib.sha256(GFILE.read_bytes()).hexdigest()==dg['old_G_certificate_sha256']
assert hashlib.sha256(FAN.read_bytes()).hexdigest()==g['fan_input_sha256']==dg['fan_sha256']
hr=next(v for v in h['time_nodes']if v['tau']=='4/125');dhr=next(v for v in dh['gradient_passes']if v['tau']=='4/125')
assert hr['complete_group_terms']==dhr['complete_group_terms']==967680
assert(arb(hr['H_arb384'])-arb(dhr['H_arb384'])).contains(0)
G=arb(fmpq(int(g['Q_over_pi5_numerator']),int(g['Q_over_pi5_denominator'])))/arb.pi()**5
DG=arb(fmpq(int(dg['gradient_Q_over_pi5_numerator']),int(dg['gradient_Q_over_pi5_denominator'])))/arb.pi()**5
Kref=G-arb(hr['H_arb384']);Dref=DG-arb(dhr['gradient_x2_to_x6'][0])
cmd=[str(HERE/'scratch/generic_short'),'full',str(FAN),'4','125',str(g['denominator'])]+list(map(str,g['point_numerators']))
tic=time.perf_counter();z=subprocess.run(cmd,capture_output=True,text=True,timeout=300);elapsed=time.perf_counter()-tic
assert z.returncode==0,z.stderr
native=json.loads(z.stdout);assert native['status']=='generic_nonzero_supported'and native['full_group_terms']==967680
(HERE/'native_result.json').write_text(json.dumps(native,indent=2)+'\n')

def fb(row):
    factor=Q(2)**int(row['common_binary_exponent']);return[Q(int(row[k]))*factor for k in ['lower_mantissa','upper_mantissa']]
def b(v):return arb(fmpq(v.numerator,v.denominator))
def fullball(row):
    lo,hi=fb(row);return arb(b((lo+hi)/2),b((hi-lo)/2))
kl,ku=fb(native['K']);dl,du=fb(native['gradient'][1]);assert b(kl)>Kref.lower()and b(ku)<Kref.upper(),'STOP: incompatible exactG-H32 K interval'
assert b(dl)>Dref.lower()and b(du)<Dref.upper(),'STOP: incompatible exactgradG-oldgradH32 x2 interval'
el,eu=fb(native['Euler_x_grad_plus5K_plusTp']);assert el<=0<=eu
euler=5*fullball(native['K'])+arb(4)/125*fullball(native['p_T'])+sum((b(Q(v,g['denominator']))*fullball(d)for v,d in zip(g['point_numerators'],native['gradient'])),arb(0));assert euler.contains(0)
pl,pu=fb(native['p_T']);assert pl>0
result={'status':'PASS_GENERIC_K_GRADIENT_P_AND_SOURCE_CROSSCHECK','full_evaluator_executions':1,'new_spatial_points':0,'T':'4/125',
    'native_exit':z.returncode,'native_external_seconds':elapsed,'native_kernel_seconds':native['kernel_seconds'],
    'native_K_strictly_inside_old_exactG_minus_H32_certificate':True,'native_gradient_x2_strictly_inside_C013_exactDG_minus_B010_DH32_certificate':True,
    'Euler_native_and_independent_saved_interval_reconstruction_contain_zero':True,'p_T_strictly_positive':True,
    'K_display_only':fullball(native['K']).str(35),'p_T_display_only':fullball(native['p_T']).str(35),
    'gradient_display_only':[fullball(d).str(35)for d in native['gradient']],
    'input_sha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in [GFILE,HFILE,DGFILE,DHFILE,FAN]},
    'source_sha256':hashlib.sha256((HERE/'generic_short.cpp').read_bytes()).hexdigest(),'binary_sha256':hashlib.sha256((HERE/'scratch/generic_short').read_bytes()).hexdigest(),
    'command':cmd,'stopped':True,'live_handles':[]}
(HERE/'validation_results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
