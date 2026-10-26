"""Registered configuration 1: no F/L/envelope evaluation, saved costs only."""
import os,sys,time,signal,resource,json,hashlib
from pathlib import Path
from fractions import Fraction as Q
sys.set_int_max_str_digits(0)
resource.setrlimit(resource.RLIMIT_CPU,(30,30));signal.alarm(170)
start=time.monotonic();cpu=time.process_time();H=Path(__file__).resolve().parent
rd=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(n,d):(H/n).write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n')
if (H/'cost_receipt.json').exists():raise RuntimeError('Single attempt only')
reg=rd(H/'registration.json');params=rd(H/'parameters.json')
assert sha(H/'parameters.json')==reg['parameters_sha256']
for p,h in reg['source_sha256'].items():assert sha(p)==h,p
scale=10**params['integer_root_decimal_scale'];certs=[]
def root(x,n,label):
    assert x>=0
    target=(x.numerator*scale**n)//x.denominator
    lo=0;hi=1<<((target.bit_length()+n-1)//n+1)
    while hi-lo>1:
        m=(lo+hi)//2
        if m**n<=target:lo=m
        else:hi=m
    a,b=Q(lo,scale),Q(lo+1,scale)
    assert a**n<=x<b**n
    certs.append({'label':label,'degree':n,'radicand':str(x),'lower':str(a),'upper':str(b)})
    return a,b
def ep(v):return Q(v['lower_mantissa'])*Q(2)**v['lower_binary_exponent'],Q(v['upper_mantissa'])*Q(2)**v['upper_binary_exponent']
def pair(a,b):return {'lower':str(a),'upper':str(b),'display_only':[float(a),float(b)]}
try:
    save('cost_receipt.json',{'status':'RUNNING','pid':os.getpid(),'configuration':1})
    old=rd(H.parent.parent/'round074/middle_long_budget/results.json')
    ledger=rd('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round062/primary_joint_account/actual_run/results.json')
    Rlo,Rhi=ep(old['long_budget_expression']);short=Q(old['comparison']['old_short_upper_read_only'])
    paid=Q(ledger['paid_C_E_H_I_upper']);allow=Q(ledger['sufficient_remaining_D_F_B_allowance'])
    g=Q(params['g0_conditional']);factor=root((6/g)**6,5,'S-per-unit-I factor')
    ILo=root((short+Rlo)**6,5,'old-middle-plan-I-lower')[0]
    IHi=root((short+Rhi)**6,5,'old-middle-plan-I-upper')[1]
    BLo=root((paid+ILo)**5,6,'old-paid-plus-middle-B-lower')[0]
    BHi=root((paid+IHi)**5,6,'old-paid-plus-middle-B-upper')[1]
    SLo=factor[0]*(paid+ILo);SHi=factor[1]*(paid+IHi)
    CLo=factor[0]*paid;CHi=factor[1]*paid
    out={'status':'PASS_SAVED_UPPER_PLAN_COST_CONVERSION','configuration':1,
         'q':'6/5','S_definition':'S=(6*B/g0)^(6/5)=6^(6/5)*I/g0^(6/5)',
         'g0':str(g),'g0_globally_proven':False,'actual_C6_bound':False,
         'interpretation':'S interval describes the existing conservative cost plan with paid C/E/H, old short+rectangle long on D, and zero far. It is not an actual-kernel lower bound.',
         'saved_long_upper_expression':pair(Rlo,Rhi),'saved_short_upper':str(short),
         'saved_paid_C_E_H_I_upper':str(paid),'saved_remaining_D_F_B_allowance':str(allow),
         'S_factor_per_unit_I':pair(*factor),'old_middle_plan_I':pair(ILo,IHi),
         'old_paid_plus_middle_plan_I':pair(paid+ILo,paid+IHi),
         'old_paid_plus_middle_plan_B':pair(BLo,BHi),
         'old_paid_plus_middle_plan_S':pair(SLo,SHi),
         'paid_C_E_H_contribution_to_S':pair(CLo,CHi),
         'root_certificates':certs,'new_envelope_endpoint_evaluations':0,
         'old_envelope_function_recalculations':0,'body_kernel_calls':0,
         'same_E_integral_lower_yet_available':False}
    save('cost_results.json',out)
    for p,h in reg['source_sha256'].items():assert sha(p)==h,p
    receipt={'status':'PASS','pid':os.getpid(),'configuration':1,'exit_code':0,
             'elapsed_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,
             'new_envelope_endpoint_evaluations':0,'body_kernel_calls':0,'results_sha256':sha(H/'cost_results.json')}
    save('cost_receipt.json',receipt)
    print(json.dumps({'receipt':receipt,'old_plan_S_display':out['old_paid_plus_middle_plan_S']['display_only']}))
except BaseException as ex:
    save('cost_receipt.json',{'status':'FAILED_STOP_NO_RETRY','pid':os.getpid(),'error':repr(ex),'elapsed_seconds':time.monotonic()-start});raise
finally:signal.alarm(0)
