"""One fixed analytic six-profile witness batch, three configurations."""
import json,hashlib,os,sys,time,signal,resource,traceback
from pathlib import Path
from fractions import Fraction as Q
from math import comb
sys.set_int_max_str_digits(0)
resource.setrlimit(resource.RLIMIT_CPU,(60,61))
def halt(sig,frame):raise TimeoutError('CPU60 or wall175 stop')
signal.signal(signal.SIGXCPU,halt);signal.signal(signal.SIGALRM,halt);signal.alarm(175)
START=time.monotonic();CPU=time.process_time();H=Path(__file__).resolve().parent
OWN=H.parents[2];MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
rd=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(n,d):(H/n).write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n')
def pair(a,b):return {'lower':str(a),'upper':str(b),'display_only':[float(a),float(b)]}
if (H/'run_receipt.json').exists():raise RuntimeError('One attempt only')
reg=rd(H/'registration.json');par=rd(H/'parameters.json');scale=10**par['integer_root_decimal_scale']
out={'status':'INCOMPLETE','configurations':[],'root_certificates':[],'profiles':{}}
receipt={'status':'RUNNING','pid':os.getpid(),'attempted_configurations':0,'completed_configurations':0,
         'old_function_replays':0,'body_kernel_calls':0,'automatic_retries':0}
def record():
    receipt.update(elapsed_seconds=time.monotonic()-START,cpu_seconds=time.process_time()-CPU)
    save('run_receipt.json',receipt);save('results.json',out)
def begin(i):
    assert i==receipt['attempted_configurations']+1<=3
    receipt['attempted_configurations']=i;record()
def done(i,name):
    out['configurations'].append({'id':i,'name':name,'pass':True})
    receipt['completed_configurations']=i;record()
def pins():
    for p,h in reg['source_sha256'].items():assert sha(p)==h,p
    assert sha(H/'parameters.json')==reg['parameters_sha256'] and sha(H/'theorem.md')==reg['theorem_sha256']
    for n,h in rd(H/'execution_binding.json').items():assert sha(H/n)==h,n
def root(x,n,label):
    assert x>=0
    if x==0:lo=hi=Q(0)
    else:
        target=(x.numerator*scale**n)//x.denominator
        a=0;b=1<<((target.bit_length()+n-1)//n+1)
        while b-a>1:
            m=(a+b)//2
            if m**n<=target:a=m
            else:b=m
        lo,hi=Q(a,scale),Q(a+1,scale)
    assert lo**n<=x<=hi**n
    out['root_certificates'].append({'label':label,'degree':n,'radicand':str(x),'lower':str(lo),'upper':str(hi)})
    return lo,hi
try:
    record();pins();begin(1)
    prior=rd(OWN/'upper/round092/middle_envelope_gap/lower_results.json')
    assert prior['status']=='STRICT_FIXED_ENVELOPE_NONCOMPETITIVE'
    q,b=Q(par['q']),Q(par['b']);assert q==Q(6,5) and b==Q(1,64) and q/4<1
    cJ,cW=Q(prior['cJ_lower']),Q(prior['cW_lower']);assert min(cJ,cW)>0
    coeff=[]
    for k in range(1,6):
        p=prior['f_coefficients'][str(k)]
        f,v,n=Q(p['uniform_F_lower']),Q(p['V_lower']),p['n']
        assert min(f,v)>0 and n==comb(5,k-1)
        coeff.append((f'k{k}',f,v,cJ,n))
    whole=('whole_union',max(x[1] for x in coeff),sum(x[2]*x[4] for x in coeff),cW,1)
    coeff.append(whole)
    assert 1/(q+1)-1/(q+2)==1/((q+1)*(q+2))
    out['source_coefficients']=[{'name':n,'f':str(f),'v':str(v),'c':str(c),'multiplicity':m} for n,f,v,c,m in coeff]
    out['full_time_integrability_exponent']=str(q/4)
    done(1,'saved_coefficients_and_abstract_witness_admissibility_gates')
    begin(2);bpow=root(b**16,5,'b-to-q-plus2');witness={}
    for name,f,v,c,mult in coeff:
        tau=root(v*f/c,2,name+'-space-cutoff')
        m=min(b,tau[0]);assert 0<m<=b
        used_volume=c*m*m/f;assert 0<used_volume<=v
        z=b-m
        fpower=root(f,5,name+'-f-to-q-minus1')
        z11=root(z**11,5,name+'-z-to-q-plus1')
        z16=root(z**16,5,name+'-z-to-q-plus2')
        bracket_lo=bpow[0]/((q+1)*(q+2))-b*z11[1]/(q+1)+z16[0]/(q+2)
        bracket_hi=bpow[1]/((q+1)*(q+2))-b*z11[0]/(q+1)+z16[1]/(q+2)
        assert 0<bracket_lo<=bracket_hi
        lo=2*c*fpower[0]*bracket_lo;hi=2*c*fpower[1]*bracket_hi
        assert 0<lo<=hi
        witness[name]=(lo,hi)
        out['profiles'][name]={'multiplicity':mult,'f':str(f),'input_v':str(v),'c':str(c),'b':str(b),
            'rational_cutoff_m':str(m),'actual_used_volume':str(used_volume),
            'admissible_height_mass_volume':True,'q_integral':pair(lo,hi),
            'scope':'actual full-time integral of the constructed abstract profile, not the real Cauchy density'}
    done(2,'six_fixed_analytic_profile_q_integrals')
    begin(3)
    tube_lo=sum(comb(5,k-1)*witness[f'k{k}'][0] for k in range(1,6))
    tube_hi=sum(comb(5,k-1)*witness[f'k{k}'][1] for k in range(1,6))
    whole_lo,whole_hi=witness['whole_union'];new_lo=min(tube_lo,whole_lo)
    B_tube=root(tube_lo**5,6,'tube-relaxation-B-lower')[0]
    B_whole=root(whole_lo**5,6,'whole-relaxation-B-lower')[0]
    B_new=root(new_lo**5,6,'new-min-relaxation-B-lower')[0]
    ledger=rd(MAIN/'upper/round062/primary_joint_account/actual_run/results.json')
    allowanceB=Q(ledger['sufficient_remaining_D_F_B_allowance'])
    allowanceI=Q(ledger['sufficient_remaining_D_F_I_allowance']);paid=Q(ledger['paid_C_E_H_I_upper'])
    g0=Q(par['g0']);target=root((g0/6)**6,5,'conditional-I-target-unit')
    target=(64*target[0],64*target[1]);Iplanlo=paid+new_lo
    factor=root((6/g0)**6,5,'conditional-U-factor')
    S_upper=target[1]-Iplanlo;U_lower=factor[0]*Iplanlo
    outcome='STRICT_PAIRED_JOINT_RELAXATION_NONCOMPETITIVE' if new_lo>allowanceI and B_new>allowanceB and S_upper<0 else 'UNRESOLVED_STOP'
    out.update(status=outcome,tube_witness_I_sum=pair(tube_lo,tube_hi),whole_witness_I=pair(whole_lo,whole_hi),
       new_min_relaxation_I_lower=str(new_lo),tube_relaxation_B_lower=str(B_tube),whole_relaxation_B_lower=str(B_whole),new_min_relaxation_B_lower=str(B_new),
       saved_D_F_I_allowance=str(allowanceI),saved_D_F_B_allowance=str(allowanceB),
       I_allowance_gap=str(new_lo-allowanceI),B_allowance_gap=str(B_new-allowanceB),
       conditional_g0=str(g0),global_g0_proven=False,paid_CEH_I_upper=str(paid),
       full_time_D_requires_no_separate_short_addition=True,conditional_plan_I_lower=str(Iplanlo),
       conditional_I_target=pair(*target),conditional_plan_U_lower=str(U_lower),conditional_plan_S_upper=str(S_upper),
       complete_abstract_full_time_inequality_proven=True,competitive_D_budget_paid=False,far_budget_paid=False,
       actual_density_witness=False,actual_model_cost_lower_bound=False,new_C6_bound=False,
       scope='Only min of separate sharp per-tube and whole-union time-rank relaxations; not optimum under all simultaneous overlap constraints.')
    done(3,'both_candidate_lower_certificates_and_conditional_cost_margin')
    pins();receipt.update(status='COMPLETED_STOPPED',exit_code=0);record()
    receipt['results_sha256']=sha(H/'results.json');save('run_receipt.json',receipt)
    print(json.dumps({'status':outcome,'pid':os.getpid(),'CPU':receipt['cpu_seconds'],'wall':receipt['elapsed_seconds'],
        'tube_B_lower_display':float(B_tube),'whole_B_lower_display':float(B_whole),'new_B_lower_display':float(B_new),
        'B_allowance_display':float(allowanceB),'U_lower_display':float(U_lower),'S_upper_display':float(S_upper),
        'root_certificates':len(out['root_certificates'])}))
except BaseException as ex:
    receipt.update(status='FAILED_STOP_NO_RETRY',exit_code=1,error=repr(ex),traceback=traceback.format_exc());record();raise
finally:signal.alarm(0)
