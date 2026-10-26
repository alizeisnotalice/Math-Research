"""One fixed rigorous long-time budget, nine configurations, no body kernels."""
import sys,time,signal,resource
START=time.monotonic();CPU=time.process_time()
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
resource.setrlimit(resource.RLIMIT_CPU,(30,31))
def halt(sig,frame):raise TimeoutError('074 CPU30s or wall175s computational limit')
signal.signal(signal.SIGXCPU,halt);signal.signal(signal.SIGALRM,halt)
signal.setitimer(signal.ITIMER_REAL,175)
from pathlib import Path
from fractions import Fraction as Q
from math import comb
import json,hashlib,os,traceback
import sympy as sp
from flint import arb,fmpq,ctx
ctx.prec=768
H=Path(__file__).resolve().parent
B=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper')
rd=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(n,o):(H/n).write_text(json.dumps(o,indent=2,ensure_ascii=False)+'\n')
def A(v):
    v=Q(v);return arb(fmpq(v.numerator,v.denominator))
def endpoints(r):return Q(r['lower_mantissa'])*Q(2)**r['lower_binary_exponent'],Q(r['upper_mantissa'])*Q(2)**r['upper_binary_exponent']
def enc(v):
    assert v.is_finite() and v>=0
    a,e=v.lower().man_exp();b,f=v.upper().man_exp()
    return {'lower_mantissa':str(a),'lower_binary_exponent':int(e),'upper_mantissa':str(b),'upper_binary_exponent':int(f),'display_only':str(v)}
def load(r):
    lo,hi=endpoints(r);assert 0<lo<=hi
    v=arb(A((lo+hi)/2),A((hi-lo)/2));assert v.contains(A(lo)) and v.contains(A(hi));return v
def upper(v):return endpoints(enc(v))[1]
def pick(values):
    j=min(range(len(values)),key=lambda i:upper(values[i]));return j,values[j]

if (H/'run_receipt.json').exists():raise RuntimeError('One attempt only; receipt exists')
reg=rd(H/'registration.json');params=rd(H/'parameters.json')
receipt={'status':'RUNNING','pid':os.getpid(),'attempted_configurations':0,'completed_configurations':0,
         'cpu_limit_seconds':30,'wall_limit_including_cleanup':180,'automatic_retries':0,
         'body32_kernel_calls':0,'old_scalar_recomputations':0,'old_grid_replays':0}
out={'status':'INCOMPLETE','precision_bits':768,'configurations':[],'cases':{},'time_slabs':[]}
def record():
    receipt['elapsed_seconds']=time.monotonic()-START;receipt['cpu_seconds']=time.process_time()-CPU
    save('run_receipt.json',receipt);save('results.json',out)
def begin(i):
    assert i==receipt['attempted_configurations']+1<=9
    receipt['attempted_configurations']=i;record()
def done(i,**kw):
    out['configurations'].append({'id':i,'pass':True,**kw});receipt['completed_configurations']=i;record()
def pins():
    assert sha(H/'parameters.json')==reg['parameters_sha256'] and sha(H/'theorem.md')==reg['theorem_sha256']
    for p,d in reg['source_sha256'].items():assert sha(p)==d,p
    for f,d in rd(H/'execution_binding.json')['code_sha256'].items():assert sha(H/f)==d,f
    for rel in ('round043/two_parameter_middle','round049/middle_joint_marginal','round051/all_unsigned_middles'):
        m=rd(B/rel/'manifest.json')
        # Check exact scientific originals used here, not archived workflow paths.
        for f in ('theorem.md','results.json'):
            assert sha(B/rel/f)==m['whitelist_sha256'][f],(rel,f)

try:
    record();pins();old=rd(B/'round051/all_unsigned_middles/results.json')
    old43=rd(B/'round043/two_parameter_middle/results.json')['configurations']
    old49=rd(B/'round049/middle_joint_marginal/results.json')['configurations'][2]
    begin(1)
    z,a=sp.symbols('z a',real=True,positive=True)
    primitive=sp.log(((z+1)**2+a*a)/((z-1)**2+a*a))/(8*(1+a*a))+(sp.atan((z+1)/a)+sp.atan((z-1)/a))/(4*a*(1+a*a))
    diff=sp.cancel(sp.diff(primitive,z)-1/(((z+1)**2+a*a)*((z-1)**2+a*a)))
    assert diff==0
    assert Q(64)/Q(1,16)==1024
    done(1,name='full_time_marginal_primitive_and_strengthened_average',primitive_derivative_residual=str(diff),
         retained_third_denominator='g_k^2+928*t^2',P_old_D_to_full_coefficient_divisor=1024)
    pi=load(old['shared_constants']['pi']);F=load(old['shared_constants']['cauchy_average_factor'])
    crem=load(old['shared_constants']['Crem']);K=A(Q(1000,729))*crem*F
    cases={}
    for cfg,k in enumerate(range(1,6),2):
        begin(cfg);m=6-k;source=old['cases'][str(k)];geom=source if k!=3 else old43[2]
        limits=source if k!=3 else old49
        c={v:load(geom[v]) for v in ('hv0','hv1','hw0','hw1','g')}
        c.update(D=load(source['D_big']),L=load(limits['z_lower']),U=load(limits['z_upper']),n=comb(5,k-1))
        c['V']=pi*pi/4*A(Q(1,64))**4*A(Q(m,k)).sqrt()
        c['Ctail']=c['D']/1024+2*K/(992*960*928)
        cases[k]=c
        out['cases'][str(k)]={'k':k,'m':m,'n':c['n'],'source_h_intervals':'051 full record, or043 configuration3 for k3',
            'old_M_recomputed':False,'V_middle':enc(c['V']),'Ctail':enc(c['Ctail']),
            'full_time_new_F_expression':'D/[1024*t^2*(1/4+t^2)*(1/256+t^2)]+K/(g^2+928*t^2)*(1/[(hv0^2+992*t^2)*(hv1^2+960*t^2)]+1/[(hw0^2+992*t^2)*(hw1^2+960*t^2)])'}
        done(cfg,name='new_long_time_type_and_volume',k=k)
    def boundF(c,t):
        u=t*t
        return c['D']/(1024*u*(A(Q(1,4))+u)*(A(Q(1,256))+u))+K/(c['g']**2+928*u)*(1/((c['hv0']**2+992*u)*(c['hv1']**2+960*u))+1/((c['hw0']**2+992*u)*(c['hw1']**2+960*u)))
    def prim(z,a):return (((z+1)**2+a*a)/((z-1)**2+a*a)).log()/(8*(1+a*a))+(((z+1)/a).atan()+((z-1)/a).atan())/(4*a*(1+a*a))
    def joint(t,L,U):
        v=128*t*t/(pi*pi)*(prim(U,32*t)-prim(L,32*t));assert v>0;return v
    begin(7);slab_total=arb(0);slab_rational=Q(0);chosen_counts=[0,0,0]
    for left,right in params['panels']:
        l,r=Q(left),Q(right);t=A(l);ratio=A(r/l);width=A(r-l)
        fs={k:boundF(cases[k],t) for k in range(1,6)}
        vol=sum((cases[k]['n']*cases[k]['V']*fs[k]**A(Q(6,5)) for k in range(1,6)),arb(0))
        lm={k:ratio**2*joint(t,cases[k]['L'],cases[k]['U']) for k in range(1,6)}
        jointsum=sum((cases[k]['n']*fs[k]**A(Q(1,5))*lm[k] for k in range(1,6)),arb(0))
        fmax=A(max(upper(v) for v in fs.values()))
        mass1=5*ratio**2*joint(t,A(Q(15,64)),A(Q(49,64)))
        mass2=ratio*t/(pi*(1+1024*t*t));_,mass=pick([mass1,mass2])
        whole=fmax**A(Q(1,5))*mass
        chosen,Ibound=pick([vol,jointsum,whole]);norm=Ibound**A(Q(5,6));pay=width*norm
        slab_total+=pay;slab_rational+=upper(pay);chosen_counts[chosen]+=1
        out['time_slabs'].append({'left':left,'right':right,'chosen':['tube_volume_union','tube_joint_mass_union','one_whole_union_mass'][chosen],
                                 'norm_upper_rectangle':enc(norm),'integrated_rectangle':enc(pay)})
    out['slab_total']=enc(slab_total);out['slab_exact_upper_sum']=str(slab_rational)
    done(7,name='fixed86_interval_union_bound',slab_count=len(params['panels']),selected_envelope_counts=chosen_counts,
         spatial_grid=False,adaptive_refinement=False,full_kernel_calls=0)
    begin(8);T=Q(params['T'])
    coefficient=sum((c['n']*c['V']*c['Ctail']**A(Q(6,5)) for c in cases.values()),arb(0))**A(Q(5,6))
    tail=coefficient/(5*A(T)**5);total=slab_total+tail
    out['infinite_tail']=enc(tail);out['long_budget_expression']=enc(total)
    rational=slab_rational+upper(tail);simple=Q((rational*10**18).__ceil__()+1,10**18)
    out['strict_long_norm_rational_upper']=str(simple)
    done(8,name='complete_long_time_to_infinity',time_domain='[1/1000,infinity)',all_actual_events_paid=True,
         overlaps_paid_positively=True,strict_rational_upper=str(simple))
    begin(9);ledger=rd(B/'round062/primary_joint_account/actual_run/results.json')
    allowance=Q(ledger['sufficient_remaining_D_F_B_allowance'])
    short=Q(old['union']['simple_rational_union_upper']);combined=short+simple
    lo,hi=endpoints(out['long_budget_expression'])
    exact_reconstruction=sum((endpoints(v['integrated_rectangle'])[1] for v in out['time_slabs']),Q(0))+endpoints(out['infinite_tail'])[1]
    assert exact_reconstruction==rational and lo<=rational<simple
    outcome='FIXED_UPPER_BUDGET_EXCEEDS_SAVED_D_PLUS_F_ALLOWANCE' if lo>allowance else ('FIXED_UPPER_BUDGET_WITHIN_SAVED_ALLOWANCE' if simple<allowance else 'UNRESOLVED_COMPARISON')
    out['comparison']={'status':outcome,'saved_D_F_allowance':str(allowance),'old_short_upper_read_only':str(short),
        'new_same_domain_full_middle_upper':str(combined),'formula_lower_exceeds_allowance':lo>allowance,
        'actual_norm_lower_bound_claimed':False,'model_excluded':False,'new_C6_bound':False}
    done(9,name='exact_saved_endpoint_comparison_and_same_domain_short_addition',comparison=outcome)
    pins();out['status']='COMPLETE_VALID_MIDDLE_LONG_BUDGET_'+outcome
    out['registration_sha256']=sha(H/'registration.json');out['parameters_sha256']=sha(H/'parameters.json')
    out['source_sha256']=reg['source_sha256'];out['body32_kernel_calls']=0;out['global_I']=False;out['global_g']=False;out['new_C6_bound']=False
    receipt.update(status='COMPLETED_STOPPED',exit_code=0);record();receipt['results_sha256']=sha(H/'results.json');save('run_receipt.json',receipt)
    print(json.dumps({'status':out['status'],'configurations':receipt['completed_configurations'],'pid':os.getpid(),
                      'wall_seconds':receipt['elapsed_seconds'],'cpu_seconds':receipt['cpu_seconds'],'strict_long_upper':str(simple),
                      'long_display_only':str(total),'slab_selection_counts':chosen_counts,'infinite_tail_display_only':str(tail)}))
except BaseException as ex:
    receipt.update(status='FAILED_STOP_NO_RETRY',exit_code=1,error=repr(ex),traceback=traceback.format_exc());record();raise
finally:
    signal.setitimer(signal.ITIMER_REAL,0)
