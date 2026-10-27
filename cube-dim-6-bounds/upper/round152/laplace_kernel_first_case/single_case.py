#!/usr/bin/env python3
"""Terra-only, dual-gated one real m2 Laplace-kernel Arb512 certificate."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,os,sys,time,traceback
sys.set_int_max_str_digits(0)
P=Path(__file__).resolve().parent
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def require(ok,msg):
    if not ok:raise ArithmeticError(msg)
def save(n,x):
    with(P/n).open('x')as f:json.dump(x,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def run():
    reg=read(P/'registration.json');inp=read(P/'frozen_sources.json');rel=read(P/'parent_release.json')
    require(rel.get('round')=='152'and rel.get('authorized')is True and rel.get('main_approved')is True and rel.get('independent_review_approved')is True,'EXACT_HASH_DUAL_GATE')
    for k,n in [('case_sha256','single_case.py'),('proof_sha256','PROOF.md'),('general_interface_sha256','GENERAL_2D_INTERFACE.md'),('frozen_sources_sha256','frozen_sources.json')]:require(reg[k]==sha(P/n)==rel.get(k),'FROZEN_'+n)
    require(rel.get('registration_sha256')==sha(P/'registration.json'),'REGISTRATION_HASH')
    require(not(P/'started.json').exists(),'ONE_SCIENTIFIC_CASE_NO_RESTART')
    for s in inp['sources']:require(sha(s['path'])==s['sha256'],'READONLY_SOURCE_HASH')
    require(inp['h']=='1/128'and inp['lambda']=='30'and inp['raw_pair_multiplicity']==1 and inp['is_full_beta2']is False,'SINGLE_REAL_FIRST_CASE')
    os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
    require(all(os.environ.get(k)=='1'for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']),'ONE_THREAD_ENV')
    save('started.json',{'round':'152','pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'workers':1,'release_sha256':sha(P/'parent_release.json'),'time_or_CPU_limit':None});start=time.monotonic();cpu=time.process_time()
    try:
        import flint
        from flint import arb,ctx
        require(flint.__version__=='0.9.0'and flint.__FLINT_VERSION__=='3.6.0','PINNED_FLINT')
        ctx.prec=512
        def ar(q):q=Q(q);return arb(q.numerator)/q.denominator
        def ends(v):
            require(v.is_finite(),'FINITE_ARB_INTERVAL');return Q(str(v.lower().fmpq())),Q(str(v.upper().fmpq()))
        def encode(v):
            a,b=ends(v);return {'lower_fraction':str(a),'upper_fraction':str(b),'ball_display_only':str(v)}
        counts={'sqrt_h':0,'erf':0,'erfc':0,'exp':0,'pi_constant':0,'closed_transform':0}
        h=ar(inp['h']);counts['sqrt_h']+=1;s=h.sqrt();counts['erf']+=1;p=s.erf()
        # This function includes analytic removable cases for static review;
        # the sole released call is lambda30, so 0/2 are not evaluated.
        def transform(lam):
            lam=Q(lam);counts['closed_transform']+=1
            if lam in (Q(0),Q(2)):
                counts['exp']+=1;counts['pi_constant']+=1;z=s*(-h).exp()/arb.pi().sqrt()
                return h*(1+p)-p/2+z if lam==0 else h*(1-p)+p/2-z
            counts['exp']+=1;counts['erfc']+=1
            numerator=(h*ar(lam*(lam-2))).exp()*(ar(lam-1)*s).erfc()-1+ar(lam-1)*p
            return numerator/ar(lam*(lam-2))
        I=transform(inp['lambda']);E=ar(inp['pair_factor'])*I**2;counts['pi_constant']+=1;pi=arb.pi();kappa=ar(inp['C130_factor_numerator'])/pi**3
        contribution=kappa*ar(inp['alpha'])*E
        iL,iU=ends(I);pL,pU=ends(p);eL,eU=ends(E);cL,cU=ends(contribution);gate=Q(1,10**80)
        require(0<iL<=iU<pL/30,'STRICT_POSITIVE_LAPLACE_AND_PROBABILITY_BOUND')
        require(0<eL<=eU<pL**2 and 0<cL<=cU,'TRUE_PAIR_PROBABILITY_AND_CONTRIBUTION')
        require(iU-iL<=gate and eU-eL<=gate and cU-cL<=gate,'FIXED_WIDTH_GATE_NO_PRECISION_RETRY')
        require(counts=={'sqrt_h':1,'erf':1,'erfc':1,'exp':1,'pi_constant':1,'closed_transform':1},'ONE_NEW_SCALAR_TRANSFORM_COUNTS')
        result={'status':'PASS_ONE_REAL_M2_CLOSED_LAPLACE_CERTIFICATE','h':'1/128','lambda':'30','cone_id':0,'body_pair_zero_based':[16,17],'alpha':'1/171992678400','p_h_interval':encode(p),'I30_interval':encode(I),'pair_expectation_interval':encode(E),'C130_factor_interval':encode(kappa),'partial_mu2_contribution_interval':encode(contribution),'raw_pair_multiplicity_applied':1,'analytic_density_mass':'1800/(2*30^2)=1','triangle_symmetry_identity':'E_pair=900*I30^2; no output-independence assumption','Arb_bits':512,'actual_scalar_calls':counts,'pointwise_k_calls':0,'numerical_quadrature_calls':0,'new_density_generations':0,'old_343456_functional_replays':0,'MC_native_calls':0,'is_full_beta2_upper':False,'new_M_upper':False,'global_g0_proved':False,'new_C6_bound':False,'case_sha256':sha(P/'single_case.py'),'release_sha256':sha(P/'parent_release.json'),'sources_sha256':sha(P/'frozen_sources.json'),'wall_before_receipt_seconds':time.monotonic()-start,'CPU_before_receipt_seconds':time.process_time()-cpu}
        save('certificate.json',result);save('exit_receipt.json',{'status':result['status'],'certificate_sha256':sha(P/'certificate.json'),'actual_scalar_calls':counts,'is_full_beta2_upper':False,'new_C6_bound':False})
        print(json.dumps({k:result[k]for k in ['status','I30_interval','pair_expectation_interval','partial_mu2_contribution_interval','is_full_beta2_upper']},indent=2))
    except BaseException:
        save('exit_receipt.json',{'status':'FAILED_STOP_NO_RETRY','failure':traceback.format_exc(),'is_full_beta2_upper':False,'new_C6_bound':False});raise
if __name__=='__main__':run()
