#!/usr/bin/env python3
"""Four new fixed122 algebra/scalar configurations; no upstream execution."""
import sys
sys.set_int_max_str_digits(0)
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[key]='1'
import hashlib,json,resource,time
from pathlib import Path
from fractions import Fraction as Q
from datetime import datetime,timezone
HERE=Path(__file__).resolve().parent
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def req(x,s):
    if not x: raise ArithmeticError(s)
def utc(): return datetime.now(timezone.utc).isoformat()
def atomic(n,d):
    p=HERE/n;t=p.with_name(p.name+'.own-tmp')
    with t.open('w') as f:
        json.dump(d,f,ensure_ascii=False,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(t,p)
def dyad(m,e): return Q(int(m))*Q(2)**int(e)
def saved(r):
    return [dyad(r[s+'_mantissa'],r[s+'_binary_exponent']) for s in ('lower','upper')]
def ss(v): return [str(x) for x in v]
def disp(v): return [format(float(Q(x)),'.17g') for x in v]
def root5(x):
    den=2**192;lo=0;hi=den
    for _ in range(256):
        if Q(hi,den)**5>=x: break
        hi*=2
    else: raise ArithmeticError('root doubling cap')
    for _ in range(512):
        if hi-lo<=1: break
        m=(hi+lo)//2
        if Q(m,den)**5<=x: lo=m
        else: hi=m
    else: raise ArithmeticError('root bisection cap')
    a,b=Q(lo,den),Q(hi,den)
    req(a**5<=x<=b**5,'exact fifth root')
    return a,b


def cfg_peaks(inp):
    from flint import arb,fmpq,ctx
    ctx.prec=512
    rows=[]
    for c in inp['peaks033']:
        r=int(c['dimension']);D=Q(c['Gram_determinant']);R=Q(c['R_exact'])
        ratio2=Q(2+Q(2*r,3))**r/(16**(2*r)*D*R**2)
        req(ratio2>1,'BL weaker than exact033 '+c['id'])
        ratio=arb(fmpq(ratio2.numerator,ratio2.denominator)).sqrt()
        rows.append({'id':c['id'],'rank':r,'ratio_BL_over_exact_squared':str(ratio2),
                     'strict_residual':str(ratio2-1),'ratio_display_only':str(ratio),
                     'best_direct_peak':'unchanged exact033'})
    req(len(rows)==9,'all1+3+5 exact projection cases')
    triples=inp['triple_classes'];counts=[Q(29,32)*int(c['unordered_count']) for c in triples]
    req(all(n.denominator==1 for n in counts),'restricted exact counts')
    req([int(n) for n in counts]==[435,1740,435,580,1305],'all no-ones classes')
    req(sum(counts)*6==26970,'all ordered no-ones triples')
    crem=24576*sum(n*Q(c['R3']) for n,c in zip(counts,triples))
    req(crem==Q(44743079,2970)==Q(29,32)*Q(inp['B3_pi6_exact']), 'complete restored Crem identity')
    return {'id':'legal_full_projection_peak_comparisons','status':'ALL_DIRECT_RANK3_4_5_REPLACEMENTS_WEAKER',
            'nine_cases':rows,'Crem_times_pi6_exact':str(crem),'no_ones_unordered_counts':[int(n) for n in counts],
            'direct_C_E_H_recovery':'0','new_peak_or_tuple_computations':0,
            'restoration_required_before_projection_BL':True}


def cfg_min_integral(inp):
    from flint import arb,fmpq,ctx
    import sympy as sp
    ctx.prec=512
    def A(q):
        q=Q(q);return arb(fmpq(q.numerator,q.denominator))
    def ep(a,side):
        m,e=(a.lower() if side==0 else a.upper()).man_exp()
        return dyad(m,e)
    d=Q(31,64);sv=[6*d*d,Q(10,3)*d*d,2*d*d];av=[992,960,928]
    C0=Q(44743079,2970);p0=Q(216,2**30)
    def poly(u):
        z=Q(1)
        for s,a in zip(sv,av): z*=s+a*u
        return C0*u**3-p0*z
    req(C0>p0*av[0]*av[1]*av[2],'strictly stronger rank6 large-time cap')
    lo,hi=Q(0),Q(1)
    req(poly(lo)<0<poly(hi),'single fixed crossover bracket')
    for _ in range(192):
        mid=(lo+hi)/2
        if poly(mid)<=0: lo=mid
        else: hi=mid
    req(0<lo<hi and poly(lo)<=0<=poly(hi),'exact new crossover bracket')
    tl,tu=ep(A(lo).sqrt(),0),ep(A(hi).sqrt(),1)
    req(tl>0 and tl*tl<=lo and tu*tu>=hi,'exact crossover square roots')
    beta=[s/a for s,a in zip(sv,av)]
    req(len(set(beta))==3,'distinct positive partial-fraction roots')
    coeff=[]
    for i,b in enumerate(beta):
        z=Q(1)
        for j,bb in enumerate(beta):
            if j!=i: z*=bb-b
        coeff.append(1/z)
    x=sp.symbols('x',real=True)
    sym=lambda v:sp.Rational(v.numerator,v.denominator)
    residual=sp.factor(sum(sym(h)/(x*x+sym(b)) for h,b in zip(coeff,beta))
                       -1/sp.prod(x*x+sym(b) for b in beta))
    req(residual==0,'new finite primitive exact partial fractions')
    def primitive(t):
        total=arb(0)
        for h,b in zip(coeff,beta):
            rt=A(b).sqrt();total+=A(h)*(A(t)/rt).atan()/rt
        return total/(av[0]*av[1]*av[2])
    fl,fu=ep(primitive(tl),0),ep(primitive(tu),1)
    req(fl>0 and fl<=fu,'positive new finite primitive enclosure')
    piL,piU=saved(inp['high062']['constants']['pi'])
    lower=C0/piU**6*fl+(p0/piU**6)/(5*tu**5)
    upper=C0/piL**6*fu+(p0/piL**6)/(5*tl**5)
    req(0<lower<upper,'positive full min time integral')
    return {'id':'actual_high_B_subdensity_full6D_min_integral','status':'COMPLETE_NEW_ALL_TIME_MIN_INTEGRAL',
            'B_actual_le_p_before_enlargement':True,'B_bound_le_p_asserted':False,
            'squares':ss(sv),'time_coefficients':av,'C0':str(C0),'p0':str(p0),
            'u_star_interval':ss([lo,hi]),'u_polynomial_endpoint_signs':ss([poly(lo),poly(hi)]),
            'T_interval':ss([tl,tu]),'T_display':disp([tl,tu]),'root_bisections':192,
            'beta':ss(beta),'partial_fraction_coefficients':ss(coeff),'partial_fraction_residual':str(residual),
            'finite_primitive_enclosure':ss([fl,fu]),'B_new_time_integral':ss([lower,upper]),
            'B_new_time_integral_display':disp([lower,upper]),
            'old_Bfull_integral_recomputed':False,'P_or_volume_or_epsilon_recomputed':False}


def cfg_high_budget(inp,cap):
    h=inp['high062'];vp=saved(h['constants']['V_power_5_over6'])
    pshort,plong=saved(h['Pshort']),saved(h['Plong'])
    jl,ju=map(Q,cap['B_new_time_integral'])
    bl,bu=vp[0]*jl,vp[1]*ju
    nl,nu=pshort[0]+plong[0]+bl,pshort[1]+plong[1]+bu
    oldN=saved(h['union_norm']);oldB=saved(h['Bfull'])
    il,iu=root5(nl**6)[0],root5(nu**6)[1]
    oldI=Q(inp['ledger']['H_root']['upper'])
    improved=nu<oldN[0] and iu<oldI
    return {'id':'new_H_budget_with_frozen_P_and_volume','status':'STRICT_FULL_H_UPPER_IMPROVEMENT' if improved else 'UNSEPARATED_STOP_FIXED_BATCH',
            'Pshort_frozen':ss(pshort),'Plong_frozen':ss(plong),'V_power_5_over6_frozen':ss(vp),
            'B_new_norm':ss([bl,bu]),'B_new_norm_display':disp([bl,bu]),
            'old_Bfull_norm':ss(oldB),'old_total_norm':ss(oldN),
            'new_H_norm':ss([nl,nu]),'new_H_norm_display':disp([nl,nu]),
            'new_H_I':ss([il,iu]),'new_H_I_display':disp([il,iu]),
            'strict_norm_improvement_residual':str(oldN[0]-nu),
            'old_H_account_I_upper':str(oldI),'strict_I_improvement_residual':str(oldI-iu),
            'strict_improvement':improved,'old_P_gap_volume_epsilon_integrals_replayed':0}


def cfg_ledger(inp,newH):
    led=inp['ledger']
    for name in ('C_root','E_root','H_root'):
        row=led[name];l,u=Q(row['lower']),Q(row['upper']);rad=Q(row['radicand']);deg=int(row['degree'])
        req(l**deg<=rad<=u**deg,'saved root certificate '+name)
    ce=Q(led['C_root']['upper'])+Q(led['E_root']['upper'])
    oldH=Q(led['H_root']['upper']);old=Q(led['paid_CEH'])
    req(ce+oldH==old==Q(inp['paid_CEH_111']),'complete exact old CEH identity')
    useH=min(oldH,Q(newH['new_H_I'][1]));new=ce+useH;recovery=old-new
    req(recovery>=0,'nonnegative recovery')
    rL,rU=Q(inp['old_far_allocation']['lower']),Q(inp['old_far_allocation']['upper'])
    return {'id':'same_disjoint_CEH_exact_account_recovery','status':'STRICT_PARTIAL_CEH_RECOVERY' if recovery>0 else 'NO_RECOVERY_FIXED_BATCH_STOP',
            'C_I_upper_unchanged':str(led['C_root']['upper']),'E_I_upper_unchanged':str(led['E_root']['upper']),
            'old_CEH_I_upper':str(old),'new_H_I_upper_used':str(useH),'new_CEH_I_upper':str(new),
            'new_CEH_I_upper_display':format(float(new),'.17g'),
            'I_recovery_exact':str(recovery),'I_recovery_display':format(float(recovery),'.17g'),
            'old_far_allocation':ss([rL,rU]),'new_far_sufficient_allocation':ss([rL+recovery,rU+recovery]),
            'new_far_sufficient_allocation_display':disp([rL+recovery,rU+recovery]),
            'g0_globally_proved':False,'actual_far_paid':False,'global_S_closed':False,'new_C6_bound':False}


def main():
    start,cpu=time.monotonic(),time.process_time()
    reg=json.loads((HERE/'registration.json').read_text());rel=json.loads((HERE/'parent_release.json').read_text())
    req(rel.get('authorized') is True and str(rel.get('round'))=='122','explicit primary release')
    req(rel['script_sha256']==sha(__file__),'released science source')
    req(rel['registration_sha256']==sha(HERE/'registration.json'),'released registration')
    for n,h in reg['prepared_pins'].items(): req(sha(HERE/n)==h,'prepared pin '+n)
    sources=json.loads((HERE/'source_manifest.json').read_text())['sources']
    for s in sources:req(sha(s['path'])==s['sha256'],'source '+s['role'])
    with (HERE/'started.json').open('x') as f:
        json.dump({'pid':os.getpid(),'utc':utc(),'release_sha256':sha(HERE/'parent_release.json')},f);f.write('\n')
    done=[];status='FAILED';failure=None
    try:
        inp=json.loads((HERE/'upstream.json').read_text())
        for fun in (lambda:cfg_peaks(inp),lambda:cfg_min_integral(inp),
                    lambda:cfg_high_budget(inp,done[1]),lambda:cfg_ledger(inp,done[2])):
            done.append(fun());atomic('checkpoint.json',{'count':len(done),'configurations':done})
        status='COMPLETE_FOUR_FIXED_CONFIGURATIONS'
        atomic('results.json',{'status':status,'configurations':done,'new_C6_bound':False,'global_g_far_S_unpaid':True,'release_sha256':sha(HERE/'parent_release.json')})
        print(json.dumps({'status':status,'pid':os.getpid(),'H_norm':done[2]['new_H_norm_display'],'CEH':done[3]['new_CEH_I_upper_display'],'recovery':done[3]['I_recovery_display']}),flush=True)
    except BaseException as e:
        failure=type(e).__name__+': '+str(e);raise
    finally:
        atomic('run_receipt.json',{'status':status,'completed_configurations':len(done),'failure':failure,'pid':os.getpid(),'utc':utc(),
            'CPU_sample_seconds':time.process_time()-cpu,'wall_sample_before_receipt':time.monotonic()-start,
            'scope':'before final receipt and exit; outer observer separately covers complete child',
            'max_rss_bytes_macos':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'children':0,'science_CPU_cores':1,'wall_limit':None,'CPU_limit':None,'native_calls':0,'old_production_programs_executed':0})
if __name__=='__main__':main()
