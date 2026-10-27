#!/usr/bin/env python3
"""One new joint-count optimization and one complete saved candidate."""
import sys
sys.set_int_max_str_digits(0)
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import gc,hashlib,json,resource,time
from pathlib import Path
from fractions import Fraction as Q
from datetime import datetime,timezone
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def req(x,s):
    if not x:raise ArithmeticError(s)
def utc():return datetime.now(timezone.utc).isoformat()
def atomic(n,d):
    p=HERE/n;t=p.with_name(p.name+'.own-tmp')
    with t.open('w')as f:
        json.dump(d,f,ensure_ascii=False,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(t,p)
def end_fraction(a,side):
    m,e=(a.lower()if side==0 else a.upper()).man_exp();return Q(int(m))*Q(2)**int(e)

def main():
    start,cpu=time.monotonic(),time.process_time()
    reg=json.loads((HERE/'registration.json').read_text());rel=json.loads((HERE/'parent_release.json').read_text())
    req(rel.get('authorized')is True and str(rel.get('round'))=='128','primary one new candidate release')
    req(rel['script_sha256']==sha(__file__) and rel['registration_sha256']==sha(HERE/'registration.json'),'released source and registration')
    for n,h in reg['prepared_pins'].items():req(sha(HERE/n)==h,'prepared pin '+n)
    sources=json.loads((HERE/'source_manifest.json').read_text())['sources']
    for s in sources:req(sha(s['path'])==s['sha256'],'source '+s['role'])
    with(HERE/'started.json').open('x')as f:
        json.dump({'pid':os.getpid(),'utc':utc(),'release_sha256':sha(HERE/'parent_release.json')},f);f.write('\n')
    status='FAILED';phase='new_optimization';failure=None;objective_calls=None;done_states=0
    try:
        inp=json.loads((HERE/'upstream.json').read_text());theta=Q(inp['threshold_M_power6_exact'])
        req(theta==Q(8,9)**6*24576**5/65000**6,'complete123 rational threshold')
        # Recheck exact provenance, permitted beta replacement and full rational warm start.
        src=json.loads((HERE/'source_manifest.json').read_text())
        previous=json.loads(Path(src['round125_input_path']).read_text())
        new126=json.loads(Path(src['round126_results_path']).read_text())
        previous_candidate=json.loads(Path(src['round125_candidate_path']).read_text())
        req(inp['warm_start_49']['numerators']==previous_candidate['numerators'] and inp['warm_start_49']['denominator']==previous_candidate['denominator'],'exact full125 warm simplex')
        for key in ('h_values','beta0','p_intervals','Phi_table','threshold_M_power6_exact','prefix_Kmax','prefix_state_count','tail_K_values','arithmetic_bits'):
            req(inp[key]==previous[key],'unchanged125 '+key)
        for j,row in enumerate(new126['rows']):
            req(Q(row['h'])==Q(inp['h_values'][j]),'same126 h')
            req(Q(inp['beta_by_h'][j][0])==Q(row['EdN_interval'][1]),'exact true126 beta1 upper')
            req(0<Q(row['EdN_interval'][0])<=Q(row['EdN_interval'][1])<Q(previous['beta_by_h'][j][0]),'strict126 cap improvement')
            req(inp['beta_by_h'][j][1:]==previous['beta_by_h'][j][1:],'unchanged125 beta2..6')

        from optimize_core import optimize
        candidate,opt=optimize(inp,lambda d:atomic('optimization_checkpoint.json',d))
        objective_calls=opt['objective_calls'];req(objective_calls<=1000,'single solver objective cap')
        atomic('candidate.json',candidate);atomic('optimization_receipt.json',opt)
        del optimize;gc.collect()
        phase='single_candidate_exact_joint_certificate'
        from joint_core import Joint
        joint=Joint(inp,candidate)
        def progress(count,chunks):
            nonlocal done_states
            done_states=count
            if chunks%10==0:print(json.dumps({'phase':phase,'prefix_states':count,'chunks':chunks}),flush=True)
        cert=joint.certify(HERE,opt['state_bytes_sha256'],atomic,progress)
        HL,HU=map(Q,cert['H_interval']);upper6=HU**5;prefix_lower6=HL**5
        from flint import arb,fmpq
        expression=arb(fmpq(upper6.numerator,upper6.denominator)).root(6)
        ml,mu=end_fraction(expression,0),end_fraction(expression,1)
        req(0<ml<=mu and ml**6<=upper6<=mu**6,'exact sixth-root expression enclosure')
        old=inp['124_exact_comparisons'];oldbest=Q(old['best_M_power6_upper']);oldmin=min(Q(v)for v in old['family_M_power6_lowers'])
        status='COMPLETE_ONE_JOINT_CANDIDATE'
        cert.update({'status':'STRICT_MOMENT_THRESHOLD_PASS'if upper6<theta else'UNSEPARATED_OR_NONCOMPETITIVE_FIXED_CANDIDATE',
            'threshold_M_power6_exact':str(theta),'M_power6_upper':str(upper6),'prefix_only_H_power5_lower':str(prefix_lower6),
            'threshold_minus_upper_exact':str(theta-upper6),'M_upper_expression_interval':[str(ml),str(mu)],
            'M_upper_expression_display':[format(float(ml),'.17g'),format(float(mu),'.17g')],
            'full_I_upper_power5_from123':str(Q(3,16)**6*upper6/384**5),
            'beats_saved124_best_upper':upper6<oldbest,'old_best_upper_minus_new_upper_exact':str(oldbest-upper6),
            'beats_all_eight_fixed_single_h_cap_families':upper6<oldmin,
            'old_minimum_family_lower_minus_new_upper_exact':str(oldmin-upper6),
            'beats_saved125_upper':upper6<Q(inp['125_exact_comparisons']['M_power6_upper']),
            'old125_upper_minus_new_upper_exact':str(Q(inp['125_exact_comparisons']['M_power6_upper'])-upper6),
            'actual_M_lower_bound':False,'joint_family_optimum_excluded':False,'global_g0_proved':False,'new_C6_bound':False})
        atomic('candidate_certificate.json',cert)
        atomic('results.json',{'status':status,'candidate_sha256':sha(HERE/'candidate.json'),
            'certificate_sha256':sha(HERE/'candidate_certificate.json'),'optimization_sha256':sha(HERE/'optimization_receipt.json'),
            'objective_calls':objective_calls,'prefix_states':735471,'tail_values':16,'new_candidates':1,
            'certificate_status':cert['status'],'M_upper_expression_display':cert['M_upper_expression_display'],
            'beats_saved124_best_upper':cert['beats_saved124_best_upper'],
            'beats_all_eight_fixed_single_h_cap_families':cert['beats_all_eight_fixed_single_h_cap_families'],
            'beats_saved125_upper':cert['beats_saved125_upper'],'true126_beta1_used':True,
            'no_tail_discard_or_renormalization':True,'global_g0_proved':False,'new_C6_bound':False})
        print(json.dumps({'status':status,'certificate_status':cert['status'],'M_upper':cert['M_upper_expression_display'],
                          'objective_calls':objective_calls,'prefix_states':735471}),flush=True)
    except BaseException as e:failure=type(e).__name__+': '+str(e);raise
    finally:
        atomic('run_receipt.json',{'status':status,'phase':phase,'failure':failure,'pid':os.getpid(),'utc':utc(),
            'objective_calls':objective_calls,'committed_prefix_states':done_states,'CPU_sample_seconds':time.process_time()-cpu,
            'wall_sample_before_receipt':time.monotonic()-start,'scope':'before final receipt and exit; observer covers full child',
            'max_rss_bytes_macos':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'science_CPU_cores':1,'children':0,
            'wall_limit':None,'CPU_limit':None,'old_program_or_optimizer_replays':0,'MC_or_native_calls':0})
if __name__=='__main__':main()
