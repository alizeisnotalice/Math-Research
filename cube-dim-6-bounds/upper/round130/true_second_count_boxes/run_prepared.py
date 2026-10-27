#!/usr/bin/env python3
"""One released130 true-second-count positive box pilot; no auto restart."""
import sys
sys.set_int_max_str_digits(0)
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[key]='1'
import json,time,resource,traceback
from pathlib import Path
from fractions import Fraction as Q
from box_core import atomic,compressed,sha,require,read_cones,Evaluator,grid
HERE=Path(__file__).resolve().parent

def main():
    start=time.monotonic();cpu=time.process_time();spec=json.loads((HERE/'inputs.json').read_text());reg=json.loads((HERE/'registration.json').read_text());rel=json.loads((HERE/'parent_release.json').read_text())
    require(rel.get('authorized')is True and str(rel.get('round'))=='130','sole primary release')
    require(rel['script_sha256']==sha(__file__) and rel['registration_sha256']==sha(HERE/'registration.json'),'released exact source/registration')
    for n,pin in reg['prepared_pins'].items():require(sha(HERE/n)==pin,'prepared pin '+n)
    for row in json.loads((HERE/'source_manifest.json').read_text())['sources']:require(sha(row['path'])==row['sha256'],'source pin '+row['role'])
    require(os.getpriority(os.PRIO_PROCESS,0)>=10,'science nice priority at least10')
    with(HERE/'started.json').open('x')as f:json.dump({'pid':os.getpid(),'release_sha256':sha(HERE/'parent_release.json'),'unix_time':time.time()},f);f.write('\n');f.flush();os.fsync(f.fileno())
    status='FAILED_NO_RETRY';phase='source_binding';failure=None;ev=None;counts=[];leaves=[];current=None;output={}
    try:
        require(spec['h']=='1/128'and spec['precision_bits']==512 and spec['leaf_cap_per_cone']==256 and spec['tail_L']==24,'registered fixed mechanism')
        old=json.loads(Path(spec['old_caps_path']).read_text());up128=json.loads(Path(spec['round128_input_path']).read_text())
        require(Q(old['h'])==Q(spec['h']) and Q(up128['h_values'][5])==Q(spec['h']),'same exact h')
        require(spec['p_interval']==up128['p_intervals'][5],'saved whole p interval')
        oldL,oldU=map(Q,spec['old_beta2_interval']);require(oldL==Q(old['caps']['b_lower'][2])and oldU==Q(old['caps']['beta_upper'][2])==Q(up128['beta_by_h'][5][1]),'complete old cap interval and exact128 beta2')
        cones,orbit=read_cones(spec);atomic(HERE/'cone_table.json',{'cones':cones,'pair_orbit_counts':orbit,'source_manifest_sha256':sha(HERE/'source_manifest.json')})
        ev=Evaluator(spec);tstar=1-Q(1,2**24);tailprob=1-tstar**6;tail=496*ev.pU**2*tailprob
        pi=ev.arb.pi();factorL,factorU=grid(ev.A(46080)/pi**3,384)
        def summarize(all_leaves):
            low=sum(Q(cones[j]['alpha'])*sum(Q(z['integral_interval'][0])for z in rows)for j,rows in enumerate(all_leaves))
            upper=sum(Q(cones[j]['alpha'])*(sum(Q(z['integral_interval'][1])for z in rows)+tail)for j,rows in enumerate(all_leaves))
            return {'pi_cubed_weighted_interval':[str(46080*low),str(46080*upper)],'moment_interval':[str(factorL*low),str(factorU*upper)],
                    'weighted_box_width':str(factorU*sum(Q(cones[j]['alpha'])*sum(Q(z['integral_interval'][1])-Q(z['integral_interval'][0])for z in rows)for j,rows in enumerate(all_leaves))),
                    'tail_upper':str(factorU*sum(Q(c['alpha'])for c in cones)*tail),'pi_factor_width':str((factorU-factorL)*low),'tail_probability':str(tailprob)}
        phase='first21_root_boxes'
        for cone in cones:
            current=cone['id'];ev.bind(cone);leaves.append([ev.box((Q(0),)*6,(tstar,)*6,'')]);counts.append(1)
        roots=summarize(leaves);trackedLow=Q(roots['moment_interval'][0]);require(trackedLow<=oldU,'CONTRADICTION: valid global lower exceeds accepted old beta2 upper');compressed(HERE/'first21_root_boxes.json.gz',{'cone_roots':leaves,'summary':roots,'kernel_calls':ev.kernel_calls})
        atomic(HERE/'root_feedback.json',{'status':'COMPLETE_FIRST21_ROOT_BOXES','evaluated_boxes':21,**roots,
              'old_beta2_interval':spec['old_beta2_interval'],'strict_stored128_cap_pass':Q(roots['moment_interval'][1])<oldU,
              'no_failure_inference_from_wide_upper':True,'kernel_calls':ev.kernel_calls})
        print(json.dumps({'phase':phase,'evaluated_boxes':21,'feedback_sha256':sha(HERE/'root_feedback.json'),'moment_interval':roots['moment_interval']}),flush=True)
        # Only a true global lower >=oldU can exclude any strict upper<oldU.
        # A large current box upper is never an impossibility test.
        excluded=Q(roots['moment_interval'][0])>=oldU
        phase='finite_refinement'
        finaldir=HERE/'leaves';finaldir.mkdir(exist_ok=False)
        for cone in cones:
            current=cone['id'];ev.bind(cone)
            while not excluded and len(leaves[current])<256:
                rows=leaves[current]
                chosen=min(range(len(rows)),key=lambda i:(-(Q(rows[i]['integral_interval'][1])-Q(rows[i]['integral_interval'][0])),rows[i]['prefix']))
                parent=rows.pop(chosen);lo=list(map(Q,parent['lower_t']));hi=list(map(Q,parent['upper_t']))
                axis=max(range(6),key=lambda j:(hi[j]-lo[j],-j));mid=(lo[axis]+hi[axis])/2
                lh=hi.copy();lh[axis]=mid;rl=lo.copy();rl[axis]=mid
                children=[ev.box(lo,lh,parent['prefix']+'0'),ev.box(rl,hi,parent['prefix']+'1')]
                trackedLow+=factorL*Q(cone['alpha'])*(sum(Q(z['integral_interval'][0])for z in children)-Q(parent['integral_interval'][0]))
                require(trackedLow<=oldU,'CONTRADICTION: valid global lower exceeds accepted old beta2 upper')
                rows.extend(children);counts[current]+=2
                require(counts[current]==2*len(rows)-1<=511,'finite per-cone box budget')
                if len(rows)%64==0:
                    compressed(HERE/'checkpoint_current_cone.json.gz',{'cone_id':current,'leaves':rows,'evaluated_boxes_by_cone':counts})
                    atomic(HERE/'progress.json',{'phase':phase,'cone_id':current,'leaves_current':len(rows),'evaluated_boxes':sum(counts),'kernel_calls':ev.kernel_calls})
            leaves[current].sort(key=lambda row:row['prefix']);p=finaldir/('cone_'+str(current).zfill(2)+'.json.gz')
            compressed(p,{'cone_id':current,'cone_table_sha256':sha(HERE/'cone_table.json'),'evaluated_boxes':counts[current],'leaves':leaves[current]})
            atomic(HERE/'checkpoint_completed.json',{'completed_cones':current+1,'evaluated_boxes_by_cone':counts,'kernel_calls':ev.kernel_calls,'last_cone_leaf_file_sha256':sha(p)})
            print(json.dumps({'phase':phase,'completed_cones':current+1,'leaves':len(leaves[current]),'evaluated_boxes':sum(counts),'kernel_calls':ev.kernel_calls}),flush=True)
        summary=summarize(leaves);muL,muU=map(Q,summary['moment_interval']);require(0<=muL<=muU,'entire positive moment interval');require(muL<=oldU,'CONTRADICTION: valid global lower exceeds accepted old beta2 upper');require(muL==trackedLow,'tracked exact global lower');require(muU-muL==Q(summary['weighted_box_width'])+Q(summary['tail_upper'])+Q(summary['pi_factor_width']),'complete exact box tail pi width decomposition')
        require(sum(counts)<=10731 and ev.kernel_calls<=686784,'registered global evaluated-box/kernel bound')
        status='COMPLETE_FIXED_FINITE_BOX_PILOT'
        output={'status':status,**summary,'old_beta2_interval':spec['old_beta2_interval'],'strict_saved128_beta2_improvement':muU<oldU,'strict_old_cap_expression_improvement':muU<oldL,
             'old128_upper_minus_new_upper':str(oldU-muU),'old_expression_lower_minus_new_upper':str(oldL-muU),
             'global_lower_excludes_strict_saved_cap_pass':muL>=oldU,'early_stop_on_global_lower':excluded,
             'leaf_counts':[len(v)for v in leaves],'evaluated_boxes_by_cone':counts,'evaluated_boxes':sum(counts),'kernel_calls':ev.kernel_calls,
             'kernel_cache_per_cone_maxsize':2048,'precision_bits':512,'h':'1/128','tail_L':24,'probability_renormalized':False,
             'pair_terms':496,'cone_count':21,'pi_factor_interval':[str(factorL),str(factorU)],'true_new_cap_is_only_upper':True,
             'new_C6_bound':False,'global_g0_proved':False,'family_failure_claim_from_nonpass':False,
             'leaf_files':[{'path':'leaves/'+p.name,'sha256':sha(p)}for p in sorted(finaldir.glob('*.json.gz'))]}
        atomic(HERE/'results.json',output)
    except BaseException as e:
        failure=type(e).__name__+': '+str(e);atomic(HERE/'failure.json',{'phase':phase,'cone_id':current,'failure':failure,'traceback':traceback.format_exc()});raise
    finally:
        atomic(HERE/'run_receipt.json',{'status':status,'phase':phase,'failure':failure,'pid':os.getpid(),'wall_before_receipt':time.monotonic()-start,'CPU_before_receipt':time.process_time()-cpu,
             'evaluated_boxes_by_cone':counts,'kernel_calls':None if ev is None else ev.kernel_calls,'children':0,'science_CPU_cores':1,'priority':os.getpriority(os.PRIO_PROCESS,0),
             'max_rss_bytes_macos':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'time_or_CPU_ceiling':None,'automatic_restarts':0,'optimizer_MC_native_MATLAB_calls':0})
    print(json.dumps({'status':status,'evaluated_boxes':sum(counts),'strict_saved128_beta2_improvement':output['strict_saved128_beta2_improvement'],'moment_interval':output['moment_interval']}),flush=True)
if __name__=='__main__':main()
