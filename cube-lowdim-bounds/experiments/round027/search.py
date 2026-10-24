#!/usr/bin/env python3
"""One 51574x288 strengthened finite LP, epsilon=1/2, 120s; import performs no experiment."""
from __future__ import annotations
from time import monotonic,time,sleep
_IMPORT_START=monotonic()
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
import numpy as np
import scipy
import scipy.sparse as sp
from scipy.io import savemat
_IMPORT_SECONDS=monotonic()-_IMPORT_START
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent/'round06';LAST=ROOT.parent/'round07'
SECONDS,LP_SECONDS,EPSILON=120.,20.,.5
K_TOL,A_TOL,GAIN=1e-9,1e-8,5e-7
DEFINITION_SHA='b5d9f5de7259b5c430ef6f0cb6e7e2c75310a41e2731497739e877d2e54cefd5'
LAST_SOURCE_SHA='c2a906680a95468f78450e01a595f44f68030b98b983ae7136f16b8c21b9b823'
BASELINE_SHA='21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
MATLAB=Path('/Applications/MATLAB_R2026a.app/bin/matlab')


class StopBudget(Exception):pass
class Rejected(Exception):pass


def dump(path,data):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+'\n');temporary.replace(path)


def finite(label,value):
    if not np.all(np.isfinite(value)):raise Rejected(label+'_nonfinite')
    return value


def digest(path,guard):
    h=sha256()
    with path.open('rb') as stream:
        while True:
            guard();block=stream.read(4*1024*1024)
            if not block:break
            h.update(block)
    return h.hexdigest()


def main():
    if (ROOT/'result.json').exists() or (ROOT/'ipc.local').exists():
        raise SystemExit('Existing finite run evidence: no repeat.')
    start=monotonic();process=None;batch_log=None;record=None;ipc=ROOT/'ipc.local'
    counts={k:0 for k in ('cache_files_verified','primitive_formula_calls','new_direction_entries',
        'baseline_slack_data_entries_checked','residual_batches_started','residual_batches_completed',
        'residual_entries_started','residual_entries_completed','original_target_decisions',
        'strengthened_target_decisions','affine_matrix_rows_checked','LP_requested','LP_started',
        'LP_completed','batch_processes_started')}
    result={'status':'started','phase':'input_integrity','counts':counts,'checks':[],
        'candidate_exported':False,'epsilon':'1/2','cache_checks':[],
        'limits':{'seconds_internal':120,'LP_calls':1,'LP_soft_seconds':20,'owned_batch_processes':1,
        'residual_batches':2,'full_scalar_values':103148,'new_primitive_calls':0,'new_direction_columns':0},
        'scope':'Finite strengthened point model only; not a continuous margin proof.'}
    def save():result['elapsed_seconds']=monotonic()-start;dump(ROOT/'result.json',result)
    def guard():
        if monotonic()-start>=SECONDS:raise StopBudget('internal_120s_deadline')
    save()
    try:
        owned_group=os.getpgrp()
        if owned_group!=os.getpid():raise Rejected('Dedicated wrapper-owned coordinator group required.')
        result['coordinator_process']={'pid':os.getpid(),'pgid':owned_group,'sid':os.getsid(0)}
        path=ROOT/'dependency_round06_search.py'
        if not path.exists():path=BASE/'search.py'
        if digest(path,guard)!=DEFINITION_SHA:raise Rejected('Reviewed arithmetic source changed.')
        module_spec=importlib.util.spec_from_file_location('round08_cached_definitions',path)
        D=importlib.util.module_from_spec(module_spec);sys.modules[module_spec.name]=D;module_spec.loader.exec_module(D)
        prior=json.loads((BASE/'result.json').read_text());last=json.loads((LAST/'result.json').read_text())
        if prior['source_sha256']!=DEFINITION_SHA or last['source_sha256']!=LAST_SOURCE_SHA or not prior['baseline_residual']['passed']:
            raise Rejected('Prior frozen cache/source chain invalid.')
        point_path=BASE/'point_manifest.json'
        if digest(point_path,guard)!=prior['dependencies']['previous_point_manifest_sha256']:
            raise Rejected('Exact fixed point manifest changed.')
        document=json.loads(point_path.read_text());keys={k:[(Q(u),Q(v)) for u,v in document[k]] for k in ('K','A')}
        if (len(keys['K']),len(keys['A']))!=(26604,24970):raise Rejected('Fixed pool dimensions changed.')
        points={k:np.asarray([(float(u),float(v)) for u,v in ps]) for k,ps in keys.items()}
        targets={'K':np.asarray([1. if u+v<=1 else 0. for u,v in keys['K']]),'A':np.zeros(len(keys['A']))}
        for entry in prior['local_cache_files']:
            cache_path=BASE/entry['file']
            if cache_path.stat().st_size!=entry['bytes'] or digest(cache_path,guard)!=entry['sha256']:
                raise Rejected('Original primitive/baseline cache corrupt or absent.')
            counts['cache_files_verified']+=1;result['cache_checks'].append({'file':'round06/'+entry['file'],'sha256':entry['sha256']});save()
        manifest_path=LAST/'local_artifact_manifest.json';manifest=json.loads(manifest_path.read_text())
        for name in ('K_directions_288.local.npy','A_directions_288.local.npy','ipc.local/result_06.json'):
            entry=manifest['files'][name];cache_path=LAST/name
            if cache_path.stat().st_size!=entry['bytes'] or digest(cache_path,guard)!=entry['sha256']:
                raise Rejected('Post-terminal final-cache manifest mismatch.')
            counts['cache_files_verified']+=1;result['cache_checks'].append({'file':'round07/'+name,'sha256':entry['sha256']});save()
        reference=json.loads((LAST/'ipc.local/result_06.json').read_text())
        if reference['exitflag']!=1 or reference['columns']!=288:raise Rejected('Frozen LP6 reference provenance invalid.')
        fine=[tuple(ij) for ij in last['fine_classes']]
        if len(fine)!=96 or len(set(fine))!=96 or any(not 0<=j<=i<=80 for i,j in fine):
            raise Rejected('Frozen 96 fine classes changed.')
        direction_order=['relative_c','relative_S']+[f'coarse_{i}_{j}' for i,j in D.COARSE]+[f'fine_{i}_{j}' for i,j in fine]
        if last['direction_order']!=direction_order:raise Rejected('Frozen final column order mismatch.')
        geometry={};base={};matrices={};new_targets={};relative_rhs={}
        for kind in ('K','A'):
            cache=BASE/'cache.local';g=object.__new__(D.Geometry)
            g.kind,g.points,g.guard=kind,points[kind],guard
            g.B={'fine':[sp.load_npz(cache/f'{kind}_fine_B{axis}.npz') for axis in (0,1)]}
            g.L={'fine':[np.load(cache/f'{kind}_fine_L{axis}.npy',mmap_mode='r') for axis in (0,1)]} if kind=='A' else {}
            g.singular=np.load(cache/f'{kind}_singular.npy',mmap_mode='r')
            base[kind]=finite('baseline_cache',np.load(cache/f'{kind}_baseline.npy',mmap_mode='r'))
            matrices[kind]=finite('frozen_288_cache',np.load(LAST/f'{kind}_directions_288.local.npy',mmap_mode='r'))
            if matrices[kind].shape!=(len(keys[kind]),288) or base[kind].shape!=(len(keys[kind]),):
                raise Rejected('Fixed final cache shape changed.')
            if not np.array_equal(matrices[kind][:,:192],np.load(cache/f'{kind}_directions_192.npy',mmap_mode='r')):
                raise Rejected('Frozen final coarse columns differ from original verified cache.')
            if g.singular.shape!=(len(keys[kind]),) or any(b.shape!=(len(keys[kind]),161) for b in g.B['fine']):
                raise Rejected('Fine primitive/singular cache shape changed.')
            if kind=='A' and any(l.shape!=(len(keys[kind]),161) for l in g.L['fine']):raise Rejected('Fine L shape changed.')
            slack=finite('baseline_slack',base[kind]-targets[kind]);negative=np.flatnonzero(slack<0.)
            result.setdefault('baseline_slack',{})[kind]={'minimum':float(slack.min()),'negative_count':len(negative),
                'zero_count':int(np.count_nonzero(slack==0.)),'negative_indices':negative.tolist()}
            counts['baseline_slack_data_entries_checked']+=len(slack);save()
            if len(negative):raise Rejected('Negative baseline slack, even if tiny; no clamp/row deletion permitted.')
            new_targets[kind]=finite('strengthened_target',targets[kind]+EPSILON*slack)
            relative_rhs[kind]=finite('relative_strengthened_rhs',(1.-EPSILON)*(targets[kind]-base[kind]))
            if np.any(new_targets[kind]<targets[kind]):raise Rejected('Strengthened finite target rounded below original.')
            result['baseline_slack'][kind]['rhs_association_max_difference']=float(np.max(abs(base[kind]+relative_rhs[kind]-new_targets[kind])))
            geometry[kind]=g;guard()
        source=ROOT.parent/'round01/scaled_certificate.json'
        if digest(source,guard)!=BASELINE_SHA:raise Rejected('Proven K1 changed.')
        spec=json.loads(source.read_text());c1=Q(spec['c']);a1={tuple(map(int,k.split(','))):Q(v) for k,v in spec['a'].items()}
        if (Q(spec['alpha']),Q(spec['R']),Q(spec['h']),int(spec['N']))!=(D.ALPHA_Q,D.R_Q,D.H_Q,D.N):raise Rejected('Fixed analytic parameters changed.')
        refinement_path=ROOT.parent/'round05/refinement_identity.json'
        if digest(refinement_path,guard)!=prior['dependencies']['refinement_identity_sha256']:raise Rejected('Refinement evidence changed.')
        mask=D.refinement(refinement_path)
        direction_maps=[D.refined_direction(ij,mask) for ij in D.COARSE]+[{ij:Q(1)} for ij in fine]
        spline_mass=sum((D.H_Q**2*len(D.orbit(*ij))*v/2 for ij,v in a1.items()),Q(0))
        singular_half=D.ALPHA*D.RADIUS**(2-D.ALPHA)/(2-D.ALPHA)
        base_cost=float(c1)*singular_half+float(spline_mass)
        objective=np.asarray([float(c1)*singular_half,float(spline_mass)]+
            [float(D.COARSE_Q**2*len(D.orbit(*ij))/2) for ij in D.COARSE]+
            [float(D.H_Q**2*len(D.orbit(*ij))/2) for ij in fine])
        lower=np.asarray([-.2,-.5]+[-.5]*190+[-.25]*96);upper=-lower
        result['dependencies']={'definition_sha256':DEFINITION_SHA,'previous_source_sha256':LAST_SOURCE_SHA,
            'post_terminal_manifest_sha256':digest(manifest_path,guard),'baseline_sha256':BASELINE_SHA,
            'point_manifest_sha256':digest(point_path,guard),'LP6_reference_sha256':manifest['files']['ipc.local/result_06.json']['sha256']}
        result['frozen_direction_order']=direction_order
        result['baseline_reuse']='Verified input cache data, not another full baseline kernel batch.'
        save();guard()

        def expand(t):
            coeff={ij:(1.+t[1])*float(v) for ij,v in a1.items()}
            for value,direction in zip(t[2:],direction_maps):
                for ij,weight in direction.items():coeff[ij]=coeff.get(ij,0.)+float(value)*float(weight)
            return (1.+t[0])*float(c1),coeff

        def cost(c,coeff):
            return float(c)*singular_half+sum(float(D.H_Q**2*len(D.orbit(*ij))/2)*float(v) for ij,v in coeff.items())

        def location_difference(kind,full,other):
            difference=abs(full-other);index=int(difference.argmax())
            scale=max(abs(float(full[index])),abs(float(other[index])))
            return {'index':index,'point':[str(x) for x in keys[kind][index]],
                'maximum_absolute_difference':float(difference[index]),'full_value':float(full[index]),
                'other_value':float(other[index]),'scale_at_argmax':scale,
                'relative_difference_at_argmax':float(difference[index])/scale if scale else None,
                'zero_scale':scale==0.}

        def full_check(c,coeff,label,t=None,reference=None):
            guard();counts['residual_batches_started']+=1
            if counts['residual_batches_started']>2:raise Rejected('Two full-value batches maximum.')
            matrix=D.full_matrix(coeff);values={};check={'label':label,'status':'started'}
            result['checks'].append(check);save()
            for kind,g in geometry.items():
                guard();counts['residual_entries_started']+=len(keys[kind]);save()
                values[kind]=finite('complete_'+kind,float(c)*g.singular+g.spline(matrix))
                tolerance=K_TOL if kind=='K' else A_TOL
                original_violation=targets[kind]-values[kind]
                strong_violation=new_targets[kind]-values[kind]
                original_bad=np.flatnonzero(original_violation>tolerance)
                strong_bad=np.flatnonzero(strong_violation>tolerance)
                item={'minimum_value':float(values[kind].min()),'tolerance':tolerance,
                    'original_max_violation':float(original_violation.max()),'original_bad_indices':original_bad.tolist(),
                    'original_bad_count':len(original_bad),'strengthened_max_violation':float(strong_violation.max()),
                    'strengthened_bad_indices':strong_bad.tolist(),'strengthened_bad_count':len(strong_bad)}
                check[kind]=item
                counts['original_target_decisions']+=len(keys[kind]);counts['strengthened_target_decisions']+=len(keys[kind])
                if t is not None:
                    increment=finite('cached_affine_increment',matrices[kind]@t)
                    model_violation=relative_rhs[kind]-increment
                    item['raw_model_bad_count']=int(np.count_nonzero(model_violation>tolerance))
                    item['raw_model_max_violation']=float(model_violation.max())
                    affine=finite('cached_affine_value',base[kind]+increment)
                    item['affine_full_difference']=location_difference(kind,values[kind],affine)
                    counts['affine_matrix_rows_checked']+=len(keys[kind])
                else:
                    item['affine_comparison']='Coefficientwise rationalization need not lie in the 288-column affine family; full coefficients are checked.'
                if reference is not None:item['rounding_full_value_difference']=location_difference(kind,values[kind],reference[kind])
                counts['residual_entries_completed']+=len(keys[kind]);save();guard()
            check.update(status='completed',passed=all(check[k]['original_bad_count']==0 and
                check[k]['strengthened_bad_count']==0 and check[k].get('raw_model_bad_count',0)==0 for k in ('K','A')),
                expanded_float_half_cost=cost(c,coeff))
            if t is not None:
                check['affine_float_half_cost']=base_cost+float(objective@t)
                check['cost_difference']=abs(check['expanded_float_half_cost']-check['affine_float_half_cost'])
            np.savez(ROOT/(label+'_full_values.npz'),K=values['K'],A=values['A'])
            dump(ROOT/(label+'_coefficients.json'),{'c':str(c),'a':{f'{i},{j}':str(v) for (i,j),v in coeff.items()}})
            counts['residual_batches_completed']+=1;save();return values,check

        # Publish one immutable all-pool model before launching a single owned batch.
        result['phase']='save_full_strengthened_model';guard();ipc.mkdir()
        dump(ipc/'control.json',{'maximum_calls':1,'maximum_soft_seconds':20,
            'deadline_unix':time()+max(0.,SECONDS-(monotonic()-start))})
        raw=sp.csr_matrix(np.vstack((matrices['K'],matrices['A'])))
        rhs=np.r_[relative_rhs['K'],relative_rhs['A']]
        finite('all_pool_raw',raw.data);finite('all_pool_rhs',rhs)
        path=ipc/'model_01.mat';temporary=ipc/'model_01.mat.tmp'
        savemat(temporary,{'A_ub':-raw,'b_ub':-rhs[:,None],'c':objective[:,None],
            'lb':lower[:,None],'ub':upper[:,None],'base_cost':base_cost},appendmat=False,do_compression=False)
        temporary.replace(path);model_hash=digest(path,guard)
        sp.save_npz(ipc/'model_01_raw.npz',raw,compressed=False)
        np.savez(ipc/'model_01_targets.npz',K_original=targets['K'],A_original=targets['A'],
            K_new=new_targets['K'],A_new=new_targets['A'],K_baseline=base['K'],A_baseline=base['A'],
            K_relative_rhs=relative_rhs['K'],A_relative_rhs=relative_rhs['A'])
        dump(ipc/'model_01_keys.json',{'rows':[[kind,i] for kind in ('K','A') for i in range(len(keys[kind]))],
            'columns':direction_order,'epsilon':'1/2','MATLAB_internal_model_not_read_back':True})
        result['model']={'rows':51574,'columns':288,'input_sha256':model_hash,
            'raw_sha256':digest(ipc/'model_01_raw.npz',guard),
            'targets_sha256':digest(ipc/'model_01_targets.npz',guard),
            'keys_sha256':digest(ipc/'model_01_keys.json',guard)}
        if not MATLAB.is_file():raise Rejected('Approved MATLAB batch binary unavailable; no fallback.')
        script=ROOT/'batch_solver.m';result['batch_source_sha256_at_launch']=digest(script,guard)
        code="run('"+str(script).replace("'","''")+"')";batch_log=(ROOT/'matlab_batch.log').open('x')
        guard();process=subprocess.Popen([str(MATLAB),'-batch',code],cwd=ROOT,stdin=subprocess.DEVNULL,
            stdout=batch_log,stderr=subprocess.STDOUT,start_new_session=False)
        counts['batch_processes_started']+=1;child_group=os.getpgid(process.pid)
        result['batch_process']={'pid':process.pid,'pgid':child_group,'coordinator_pgid':owned_group,
            'same_owned_group':child_group==owned_group,'shared_GUI_used':False,'started_seconds':monotonic()-start}
        save()
        if child_group!=owned_group:raise Rejected('Owned MATLAB child process group mismatch.')
        record={'iteration':1,'status':'queued','input_sha256':model_hash,'rows':51574,'columns':288,
            'basis_warm_start':False,'requested_seconds':monotonic()-start}
        result['solve']=record;counts['LP_requested']+=1;save()
        soft=min(20.,max(0.,SECONDS-(monotonic()-start)))
        if soft<=0:raise StopBudget('No remaining time for the one LP.')
        dump(ipc/'request_01.ready.json',{'iteration':1,'model_file':'model_01.mat',
            'input_sha256':model_hash,'soft_seconds':soft})

        def observe_started():
            marker=ipc/'started_01.json'
            if marker.exists() and not record.get('worker_started_observed'):
                state=json.loads(marker.read_text())
                if state['iteration']!=1 or state['input_sha256']!=model_hash:raise Rejected('Started model binding mismatch.')
                record.update(worker_started_observed=True,worker_started=state);counts['LP_started']+=1;save()

        ready=ipc/'result_01.ready.json'
        while not ready.exists():
            guard();observe_started()
            if (ipc/'worker_error.json').exists():
                record['worker_error']=json.loads((ipc/'worker_error.json').read_text())
                raise Rejected('Single batch worker failed; no second batch.')
            if process.poll() is not None:raise Rejected('Batch exited before ready result.')
            sleep(min(.05,max(.001,SECONDS-(monotonic()-start))))
        observe_started();marker=json.loads(ready.read_text());output_path=ipc/'result_01.json'
        if (marker['iteration'],marker['input_sha256'],marker['result_file'])!=(1,model_hash,output_path.name):
            raise Rejected('Ready result bound to wrong model.')
        if digest(output_path,guard)!=marker['result_sha256']:raise Rejected('Result SHA mismatch.')
        output=json.loads(output_path.read_text())
        if output['iteration']!=1 or output['input_sha256']!=model_hash or not record.get('worker_started_observed'):
            raise Rejected('Returned result/started ledger binding mismatch.')
        record.update(status='completed',worker_result=output,result_sha256=marker['result_sha256'])
        counts['LP_completed']+=1;save();guard()
        if output.get('status')!='returned' or output.get('exitflag')!=1 or not output.get('finite_primal'):
            raise Rejected('The only LP is not optimal/finite; no retry or candidate.')
        t=finite('returned_vector',np.asarray(output['x'],dtype=float).reshape(-1))
        if t.shape!=(288,):raise Rejected('Returned dimension mismatch.')
        bound_error=max(0.,float(np.max(lower-t)),float(np.max(t-upper)))
        record['variable_bound_error']=bound_error
        if bound_error>1e-10:raise Rejected('Fixed variable box violated.')
        result['phase']='float_full_pool';c,coeff=expand(t)
        values,check=full_check(c,coeff,'floating',t)
        if not check['passed']:raise Rejected('Original/strengthened/raw/full finite gate failed; no rhs repair.')
        if check['expanded_float_half_cost']<base_cost-GAIN:
            den=2**50;cq=Q(round(c*den),den);aq={ij:Q(round(v*den),den) for ij,v in coeff.items()}
            if cq<=0:raise Rejected('Rationalized c nonpositive.')
            result['phase']='unique_rationalized_full_pool'
            _,rational=full_check(cq,aq,'rationalized',reference=values)
            if not rational['passed']:raise Rejected('Unique rationalized gate failed; no repair.')
            nz={ij:v for ij,v in aq.items() if v};actual_U=(max(i for i,j in nz)+2)*D.H_Q
            if not D.R_Q<=actual_U<=Q(41,16):raise Rejected('Rationalized support scan bound invalid.')
            target=ROOT/'candidate.json'
            if target.exists():raise Rejected('Second candidate export forbidden.')
            guard();dump(target,{'alpha':str(D.ALPHA_Q),'R':str(D.R_Q),'h':str(D.H_Q),'N':D.N,
                'c':str(cq),'a':{f'{i},{j}':str(v) for (i,j),v in sorted(nz.items())}})
            result.update(candidate_exported=True,candidate_sha256=digest(target,guard),
                candidate_float_half_cost=cost(cq,aq),support_U=str(actual_U),rationalization_denominator=str(den),
                next='Lead schedules exact cost and both planned full-domain checks; no theorem claim yet.')
        else:result['termination']='no_significant_floating_gain'
        result['status']='single_strengthened_all_pool_LP_closed'
    except StopBudget as exc:result.update(status='interrupted_saved_checkpoint',reason=str(exc))
    except Exception as exc:result.update(status='failed_or_rejected',reason=f'{type(exc).__name__}: {exc}')
    finally:
        if ipc.exists() and not (ipc/'stop.json').exists():dump(ipc/'stop.json',{'stop':True,'reason':'coordinator_finally'})
        if process is not None:
            cleanup={'owned_pid':process.pid,'stop_requested':True,'started_seconds':monotonic()-start}
            try:process.wait(timeout=min(5.,max(.05,SECONDS-(monotonic()-start))))
            except subprocess.TimeoutExpired:
                cleanup['terminate_sent']=True;process.terminate()
                try:process.wait(timeout=min(2.,max(.05,SECONDS-(monotonic()-start))))
                except subprocess.TimeoutExpired:cleanup['kill_sent']=True;process.kill();process.wait(timeout=1.)
            cleanup.update(returncode=process.returncode,completed_seconds=monotonic()-start,
                leftover_cleanup='Outer wrapper cleans only its owned group; shared GUI untouched.')
            result['batch_cleanup']=cleanup
        if batch_log is not None:batch_log.close()
        if record is not None:
            try:
                start_path=ipc/'started_01.json'
                if start_path.exists() and not record.get('worker_started_observed'):
                    state=json.loads(start_path.read_text())
                    if state['iteration']!=1 or state['input_sha256']!=record['input_sha256']:raise Rejected('Cleanup start mismatch.')
                    record.update(worker_started_observed=True,worker_started=state);counts['LP_started']+=1
                if (ipc/'result_01.ready.json').exists() and record['status']!='completed':
                    marker=json.loads((ipc/'result_01.ready.json').read_text());p=ipc/'result_01.json'
                    if (marker['iteration'],marker['input_sha256'],marker['result_file'])!=(1,record['input_sha256'],p.name):raise Rejected('Cleanup result mismatch.')
                    if sha256(p.read_bytes()).hexdigest()!=marker['result_sha256'] or not record.get('worker_started_observed'):raise Rejected('Cleanup result hash/start mismatch.')
                    record.update(status='completed_during_cleanup_not_used',worker_result=json.loads(p.read_text()),result_sha256=marker['result_sha256']);counts['LP_completed']+=1
                elif record['status']!='completed':record['status']='interrupted_in_flight' if record.get('worker_started_observed') else 'requested_not_started'
            except Exception as exc:record['ledger_reconciliation_error']=f'{type(exc).__name__}: {exc}'
    assert counts['LP_completed']<=counts['LP_started']<=counts['LP_requested']<=1 and counts['batch_processes_started']<=1
    assert counts['residual_batches_started']<=2 and counts['residual_entries_started']<=103148
    assert counts['primitive_formula_calls']==counts['new_direction_entries']==0
    result.update(source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        batch_source_sha256_final=sha256((ROOT/'batch_solver.m').read_bytes()).hexdigest(),
        module_import_seconds_separate=_IMPORT_SECONDS,environment={'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__})
    save();print(json.dumps({'status':result['status'],'counts':counts,'seconds':result['elapsed_seconds'],
                           'candidate_exported':result['candidate_exported']},ensure_ascii=False),flush=True)


if __name__=='__main__':main()
