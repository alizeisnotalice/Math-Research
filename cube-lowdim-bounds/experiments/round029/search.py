#!/usr/bin/env python3
"""Three fixed-beta conservative 51574x128 LPs; <=180s/6 full batches.

Import performs no experiment. Old B/L/alpha-generator values are float caches,
not outward interval bounds. No beta hypergeometric formula is called here.
"""
from __future__ import annotations
from time import monotonic, time, sleep
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
SECONDS,LP_SECONDS,EPSILON=180.,20.,.5
BETAS=(Q(47,40),Q(6,5),Q(5,4))
K_TOL,A_TOL,GAIN=1e-9,1e-8,5e-7
DEFINITION_SHA='b5d9f5de7259b5c430ef6f0cb6e7e2c75310a41e2731497739e877d2e54cefd5'
LAST_SOURCE_SHA='c2a906680a95468f78450e01a595f44f68030b98b983ae7136f16b8c21b9b823'
BASELINE_SHA='21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
MATLAB=Path('/Applications/MATLAB_R2026a.app/bin/matlab')


class StopBudget(Exception):pass
class Rejected(Exception):pass


def dump(path,data):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    temporary.replace(path)


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
        raise SystemExit('Existing run evidence: no automatic repetition.')
    start=monotonic();process=None;batch_log=None;ipc=ROOT/'ipc.local'
    counts={k:0 for k in ('cache_files_verified','primitive_formula_calls','new_T_columns',
        'transformed_singular_entries_started','transformed_singular_entries_completed',
        'baseline_slack_data_entries_checked','baseline_model_rows_checked',
        'residual_batches_started','residual_batches_completed','residual_entries_started',
        'residual_entries_completed','original_target_decisions','strengthened_target_decisions',
        'affine_matrix_rows_checked','LP_requested','LP_started','LP_completed',
        'MATLAB_matrix_rows_checked','batch_processes_started','rationalizations_started')}
    result={'status':'started','phase':'input_integrity','counts':counts,'pairs':[],
        'checks':[],'cache_checks':[],'candidate_exported':False,'epsilon':'1/2',
        'limits':{'wall_seconds_hard_wrapper':180,'LP_calls':3,'LP_soft_seconds':20,
            'batch_processes':1,'full_batches':6,'full_scalar_values':309444,
            'new_primitive_calls':0,'new_T_columns':0,'transformed_singular_entries':154722},
        'scope':'Float conservative point model; A_beta is not directly evaluated or certified.'}
    def save():result['elapsed_seconds']=monotonic()-start;dump(ROOT/'result.json',result)
    def guard():
        if monotonic()-start>=SECONDS:raise StopBudget('internal_180s_deadline')
    def verify(path,expected,size=None):
        if size is not None and path.stat().st_size!=size:raise Rejected('Input byte count mismatch: '+path.name)
        actual=digest(path,guard)
        if actual!=expected:raise Rejected('Input SHA mismatch: '+path.name)
        counts['cache_files_verified']+=1
        result['cache_checks'].append({'file':str(path.relative_to(ROOT.parent)),'sha256':actual});save()
    save()
    try:
        owned_group=os.getpgrp()
        if owned_group!=os.getpid():raise Rejected('Dedicated wrapper-owned coordinator group required.')
        result['coordinator_process']={'pid':os.getpid(),'pgid':owned_group,'sid':os.getsid(0)}
        source_path=ROOT/'dependency_round06_search.py'
        if not source_path.exists():source_path=BASE/'search.py'
        verify(source_path,DEFINITION_SHA)
        module_spec=importlib.util.spec_from_file_location('round10_reviewed_cache_definitions',source_path)
        D=importlib.util.module_from_spec(module_spec);sys.modules[module_spec.name]=D;module_spec.loader.exec_module(D)
        prior=json.loads((BASE/'result.json').read_text());last=json.loads((LAST/'result.json').read_text())
        if prior['source_sha256']!=DEFINITION_SHA or last['source_sha256']!=LAST_SOURCE_SHA or not prior['baseline_residual']['passed']:
            raise Rejected('Old cache/source provenance invalid.')
        point_path=BASE/'point_manifest.json'
        verify(point_path,prior['dependencies']['previous_point_manifest_sha256'])
        document=json.loads(point_path.read_text())
        keys={k:[(Q(u),Q(v)) for u,v in document[k]] for k in ('K','A')}
        if (len(keys['K']),len(keys['A']))!=(26604,24970):raise Rejected('Fixed pool counts changed.')
        if any(not (0<=v<=u and u+v>0) for u,v in keys['K']):raise Rejected('K point domain changed.')
        if any(not (0<v<=u and u+v!=D.R_Q) for u,v in keys['A']):raise Rejected('A point domain changed.')
        points={k:np.asarray([(float(u),float(v)) for u,v in p]) for k,p in keys.items()}
        targets={'K':np.asarray([1. if u+v<=1 else 0. for u,v in keys['K']]),'A':np.zeros(len(keys['A']))}
        for entry in prior['local_cache_files']:verify(BASE/entry['file'],entry['sha256'],entry['bytes'])
        manifest_path=LAST/'local_artifact_manifest.json';manifest=json.loads(manifest_path.read_text())
        for name in ('K_directions_288.local.npy','A_directions_288.local.npy'):
            entry=manifest['files'][name];verify(LAST/name,entry['sha256'],entry['bytes'])
        baseline_path=ROOT.parent/'round01/scaled_certificate.json';verify(baseline_path,BASELINE_SHA)
        spec=json.loads(baseline_path.read_text());c1=Q(spec['c'])
        a1={tuple(map(int,k.split(','))):Q(v) for k,v in spec['a'].items()}
        if (Q(spec['alpha']),Q(spec['R']),Q(spec['h']),int(spec['N']))!=(D.ALPHA_Q,D.R_Q,D.H_Q,D.N):
            raise Rejected('Fixed alpha/R/h/N changed.')
        arithmetic_path=ROOT.parent/'round01/scaling_arithmetic.json'
        arithmetic=json.loads(arithmetic_path.read_text())
        if arithmetic['new_certificate_sha256']!=BASELINE_SHA:raise Rejected('B1 exact cost provenance changed.')
        b1_cost=float(Q(arithmetic['new_cost_lower']))
        fine=[tuple(ij) for ij in last['fine_classes']]
        coarse=[ij for ij in D.COARSE if sum(ij)<=9]
        if len(coarse)!=30 or len(fine)!=96 or len(set(fine))!=96 or any(not 0<=j<=i or i+j>48 for i,j in fine):
            raise Rejected('Frozen inner-support direction set invalid.')
        old_order=['relative_c','relative_S']+[f'coarse_{i}_{j}' for i,j in D.COARSE]+[f'fine_{i}_{j}' for i,j in fine]
        if last['direction_order']!=old_order:raise Rejected('Frozen 288-column order changed.')
        selected=[1]+[2+D.COARSE.index(ij) for ij in coarse]+list(range(192,288))
        order=['absolute_c','absolute_s']+[f'coarse_{i}_{j}' for i,j in coarse]+[f'fine_{i}_{j}' for i,j in fine]
        refinement_path=ROOT.parent/'round05/refinement_identity.json'
        verify(refinement_path,prior['dependencies']['refinement_identity_sha256'])
        mask=D.refinement(refinement_path)
        direction_maps=[D.refined_direction(ij,mask) for ij in coarse]+[{ij:Q(1)} for ij in fine]
        if any(i+j>48 for direction in direction_maps for (i,j),v in direction.items() if v):
            raise Rejected('T refinement escaped strict inner support.')
        for ij,direction in zip(coarse,direction_maps):
            mass=sum((D.H_Q**2*len(D.orbit(*kl))*v/2 for kl,v in direction.items()),Q(0))
            if mass!=D.COARSE_Q**2*len(D.orbit(*ij))/2:raise Rejected('Exact coarse half mass mismatch.')
        spline_mass=sum((D.H_Q**2*len(D.orbit(*ij))*v/2 for ij,v in a1.items()),Q(0))
        mass_tail=[float(spline_mass)]+[float(D.COARSE_Q**2*len(D.orbit(*ij))/2) for ij in coarse]+[
            float(D.H_Q**2*len(D.orbit(*ij))/2) for ij in fine]
        geometry={};base={};spline_columns={};new_targets={};radii={}
        for kind in ('K','A'):
            cache=BASE/'cache.local';g=object.__new__(D.Geometry)
            g.kind,g.points,g.guard=kind,points[kind],guard
            g.B={'fine':[sp.load_npz(cache/f'{kind}_fine_B{axis}.npz') for axis in (0,1)]}
            g.L={'fine':[np.load(cache/f'{kind}_fine_L{axis}.npy',mmap_mode='r') for axis in (0,1)]} if kind=='A' else {}
            g.singular=finite('cached_alpha_singular',np.load(cache/f'{kind}_singular.npy',mmap_mode='r'))
            base[kind]=finite('cached_baseline',np.load(cache/f'{kind}_baseline.npy',mmap_mode='r'))
            matrix=finite('frozen_288',np.load(LAST/f'{kind}_directions_288.local.npy',mmap_mode='r'))
            if matrix.shape!=(len(keys[kind]),288):
                raise Rejected('Cached fixed shape invalid.')
            if base[kind].shape!=(len(keys[kind]),) or g.singular.shape!=(len(keys[kind]),):raise Rejected('Cached scalar shape invalid.')
            if not np.array_equal(matrix[:,:192],np.load(cache/f'{kind}_directions_192.npy',mmap_mode='r')):
                raise Rejected('Frozen first 192 directions changed.')
            if any(b.shape!=(len(keys[kind]),161) for b in g.B['fine']):raise Rejected('B cache shape invalid.')
            if kind=='A' and any(l.shape!=(len(keys[kind]),161) for l in g.L['fine']):raise Rejected('L cache shape invalid.')
            spline_columns[kind]=finite('selected_S1_T_columns',np.asarray(matrix[:,selected]))
            slack=finite('old_baseline_slack',base[kind]-targets[kind]);negative=np.flatnonzero(slack<0.)
            result.setdefault('old_baseline_slack',{})[kind]={'minimum':float(slack.min()),
                'negative_indices':negative.tolist(),'zero_count':int(np.count_nonzero(slack==0.))}
            counts['baseline_slack_data_entries_checked']+=len(slack);save()
            if len(negative):raise Rejected('Negative baseline slack; no clamp or deletion.')
            new_targets[kind]=finite('strengthened_target',targets[kind]+EPSILON*slack)
            if np.any(new_targets[kind]<targets[kind]):raise Rejected('Strengthened target rounded below original.')
            radii[kind]=points[kind].sum(axis=1);geometry[kind]=g;guard()
        result['dependencies']={'definitions_sha256':DEFINITION_SHA,'last_source_sha256':LAST_SOURCE_SHA,
            'post_terminal_manifest_sha256':digest(manifest_path,guard),'baseline_sha256':BASELINE_SHA,
            'point_manifest_sha256':digest(point_path,guard),'refinement_sha256':digest(refinement_path,guard),
            'B1_arithmetic_sha256':digest(arithmetic_path,guard)}
        result['direction_order']=order;result['B1_float_reference']=b1_cost
        result['method']='chi(r)*A_alpha_psi_alpha lower model; not A(chi*psi), not the true beta generator value'
        save();guard()

        def expand(parameters,exact=False):
            c,s=parameters[:2];T={}
            for value,direction in zip(parameters[2:],direction_maps):
                for ij,weight in direction.items():
                    w=weight if exact else float(weight)
                    T[ij]=T.get(ij,Q(0) if exact else 0.)+value*w
            coeff={ij:s*(v if exact else float(v)) for ij,v in a1.items()}
            for ij,v in T.items():coeff[ij]=coeff.get(ij,Q(0) if exact else 0.)+v
            return c,s,{ij:v for ij,v in T.items() if v},{ij:v for ij,v in coeff.items() if v}

        def location(kind,amount,**arrays):
            index=int(np.argmax(amount))
            return {'index':index,'point':[str(x) for x in keys[kind][index]],'amount':float(amount[index]),
                    **{name:float(value[index]) for name,value in arrays.items()}}

        def full_check(parameters,beta,singular_columns,matrices,objective,label,pair):
            guard();counts['residual_batches_started']+=1
            if counts['residual_batches_started']>6:raise Rejected('Six full batches cap.')
            c,s,T,coeff=expand(parameters,exact=isinstance(parameters[0],Q))
            matrix=D.full_matrix(coeff);values={};check={'label':label,'beta':str(beta),'status':'started'}
            result['checks'].append(check);pair.setdefault('checks',[]).append(label);save()
            vector=np.asarray([float(x) for x in parameters]);finite('checked_parameters',vector)
            for kind,g in geometry.items():
                guard();counts['residual_entries_started']+=len(keys[kind]);save()
                values[kind]=finite('complete_conservative_'+kind,float(c)*singular_columns[kind]+g.spline(matrix))
                tolerance=K_TOL if kind=='K' else A_TOL
                original=targets[kind]-values[kind];strong=new_targets[kind]-values[kind]
                affine=finite('cached_affine_model',matrices[kind]@vector)
                raw_violation=new_targets[kind]-affine
                difference=abs(values[kind]-affine)
                item={'tolerance':tolerance,'minimum':float(values[kind].min()),
                    'original_bad_count':int(np.count_nonzero(original>tolerance)),
                    'strengthened_bad_count':int(np.count_nonzero(strong>tolerance)),
                    'raw_bad_count':int(np.count_nonzero(raw_violation>tolerance)),
                    'original_worst':location(kind,original,full=values[kind],target=targets[kind]),
                    'strengthened_worst':location(kind,strong,full=values[kind],target=new_targets[kind]),
                    'raw_worst':location(kind,raw_violation,affine=affine,target=new_targets[kind]),
                    'affine_full_difference':location(kind,difference,full=values[kind],affine=affine)}
                index=item['affine_full_difference']['index'];scale=max(abs(values[kind][index]),abs(affine[index]))
                item['affine_full_difference']['scale_at_argmax']=float(scale)
                item['affine_full_difference']['relative_at_argmax']=float(difference[index]/scale) if scale else None
                check[kind]=item;counts['affine_matrix_rows_checked']+=len(keys[kind])
                counts['original_target_decisions']+=len(keys[kind]);counts['strengthened_target_decisions']+=len(keys[kind])
                counts['residual_entries_completed']+=len(keys[kind]);save();guard()
            half=float(beta)*D.RADIUS**(2-float(beta))/(2-float(beta))
            full_cost=float(c)*half+sum(float(D.H_Q**2*len(D.orbit(*ij))/2)*float(v) for ij,v in coeff.items())
            check.update(status='completed',passed=all(check[k][name]==0 for k in ('K','A') for name in (
                'original_bad_count','strengthened_bad_count','raw_bad_count')),expanded_float_half_cost=full_cost,
                affine_float_half_cost=float(objective@vector),cost_difference=abs(full_cost-float(objective@vector)))
            np.savez(ROOT/(label+'_full_values.npz'),K=values['K'],A_lower_model=values['A'])
            dump(ROOT/(label+'_parameters.json'),{'beta':str(beta),'parameters':[str(x) for x in parameters],
                'c':str(c),'s':str(s),'T':{f'{i},{j}':str(v) for (i,j),v in sorted(T.items())},
                'a':{f'{i},{j}':str(v) for (i,j),v in sorted(coeff.items())}})
            counts['residual_batches_completed']+=1;save();return check,(c,s,T,coeff)

        def observe(pair,cleanup=False):
            iteration=pair['iteration'];prefix=f'{iteration:02d}'
            path=ipc/f'started_{prefix}.json'
            if path.exists() and not pair.get('worker_started_observed'):
                state=json.loads(path.read_text())
                if (state['iteration'],state['input_sha256'])!=(iteration,pair['input_sha256']):raise Rejected('Started binding mismatch.')
                pair.update(worker_started_observed=True,worker_started=state);counts['LP_started']+=1
            ready=ipc/f'result_{prefix}.ready.json'
            if ready.exists() and not pair.get('result_observed'):
                marker=json.loads(ready.read_text());path=ipc/f'result_{prefix}.json'
                if (marker['iteration'],marker['input_sha256'],marker['result_file'])!=(iteration,pair['input_sha256'],path.name):
                    raise Rejected('Result-ready binding mismatch.')
                if sha256(path.read_bytes()).hexdigest()!=marker['result_sha256']:raise Rejected('Result SHA mismatch.')
                output=json.loads(path.read_text())
                if (output['iteration'],output['input_sha256'])!=(iteration,pair['input_sha256']) or not pair.get('worker_started_observed'):
                    raise Rejected('Result body/started binding mismatch.')
                pair.update(result_observed=True,worker_result=output,result_sha256=marker['result_sha256'],
                    status='returned_not_used_cleanup' if cleanup else 'returned')
                counts['LP_completed']+=1
                counts['MATLAB_matrix_rows_checked']+=int(output.get('existing_matrix_residual_entries',0))
            if not cleanup:save()

        ipc.mkdir();dump(ipc/'control.json',{'maximum_calls':3,'maximum_soft_seconds':20,
            'rows':51574,'columns':128,'deadline_unix':time()+max(0.,SECONDS-(monotonic()-start))})
        best=None
        for iteration,beta in enumerate(BETAS,1):
            guard();pair={'iteration':iteration,'beta':str(beta),'status':'building_model','checks':[]}
            result['pairs'].append(pair);result['phase']=f'beta_{iteration:02d}_build';save()
            bf=float(beta);power=bf/D.ALPHA;mu=power*D.RADIUS**(D.ALPHA-bf)
            c0=float(c1)/mu;half=bf*D.RADIUS**(2-bf)/(2-bf)
            finite('beta_scalars',np.asarray([power,mu,c0,half]))
            if not mu>0:raise Rejected('mu must be positive; no mu>=1 requirement.')
            singular_columns={};matrices={};baseline_diagnostics={}
            for kind,g in geometry.items():
                guard();counts['transformed_singular_entries_started']+=len(keys[kind]);save()
                if counts['transformed_singular_entries_started']>154722:raise Rejected('Transformed column cap.')
                r=radii[kind]
                if kind=='K':
                    column=np.zeros(len(r));inside=r<D.RADIUS
                    column[inside]=D.RADIUS**(-bf)*np.expm1(-bf*np.log(r[inside]/D.RADIUS))
                else:
                    chi=power*np.minimum(r,D.RADIUS)**(D.ALPHA-bf)
                    finite('positive_chi',chi)
                    if np.any(chi<=0):raise Rejected('Nonpositive chi.')
                    column=chi*g.singular
                singular_columns[kind]=finite('transformed_singular',column)
                counts['transformed_singular_entries_completed']+=len(column)
                matrices[kind]=finite('all_pool_matrix',np.column_stack((column,spline_columns[kind])))
                baseline_model=finite('new_beta_baseline_model',c0*column+spline_columns[kind][:,0])
                counts['baseline_model_rows_checked']+=len(column)
                baseline_diagnostics[kind]={'original_worst':location(kind,targets[kind]-baseline_model),
                    'strengthened_worst':location(kind,new_targets[kind]-baseline_model),
                    'versus_old_baseline':location(kind,abs(baseline_model-base[kind]))}
                if beta==D.ALPHA_Q:
                    baseline_diagnostics[kind]['alpha_control_column_difference']=location(kind,abs(column-g.singular))
                np.save(ROOT/f'beta_{iteration:02d}_{kind}_singular.local.npy',column);save();guard()
            objective=finite('objective',np.asarray([half]+mass_tail))
            lower=np.asarray([c0/2,.5]+[-.5]*30+[-.25]*96);upper=np.asarray([3*c0/2,1.5]+[.5]*30+[.25]*96)
            raw=sp.csr_matrix(np.vstack((matrices['K'],matrices['A'])))
            rhs=np.r_[new_targets['K'],new_targets['A']]
            pair.update(mu=mu,chi_formula='(beta/alpha)*min(r,R)^(alpha-beta)',baseline_c=c0,
                baseline_model_diagnostics=baseline_diagnostics,singular_half_mass=half,
                rows=51574,columns=128,basis_warm_start=False)
            path=ipc/f'model_{iteration:02d}.mat';temporary=path.with_suffix('.mat.tmp')
            savemat(temporary,{'A_ub':-raw,'b_ub':-rhs[:,None],'c':objective[:,None],
                'lb':lower[:,None],'ub':upper[:,None],'base_cost':0.},appendmat=False,do_compression=False)
            temporary.replace(path);pair['input_sha256']=digest(path,guard)
            np.savez(ipc/f'model_{iteration:02d}_targets.npz',K_original=targets['K'],A_original=targets['A'],
                K_new=new_targets['K'],A_new=new_targets['A'])
            dump(ipc/f'model_{iteration:02d}_keys.json',{'rows':[[k,i] for k in ('K','A') for i in range(len(keys[k]))],
                'columns':order,'beta':str(beta),'epsilon':'1/2','A_is_lower_model_not_true_beta_value':True})
            pair['target_sha256']=digest(ipc/f'model_{iteration:02d}_targets.npz',guard)
            pair['keys_sha256']=digest(ipc/f'model_{iteration:02d}_keys.json',guard);save()
            if process is None:
                if not MATLAB.is_file():raise Rejected('Approved dedicated MATLAB unavailable; no fallback.')
                script=ROOT/'batch_solver.m';result['batch_source_sha256_at_launch']=digest(script,guard)
                batch_log=(ROOT/'matlab_batch.log').open('x')
                code="run('"+str(script).replace("'","''")+"')"
                guard();process=subprocess.Popen([str(MATLAB),'-batch',code],cwd=ROOT,stdin=subprocess.DEVNULL,
                    stdout=batch_log,stderr=subprocess.STDOUT,start_new_session=False)
                counts['batch_processes_started']+=1;child_group=os.getpgid(process.pid)
                result['batch_process']={'pid':process.pid,'pgid':child_group,'coordinator_pgid':owned_group,
                    'same_owned_group':child_group==owned_group,'shared_GUI_used':False}
                save()
                if child_group!=owned_group:raise Rejected('Owned MATLAB PGID mismatch.')
            soft=min(LP_SECONDS,max(0.,SECONDS-(monotonic()-start)))
            if soft<=0:raise StopBudget('No remaining time for next predeclared LP.')
            pair.update(status='requested',requested_seconds=monotonic()-start);counts['LP_requested']+=1;save()
            dump(ipc/f'request_{iteration:02d}.ready.json',{'iteration':iteration,'model_file':path.name,
                'input_sha256':pair['input_sha256'],'soft_seconds':soft})
            while not pair.get('result_observed'):
                guard();observe(pair)
                if pair.get('result_observed'):break
                if (ipc/'worker_error.json').exists():raise Rejected('Dedicated batch worker error; no second batch.')
                if process.poll() is not None:raise Rejected('Dedicated batch exited before result.')
                sleep(.05)
            output=pair['worker_result'];guard()
            if output.get('status')!='returned':raise Rejected('Solver execution error; no fallback.')
            if output.get('exitflag')!=1 or not output.get('finite_primal'):
                pair.update(status='rejected_solver_status',eligible=False);save();continue
            vector=finite('returned_vector',np.asarray(output['x'],dtype=float).reshape(-1))
            if vector.shape!=(128,):raise Rejected('Returned vector dimension invalid.')
            bound_error=max(0.,float(np.max(lower-vector)),float(np.max(vector-upper)))
            pair['variable_bound_error']=bound_error
            if bound_error>1e-10:
                pair.update(status='rejected_variable_box',eligible=False);save();continue
            label=f'beta_{iteration:02d}_floating'
            check,_=full_check(vector,beta,singular_columns,matrices,objective,label,pair)
            if not check['passed'] or not check['expanded_float_half_cost']<b1_cost-GAIN:
                pair.update(status='rejected_float_gate_or_gain',eligible=False);save();continue
            counts['rationalizations_started']+=1;den=2**50
            rational=[Q(round(float(x)*den),den) for x in vector]
            if rational[0]<=0 or rational[1]<0:raise Rejected('Rationalized c/s structural sign failed.')
            qvector=np.asarray([float(x) for x in rational])
            if max(float(np.max(lower-qvector)),float(np.max(qvector-upper)))>1e-10:
                pair.update(status='rejected_rational_box',eligible=False);save();continue
            check,expanded=full_check(rational,beta,singular_columns,matrices,objective,
                f'beta_{iteration:02d}_rationalized',pair)
            c,s,T,coeff=expanded
            if not check['passed'] or not check['expanded_float_half_cost']<b1_cost-GAIN:
                pair.update(status='rejected_unique_rational_gate',eligible=False);save();continue
            if any(i+j>48 for (i,j),v in T.items() if v):raise Rejected('Rationalized T support failed.')
            if any(coeff.get(ij,Q(0))!=s*a1.get(ij,Q(0))+T.get(ij,Q(0)) for ij in set(coeff)|set(a1)|set(T)):
                raise Rejected('Exact structural reconstruction failed.')
            pair.update(status='finite_candidate_eligible',eligible=True,float_half_cost=check['expanded_float_half_cost'])
            candidate={'alpha':str(D.ALPHA_Q),'beta':str(beta),'R':str(D.R_Q),'h':str(D.H_Q),'N':D.N,
                'c':str(c),'a':{f'{i},{j}':str(v) for (i,j),v in sorted(coeff.items())},
                'structure':{'s':str(s),'T':{f'{i},{j}':str(v) for (i,j),v in sorted(T.items())},'baseline_sha256':BASELINE_SHA}}
            if best is None or pair['float_half_cost']<best[0]:best=(pair['float_half_cost'],candidate,iteration)
            save()
        if best is not None:
            guard();path=ROOT/'candidate.json'
            if path.exists():raise Rejected('Only one final candidate export allowed.')
            dump(path,best[1]);result.update(candidate_exported=True,candidate_sha256=digest(path,guard),
                candidate_beta=best[1]['beta'],candidate_iteration=best[2],candidate_float_half_cost=best[0],
                next='Lead checks exact beta cost/structure and schedules beta-aware global certificate; no new theorem claim.')
        result['status']='three_fixed_beta_finite_models_closed'
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
        if 'observe' in locals():
            for pair in result['pairs']:
                if 'input_sha256' not in pair or pair.get('status')=='building_model':continue
                try:
                    observe(pair,cleanup=True)
                    if not pair.get('result_observed'):
                        pair['status']='interrupted_in_flight' if pair.get('worker_started_observed') else 'requested_not_started'
                except Exception as exc:pair['ledger_reconciliation_error']=f'{type(exc).__name__}: {exc}'
    assert counts['LP_completed']<=counts['LP_started']<=counts['LP_requested']<=3
    assert counts['batch_processes_started']<=1 and counts['residual_batches_started']<=6
    assert counts['residual_entries_started']<=309444 and counts['rationalizations_started']<=3
    assert counts['transformed_singular_entries_started']<=154722
    assert counts['primitive_formula_calls']==counts['new_T_columns']==0
    result.update(source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        batch_source_sha256_final=sha256((ROOT/'batch_solver.m').read_bytes()).hexdigest(),
        module_import_seconds_separate=_IMPORT_SECONDS,
        environment={'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__})
    save();print(json.dumps({'status':result['status'],'counts':counts,'seconds':result['elapsed_seconds'],
        'candidate_exported':result['candidate_exported']},ensure_ascii=False),flush=True)


if __name__=='__main__':main()
