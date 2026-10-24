#!/usr/bin/env python3
"""Cached finite exchange: seed + <=6 MATLAB LPs, one owned batch, 300s; no import work."""
from __future__ import annotations
from time import monotonic, time, sleep
_IMPORT_START=monotonic()
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
import math
import os
from pathlib import Path
import signal
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
PREVIOUS=ROOT.parent/'round06'
SECONDS,MAX_LP,LP_SECONDS,MAX_FINE=300.,6,20.,96
K_TOL,A_TOL,GAIN=1e-9,1e-8,5e-7
MATLAB=Path('/Applications/MATLAB_R2026a.app/bin/matlab')
DEFINITION_SHA='b5d9f5de7259b5c430ef6f0cb6e7e2c75310a41e2731497739e877d2e54cefd5'
BASELINE_SHA='21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
SEED_SHA='58babd6a834bac16eddfe0cbbd67607151904aff2dab646bdd4bf0b866bd7bdf'


class StopBudget(Exception):pass
class Rejected(Exception):pass


def dump(path,data):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    temporary.replace(path)


def digest(path,guard):
    out=sha256()
    with path.open('rb') as stream:
        while True:
            guard();block=stream.read(4*1024*1024)
            if not block:break
            out.update(block)
    return out.hexdigest()


def finite(name,value):
    if not np.all(np.isfinite(value)):raise Rejected(name+'_nonfinite')
    return value


def definitions(guard):
    path=ROOT/'dependency_round06_search.py'
    if not path.exists():path=PREVIOUS/'search.py'
    if digest(path,guard)!=DEFINITION_SHA:raise Rejected('Reviewed cache-arithmetic definition hash changed.')
    spec=importlib.util.spec_from_file_location('round07_cached_definitions',path)
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module
    spec.loader.exec_module(module)  # No main/B/L/psi evaluation or Highs construction.
    return module


def main():
    if (ROOT/'result.json').exists() or (ROOT/'ipc.local').exists():
        raise SystemExit('Existing evidence/IPC: no repeated finite round07 run.')
    start=monotonic();process=None;batch_log=None
    counts={k:0 for k in ('cache_files_verified','primitive_formula_calls',
        'residual_batches_started','residual_batches_completed','residual_entries_started',
        'residual_entries_completed','affine_matrix_rows_checked','direction_entries_started',
        'direction_entries_completed','LP_requested','LP_started','LP_completed','batch_processes_started')}
    result={'status':'started','phase':'cache_integrity','counts':counts,'solves':[],
        'candidate_exported':False,'checks':[],'cache_checks':[],
        'limits':{'seconds_internal':SECONDS,'LP_calls':MAX_LP,'LP_soft_seconds':LP_SECONDS,
        'owned_batch_processes':1,'new_primitive_formula_calls':0,'fine_classes':MAX_FINE,
        'variables':288,'residual_batches':8,'residual_scalar_checks':480000,
        'new_fine_direction_entries':5760000,'affine_existing_matrix_checks':420000},
        'LP_basis_warm_start':'Not available/claimed. A single MATLAB process runs separate linprog calls.',
        'clock_scope':'External 300s process-group cap includes imports, cache IO, MATLAB startup, IPC and output.'}
    def save():
        result['elapsed_seconds']=monotonic()-start;dump(ROOT/'result.json',result)
    def guard():
        if monotonic()-start>=SECONDS:raise StopBudget('internal_300s_deadline')
    ipc=ROOT/'ipc.local'
    save()
    try:
        # The lead-owned wrapper must create this session; no shared GUI process is controlled.
        owned_group=os.getpgrp()
        if owned_group!=os.getpid():raise Rejected('Coordinator must be leader of its own wrapper-created group.')
        result['coordinator_process']={'pid':os.getpid(),'pgid':owned_group,'sid':os.getsid(0)}
        D=definitions(guard)
        prior_path=PREVIOUS/'result.json';prior=json.loads(prior_path.read_text())
        if prior['source_sha256']!=DEFINITION_SHA or not prior['baseline_residual']['passed']:
            raise Rejected('Provenance of prior cached baseline invalid.')
        seed_path=PREVIOUS/'matlab_result.json';seed=json.loads(seed_path.read_text())
        if digest(seed_path,guard)!=SEED_SHA:raise Rejected('The approved saved MATLAB seed hash changed.')
        if seed['exitflag']!=1 or seed['columns']!=192 or seed['rows']!=30968:
            raise Rejected('Only the fixed saved MATLAB optimal vector is the approved seed.')
        seed_t=finite('saved_seed',np.asarray(seed['x'],dtype=float))
        if seed_t.shape!=(192,):raise Rejected('Seed dimension changed.')
        source=ROOT.parent/'round01/scaled_certificate.json'
        if digest(source,guard)!=BASELINE_SHA:raise Rejected('Proven K1 input changed.')
        spec=json.loads(source.read_text());c1=Q(spec['c'])
        a1={tuple(map(int,k.split(','))):Q(v) for k,v in spec['a'].items()}
        if (Q(spec['alpha']),Q(spec['R']),Q(spec['h']),int(spec['N']))!=(D.ALPHA_Q,D.R_Q,D.H_Q,D.N):
            raise Rejected('Fixed alpha/R/h/N changed.')
        point_path=PREVIOUS/'point_manifest.json'
        if digest(point_path,guard)!=prior['dependencies']['previous_point_manifest_sha256']:
            raise Rejected('Cached point order/source hash mismatch.')
        point_document=json.loads(point_path.read_text())
        keys={kind:[(Q(u),Q(v)) for u,v in point_document[kind]] for kind in ('K','A')}
        points={kind:np.asarray([(float(u),float(v)) for u,v in keys[kind]]) for kind in keys}
        if (len(keys['K']),len(keys['A']))!=(26604,24970):raise Rejected('Fixed pool dimensions changed.')
        targets={'K':np.asarray([1. if u+v<=1 else 0. for u,v in keys['K']]),'A':np.zeros(len(keys['A']))}
        for entry in prior['local_cache_files']:
            path=PREVIOUS/entry['file']
            if path.stat().st_size!=entry['bytes'] or digest(path,guard)!=entry['sha256']:
                raise Rejected('Prior local cache missing/corrupt: '+path.name)
            counts['cache_files_verified']+=1
            result['cache_checks'].append({'file':entry['file'],'sha256':entry['sha256'],'passed':True});save()
        cache=PREVIOUS/'cache.local'
        geometry={};base={};matrices={}
        for kind in ('K','A'):
            g=object.__new__(D.Geometry)
            g.kind,g.points,g.guard,g.counts,g.save=kind,points[kind],guard,counts,save
            g.B={'fine':[sp.load_npz(cache/f'{kind}_fine_B{axis}.npz') for axis in (0,1)]}
            g.L={'fine':[np.load(cache/f'{kind}_fine_L{axis}.npy',mmap_mode='r') for axis in (0,1)]} if kind=='A' else {}
            g.singular=np.load(cache/f'{kind}_singular.npy',mmap_mode='r')
            base[kind]=np.load(cache/f'{kind}_baseline.npy',mmap_mode='r')
            matrices[kind]=np.load(cache/f'{kind}_directions_192.npy',mmap_mode='r')
            if matrices[kind].shape!=(len(points[kind]),192) or base[kind].shape!=(len(points[kind]),):
                raise Rejected('Cached column/baseline dimensions incorrect.')
            if g.singular.shape!=(len(points[kind]),) or any(b.shape!=(len(points[kind]),161) for b in g.B['fine']):
                raise Rejected('Cached singular/B dimensions incorrect.')
            if kind=='A' and any(l.shape!=(len(points[kind]),161) for l in g.L['fine']):
                raise Rejected('Cached L dimensions incorrect.')
            g.columns=[matrices[kind][:,j] for j in range(192)]
            g.column_keys={('global_c',),('global_S',)}|{('coarse',*ij) for ij in D.COARSE}
            geometry[kind]=g
        expected_order=['relative_c','relative_S']+[f'coarse_{i}_{j}' for i,j in D.COARSE]
        if prior['direction_order']!=expected_order:raise Rejected('Cached direction order changed.')
        rows_path=PREVIOUS/'interface_01_row_keys.json'
        if digest(rows_path,guard)!=prior['load_snapshots'][0]['input_files']['_row_keys.json']:
            raise Rejected('Initial intended row-key source changed.')
        active=[(str(k),int(i)) for k,i in json.loads(rows_path.read_text())['rows']]
        if len(active)!=30968 or len(set(active))!=len(active):raise Rejected('Initial active rows invalid.')
        if any(k not in keys or not 0<=i<len(keys[k]) for k,i in active):raise Rejected('Invalid initial row index.')
        active_set=set(active)
        refinement_path=ROOT.parent/'round05/refinement_identity.json'
        if digest(refinement_path,guard)!=prior['dependencies']['refinement_identity_sha256']:
            raise Rejected('Independent refinement evidence changed.')
        mask=D.refinement(refinement_path)
        direction_maps=[D.refined_direction(ij,mask) for ij in D.COARSE]
        fine=[]
        spline_mass=sum((D.H_Q**2*len(D.orbit(*ij))*v/2 for ij,v in a1.items()),Q(0))
        singular_half=D.ALPHA*D.RADIUS**(2-D.ALPHA)/(2-D.ALPHA)
        base_cost=float(c1)*singular_half+float(spline_mass)
        objective=[float(c1)*singular_half,float(spline_mass)]+[float(D.COARSE_Q**2*len(D.orbit(*ij))/2) for ij in D.COARSE]
        lower=[-.2,-.5]+[-.5]*190;upper=[.2,.5]+[.5]*190
        result['dependencies']={'round06_definition_sha256':DEFINITION_SHA,
            'prior_result_sha256':digest(prior_path,guard),'seed_sha256':digest(seed_path,guard),
            'baseline_sha256':BASELINE_SHA,'point_manifest_sha256':digest(point_path,guard),
            'refinement_sha256':digest(refinement_path,guard)}
        result['baseline_check_reused']='Prior passed baseline with verified cache/source/point hashes; no new baseline batch.'
        result['pool_counts']={k:len(v) for k,v in keys.items()}
        ipc.mkdir();dump(ipc/'control.json',{'maximum_calls':MAX_LP,'maximum_soft_seconds':LP_SECONDS,
            'deadline_unix':time()+max(0.,SECONDS-(monotonic()-start))})
        save();guard()

        def expand(t):
            coeff={ij:(1.+t[1])*float(v) for ij,v in a1.items()}
            for value,direction in zip(t[2:],direction_maps):
                for ij,weight in direction.items():coeff[ij]=coeff.get(ij,0.)+float(value)*float(weight)
            return (1.+t[0])*float(c1),coeff

        def cost(c,coeff):
            return float(c)*singular_half+sum(float(D.H_Q**2*len(D.orbit(*ij))/2)*float(v) for ij,v in coeff.items())

        def check_full(c,coeff,label,t=None):
            guard();counts['residual_batches_started']+=1
            if counts['residual_batches_started']>8:raise Rejected('Eight residual-batch cap.')
            matrix=D.full_matrix(coeff);values={};violations={};record={'label':label,'status':'started'}
            result['checks'].append(record);save()
            for kind,g in geometry.items():
                guard();counts['residual_entries_started']+=len(keys[kind]);save()
                values[kind]=finite('full_'+kind,float(c)*g.singular+g.spline(matrix))
                violations[kind]=targets[kind]-values[kind]
                bad=np.flatnonzero(violations[kind]>(K_TOL if kind=='K' else A_TOL))
                record[kind]={'minimum_value':float(values[kind].min()),'maximum_violation':float(violations[kind].max()),
                    'bad_indices':bad.tolist(),'bad_count':len(bad)}
                counts['residual_entries_completed']+=len(keys[kind])
                if t is not None:
                    affine=finite('affine_'+kind,base[kind]+matrices[kind]@t)
                    record[kind]['affine_full_max_difference']=float(np.max(abs(affine-values[kind])))
                    indices=np.asarray([i for k,i in active if k==kind],dtype=int)
                    tolerance=K_TOL if kind=='K' else A_TOL
                    record[kind]['active_affine_bad_count']=int(np.count_nonzero(targets[kind][indices]-affine[indices]>tolerance))
                    record[kind]['active_full_bad_count']=int(np.count_nonzero(violations[kind][indices]>tolerance))
                    counts['affine_matrix_rows_checked']+=len(keys[kind])
                guard();save()
            record.update(status='completed',passed=all(record[k]['bad_count']==0 for k in ('K','A')),
                float_half_cost=cost(c,coeff))
            if t is not None:
                record['affine_half_cost']=base_cost+float(np.asarray(objective)@t)
                record['cost_difference']=abs(record['float_half_cost']-record['affine_half_cost'])
            np.savez(ROOT/(label+'_full_values.npz'),K=values['K'],A=values['A'])
            dump(ROOT/(label+'_coefficients.json'),{'c':str(c),'a':{f'{i},{j}':str(v) for (i,j),v in coeff.items()}})
            counts['residual_batches_completed']+=1;save();return values,violations,record

        def separate(violations,record):
            worst=[]
            for kind,cap in (('K',256),('A',512)):
                omitted=[i for i in record[kind]['bad_indices'] if (kind,i) not in active_set]
                omitted.sort(key=lambda i:(-float(violations[kind][i]),*keys[kind][i]))
                worst.extend((kind,i) for i in omitted[:cap])
            if not worst:return False
            witnesses=sorted(worst,key=lambda row:(-float(violations[row[0]][row[1]])/(K_TOL if row[0]=='K' else A_TOL),
                *keys[row[0]][row[1]],row[0]))
            proposed=[];existing=set(fine);limit=min(24,MAX_FINE-len(fine))
            for kind,index in (witnesses if limit else []):
                u,v=keys[kind][index];iu,iv=int(u//D.H_Q),int(v//D.H_Q)
                for i in range(iu-1,iu+3):
                    for j in range(iv-1,iv+3):
                        ij=D.canonical(i,j)
                        if ij[0]<=D.N and ij not in existing:
                            proposed.append(ij);existing.add(ij)
                            if len(proposed)>=limit:break
                    if len(proposed)>=limit:break
                if len(proposed)>=limit:break
            for ij in proposed:
                for g in geometry.values():
                    if counts['direction_entries_started']+len(g.points)>5760000:raise Rejected('Fine direction-entry cap.')
                    g.orbit_column('fine',ij)
                fine.append(ij);direction_maps.append({ij:Q(1)})
                objective.append(float(D.H_Q**2*len(D.orbit(*ij))/2));lower.append(-.25);upper.append(.25)
            if proposed:
                for kind,g in geometry.items():
                    matrices[kind]=g.matrix();path=ROOT/f'{kind}_directions_{len(objective)}.local.npy'
                    np.save(path,matrices[kind])
            active.extend(worst);active_set.update(worst)
            record['separation']={'new_K_rows':sum(k=='K' for k,i in worst),'new_A_rows':sum(k=='A' for k,i in worst),
                'new_fine_classes':[list(ij) for ij in proposed]};save();guard();return True

        def write_model(iteration):
            guard();raw=sp.csr_matrix(np.vstack([matrices[k][i] for k,i in active]))
            rhs=np.asarray([targets[k][i]-base[k][i] for k,i in active])
            finite('raw_model',raw.data);finite('raw_rhs',rhs)
            path=ipc/f'model_{iteration:02d}.mat'
            if path.exists():raise Rejected('Repeated immutable model number.')
            temporary=path.with_suffix('.mat.tmp')
            savemat(temporary,{'A_ub':-raw,'b_ub':-rhs[:,None],'c':np.asarray(objective)[:,None],
                'lb':np.asarray(lower)[:,None],'ub':np.asarray(upper)[:,None],'base_cost':base_cost},
                appendmat=False,do_compression=False)
            temporary.replace(path)
            source_hash=digest(path,guard)
            sp.save_npz(ipc/f'model_{iteration:02d}_raw.npz',raw,compressed=False)
            dump(ipc/f'model_{iteration:02d}_keys.json',{'rows':[[k,i] for k,i in active],
                'columns':['relative_c','relative_S']+[f'coarse_{i}_{j}' for i,j in D.COARSE]+[f'fine_{i}_{j}' for i,j in fine],
                'raw_model_sha256':source_hash,'solver_loaded_model_not_exposed_by_MATLAB':True})
            guard();return path,source_hash

        def start_batch():
            nonlocal process,batch_log
            if counts['batch_processes_started'] or process is not None:raise Rejected('Second MATLAB batch forbidden.')
            if not MATLAB.is_file():raise Rejected('Approved MATLAB batch binary unavailable; no fallback.')
            script=ROOT/'batch_solver.m'
            result['batch_source_sha256_at_launch']=digest(script,guard)
            code="run('"+str(script).replace("'","''")+"')"
            batch_log=(ROOT/'matlab_batch.log').open('x')
            guard()
            process=subprocess.Popen([str(MATLAB),'-batch',code],cwd=ROOT,stdin=subprocess.DEVNULL,
                stdout=batch_log,stderr=subprocess.STDOUT,start_new_session=False)
            counts['batch_processes_started']+=1
            child_group=os.getpgid(process.pid)
            result['batch_process']={'pid':process.pid,'pgid':child_group,'coordinator_pgid':owned_group,
                'same_owned_group':child_group==owned_group,'shared_GUI_used':False,'started_seconds':monotonic()-start}
            save()
            if child_group!=owned_group:raise Rejected('MATLAB launcher did not inherit the owned process group.')

        def update_started(record):
            path=ipc/f"started_{record['iteration']:02d}.json"
            if path.exists() and not record.get('worker_started_observed'):
                state=json.loads(path.read_text())
                if state['iteration']!=record['iteration'] or state['input_sha256']!=record['input_sha256']:
                    raise Rejected('Started ledger not bound to this immutable model.')
                record['worker_started_observed']=True;record['worker_started']=state
                counts['LP_started']+=1;save()

        def request_solve(iteration):
            path,source_hash=write_model(iteration)
            if process is None:start_batch()
            soft=min(LP_SECONDS,max(0.,SECONDS-(monotonic()-start)))
            if soft<=0:raise StopBudget('No remaining time for queued LP.')
            record={'iteration':iteration,'status':'queued','input_sha256':source_hash,
                'variables':len(objective),'rows':len(active),'soft_seconds':soft,
                'requested_seconds':monotonic()-start,'basis_warm_start':False}
            result['solves'].append(record);counts['LP_requested']+=1;save()
            ready=ipc/f'request_{iteration:02d}.ready.json'
            if ready.exists():raise Rejected('Repeated request ready marker.')
            dump(ready,{'iteration':iteration,'model_file':path.name,'input_sha256':source_hash,'soft_seconds':soft})
            response=ipc/f'result_{iteration:02d}.ready.json'
            while not response.exists():
                guard();update_started(record)
                if (ipc/'worker_error.json').exists():
                    record['worker_error']=json.loads((ipc/'worker_error.json').read_text())
                    raise Rejected('Owned batch worker failed; no second batch.')
                if process.poll() is not None:raise Rejected('Owned batch exited before its result marker.')
                sleep(min(.05,max(.001,SECONDS-(monotonic()-start))))
            update_started(record)
            marker=json.loads(response.read_text())
            expected=f'result_{iteration:02d}.json'
            if (marker['iteration'],marker['input_sha256'],marker['result_file'])!=(iteration,source_hash,expected):
                raise Rejected('Ready result is stale or bound to a different model.')
            result_path=ipc/expected
            if digest(result_path,guard)!=marker['result_sha256']:raise Rejected('Result SHA mismatch.')
            output=json.loads(result_path.read_text())
            if output['iteration']!=iteration or output['input_sha256']!=source_hash:
                raise Rejected('Worker result model binding changed.')
            if not record.get('worker_started_observed'):raise Rejected('Completed result without started ledger.')
            record.update(status='completed',worker_result=output,result_sha256=marker['result_sha256'],
                observed_completed_seconds=monotonic()-start)
            counts['LP_completed']+=1;save();guard()
            if output.get('status')!='returned' or output.get('exitflag')!=1 or not output.get('finite_primal'):
                raise Rejected('MATLAB not optimal/finite; retained output but no fallback or candidate.')
            t=finite('MATLAB_vector',np.asarray(output['x'],dtype=float).reshape(-1))
            if t.shape!=(len(objective),):raise Rejected('Worker vector dimension mismatch.')
            return t,record

        def variable_gate(t):
            violation=max(0.,float(np.max(np.asarray(lower)-t)),float(np.max(t-np.asarray(upper))))
            if violation>1e-10:raise Rejected('Vector outside the fixed direction box.')
            return violation

        result['phase']='seed_full_pool';result['seed_variable_bound_error']=variable_gate(seed_t)
        c,coeff=expand(seed_t)
        values,violations,check=check_full(c,coeff,'seed',seed_t)
        result['seed_check']=check;best=None
        if check['passed']:
            best=(c,coeff,check['float_half_cost']);result['termination']='seed_full_pool_feasible'
        else:
            for iteration in range(1,MAX_LP+1):
                if not separate(violations,check):
                    result['termination']='no_new_violated_rows';break
                result['phase']='MATLAB_LP_'+str(iteration);save()
                t,record=request_solve(iteration);record['variable_bound_error']=variable_gate(t)
                c,coeff=expand(t)
                values,violations,check=check_full(c,coeff,f'LP_{iteration:02d}',t)
                record['complete_check']=check
                if any(check[k]['active_affine_bad_count'] or check[k]['active_full_bad_count'] for k in ('K','A')):
                    raise Rejected('Optimal vector violates intended active raw/full rows; no rhs repair.')
                if check['passed']:
                    best=(c,coeff,check['float_half_cost']);result['termination']='LP_full_pool_feasible';break
                if iteration==MAX_LP:result['termination']='six_LP_cap_not_converged'
                save();guard()
        result['fine_classes']=[list(ij) for ij in fine]
        result['direction_order']=['relative_c','relative_S']+[f'coarse_{i}_{j}' for i,j in D.COARSE]+[f'fine_{i}_{j}' for i,j in fine]
        if best is not None and best[2]<base_cost-GAIN:
            c,coeff,value=best;den=2**50
            cq=Q(round(c*den),den);aq={ij:Q(round(v*den),den) for ij,v in coeff.items()}
            if cq<=0:raise Rejected('Rationalized singular coefficient nonpositive.')
            _,_,rational_check=check_full(cq,aq,'rationalized')
            result['rationalized_check']=rational_check;save();guard()
            if not rational_check['passed']:raise Rejected('One rationalized candidate failed; no repair.')
            nz={ij:v for ij,v in aq.items() if v}
            actual_U=(max(i for i,j in nz)+2)*D.H_Q
            if not D.R_Q<=actual_U<=Q(41,16):raise Rejected('Actual exported support scan window invalid.')
            target=ROOT/'candidate.json'
            if target.exists():raise Rejected('Second rationalization forbidden.')
            dump(target,{'alpha':str(D.ALPHA_Q),'R':str(D.R_Q),'h':str(D.H_Q),'N':D.N,'c':str(cq),
                'a':{f'{i},{j}':str(v) for (i,j),v in sorted(nz.items())}})
            result.update(candidate_exported=True,candidate_sha256=digest(target,guard),
                candidate_float_half_cost=cost(cq,aq),support_U=str(actual_U),rationalization_denominator=str(den),
                next='Lead schedules exact cost and the two previously planned full-domain checks; no automatic theorem claim.')
        result['status']='finite_cached_MATLAB_exchange_closed'
    except StopBudget as exc:result.update(status='interrupted_saved_checkpoint',reason=str(exc))
    except Exception as exc:result.update(status='failed_or_rejected',reason=f'{type(exc).__name__}: {exc}')
    finally:
        if ipc.exists() and not (ipc/'stop.json').exists():
            dump(ipc/'stop.json',{'stop':True,'reason':'coordinator_finally','seconds':monotonic()-start})
        if process is not None:
            cleanup={'owned_pid':process.pid,'started_seconds':monotonic()-start,'stop_requested':True}
            try:
                process.wait(timeout=min(5.,max(.05,SECONDS-(monotonic()-start))))
            except subprocess.TimeoutExpired:
                cleanup['terminate_sent']=True;process.terminate()
                try:process.wait(timeout=min(2.,max(.05,SECONDS-(monotonic()-start))))
                except subprocess.TimeoutExpired:cleanup['kill_sent']=True;process.kill();process.wait(timeout=1.)
            cleanup.update(returncode=process.returncode,completed_seconds=monotonic()-start,
                leftover_group_cleanup='Lead-owned outer wrapper cleans only the group it created, never the shared GUI.')
            result['batch_cleanup']=cleanup
        if batch_log is not None:batch_log.close()
        # Reconcile tiny immutable IPC ledgers after cleanup, even when the deadline
        # interrupted polling. These completions are evidence only, never reused as candidates.
        for record in result['solves']:
            try:
                number=record['iteration'];started_path=ipc/f'started_{number:02d}.json'
                if started_path.exists() and not record.get('worker_started_observed'):
                    state=json.loads(started_path.read_text())
                    if state['iteration']!=number or state['input_sha256']!=record['input_sha256']:
                        raise Rejected('Cleanup started ledger binding mismatch.')
                    record['worker_started_observed']=True;record['worker_started']=state;counts['LP_started']+=1
                ready=ipc/f'result_{number:02d}.ready.json'
                if ready.exists() and record['status']!='completed':
                    marker=json.loads(ready.read_text());path=ipc/f'result_{number:02d}.json'
                    if (marker['iteration'],marker['input_sha256'],marker['result_file'])!=(number,record['input_sha256'],path.name):
                        raise Rejected('Cleanup result marker binding mismatch.')
                    if sha256(path.read_bytes()).hexdigest()!=marker['result_sha256']:
                        raise Rejected('Cleanup result hash mismatch.')
                    if not record.get('worker_started_observed'):raise Rejected('Cleanup completion without start.')
                    record.update(status='completed_during_cleanup_not_used',worker_result=json.loads(path.read_text()),
                                  result_sha256=marker['result_sha256']);counts['LP_completed']+=1
                elif record['status']!='completed':
                    record['status']='interrupted_in_flight' if record.get('worker_started_observed') else 'requested_not_started'
            except Exception as exc:
                record['ledger_reconciliation_error']=f'{type(exc).__name__}: {exc}'
    assert counts['batch_processes_started']<=1 and counts['LP_requested']<=6
    assert counts['LP_completed']<=counts['LP_started']<=counts['LP_requested']
    assert counts['residual_batches_started']<=8 and counts['residual_entries_started']<=480000
    assert counts['direction_entries_started']<=5760000 and counts['affine_matrix_rows_checked']<=420000
    result.update(source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        batch_source_sha256_final=sha256((ROOT/'batch_solver.m').read_bytes()).hexdigest(),
        module_import_seconds_separate=_IMPORT_SECONDS,
        environment={'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__},
        scope='Finite numerical model only; no new B/L/psi formulas, full-domain checker, exact cost proof or automatic restart.')
    save()
    print(json.dumps({'status':result['status'],'counts':counts,'seconds':result['elapsed_seconds'],
                      'candidate_exported':result['candidate_exported']},ensure_ascii=False),flush=True)


if __name__=='__main__':main()
