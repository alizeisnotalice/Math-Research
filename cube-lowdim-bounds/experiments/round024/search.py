#!/usr/bin/env python3
"""Round05 finite coarse/fine exchange: <=6 LP, 300s; import runs no experiment."""
from __future__ import annotations
from time import monotonic
_IMPORT_START = monotonic()
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import itertools
import json
import math
import os
from pathlib import Path
import sys
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
import numpy as np
import scipy
import scipy.sparse as sp
import highspy
_IMPORT_SECONDS = monotonic()-_IMPORT_START
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
ALPHA_Q, R_Q, H_Q, COARSE_Q = Q(47,40), Q(33,20), Q(1,32), Q(1,8)
ALPHA, RADIUS, N = float(ALPHA_Q), float(R_Q), 80
SECONDS, LP_SECONDS, MAX_LP, MAX_FINE = 300., 20., 6, 96
K_TOL, A_TOL, GAIN = 1e-9, 1e-8, 5e-7
COARSE = [(i,j) for i in range(19) for j in range(i+1)]


class BudgetStop(Exception): pass
class Rejected(Exception): pass


def dump(path, data):
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False)+'\n')
    temporary.replace(path)


def finite(name, data):
    if not np.all(np.isfinite(data)): raise Rejected(name+'_nonfinite')
    return data


def orbit(i,j):
    return sorted({(x,y) for x in (i,-i) for y in (j,-j)} |
                  {(y,x) for x in (i,-i) for y in (j,-j)})


def canonical(i,j): return max(abs(i),abs(j)), min(abs(i),abs(j))


def load_helper():
    for path in (ROOT/'dependency_round02_search.py', ROOT.parent/'round02/search.py'):
        if path.is_file():
            spec = importlib.util.spec_from_file_location('round05_reviewed_round02', path)
            module = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = module; spec.loader.exec_module(module)
            return module, path
    raise FileNotFoundError('Reviewed round02 definitions not found.')


def legacy_sets():
    """Exact coordinates of round03 unions; no formula/LP evaluations."""
    U=Q(9,4)
    def grid(den):
        return {((Q(m)+Q(1,3))/den,(Q(n)+Q(1,3))/den)
                for m in range(int(U*den)+1) if (Q(m)+Q(1,3))/den<=U
                for n in range(m+1)}
    k0=grid(32)|{(Q(m,32),Q(n,32)) for m in range(73) for n in range(m+1) if m or n}
    k0|={(1-Q(j,64),Q(j,64)) for j in range(33)}
    k1=grid(64)|{(1-Q(j,128),Q(j,128)) for j in range(65)}
    origin={(Q(1,2**k)*(1-f),Q(1,2**k)*f) for k in range(3,19)
            for f in (Q(1,2),Q(1,4),Q(1,16),Q(1,256))}
    valley={(Q('1.541724')+Q(i,4096),Q('.132603')+Q(j,4096))
            for i in range(-2,3) for j in range(-2,3)}
    a0,a1=grid(32)|origin|valley,grid(64)|origin|valley
    for dst,ks,den,js in ((a0,range(6,15),64,range(1,33)),
                          (a1,range(4,16),128,range(1,65))):
        for sign,k,j in itertools.product((-1,1),ks,js):
            r=R_Q+sign*Q(1,2**k); v=r*Q(j,den); dst.add((r-v,v))
    aa={(u,v) for u,v in a0|a1 if 0<v<=u<=U and u+v!=R_Q}
    if len(k0|k1)!=15816 or len(aa)!=14693: raise Rejected('Legacy exact pool count changed.')
    return k0|k1,aa


def pools(points_path):
    ki,ai=legacy_sets()
    manifest=json.loads(points_path.read_text())
    if len(manifest['points'])!=12: raise Rejected('Frozen diagnostic point schedule changed.')
    for p in manifest['points']:
        (ki if p['quantity']=='K' else ai).add((Q(p['u']),Q(p['v'])))
    ki|={(Q(j,1024),1-Q(j,1024)) for j in range(512,1025)}
    kv,av=set(ki),set(ai)
    radial={(r*(1+q)/2,r*(1-q)/2) for r in (Q(13*k,640) for k in range(1,129))
            for q in (Q(2*j+1,128) for j in range(64))}
    square={(Q(41*i,1024),Q(41*j,1024)) for i in range(65) for j in range(i+1)}
    kv|=radial|square; av|=radial|square
    kv={(u,v) for u,v in kv if 0<u+v and 0<=v<=u}
    av={(u,v) for u,v in av if 0<v<=u and u+v!=R_Q}
    if len(kv)>30000 or len(av)>30000: raise Rejected('Fixed validation pool cap exceeded.')
    keys={'K':sorted(kv),'A':sorted(av)}
    initial={'K':ki,'A':ai}
    arrays={kind:np.asarray([(float(u),float(v)) for u,v in ps]) for kind,ps in keys.items()}
    selected={kind:[i for i,p in enumerate(ps) if p in initial[kind]] for kind,ps in keys.items()}
    targets={'K':np.asarray([1. if u+v<=1 else 0. for u,v in keys['K']]),
             'A':np.zeros(len(keys['A']))}
    return keys,arrays,selected,targets


def full_matrix(coefficients):
    matrix=np.zeros((2*N+1,2*N+1))
    for (i,j),value in coefficients.items():
        for x,y in orbit(i,j): matrix[x+N,y+N]=float(value)
    return matrix


def refinement(path):
    proof=json.loads(path.read_text())
    if not proof.get('passed'): raise Rejected('Independent refinement identity unavailable.')
    mask=[Q(x) for x in proof['mask']]
    if len(mask)!=13 or proof['indices']!=list(range(-6,7)) or sum(mask)!=4:
        raise Rejected('Refinement mask manifest inconsistent.')
    # Only read the independent symbolic evidence; never execute its generator.
    return list(zip(range(-6,7),mask))


def refined_direction(ij,mask):
    full={}
    for i,j in orbit(*ij):
        for k,x in mask:
            for ell,y in mask:
                key=(4*i+k,4*j+ell)
                full[key]=full.get(key,Q(0))+x*y
    out={}
    for (i,j),value in full.items():
        c=canonical(i,j)
        if max(c)>N: raise Rejected('Refinement exceeded declared N window.')
        if c in out and out[c]!=value: raise Rejected('Refinement lost D4 symmetry.')
        out[c]=value
    return out


class Geometry:
    """B/L primitives cached once; every direction column is built at most once."""
    def __init__(self,kind,points,bspl,singular,guard,counts,save):
        self.kind,self.points,self.guard,self.counts,self.save=kind,points,guard,counts,save
        self.B={}; self.L={}; self.column_keys=set(); self.columns=[]
        for label,h,n in (('fine',float(H_Q),N),('coarse',float(COARSE_Q),18)):
            bb=[[],[]]; ll=[[],[]]
            indices=np.arange(-n,n+1)
            for first in range(0,len(points),256):
                guard(); p=points[first:first+256]
                for axis in (0,1):
                    arg=p[:,axis,None]/h-indices[None,:]
                    counts['primitive_B_entries_started']+=arg.size
                    b=finite('B',bspl.Bv(arg)); bb[axis].append(sp.csr_matrix(b))
                    counts['primitive_B_entries_completed']+=arg.size
                    if kind=='A':
                        counts['primitive_L_entries_started']+=arg.size
                        l=finite('true_LB',bspl.LBv(arg,ALPHA)*h**(-ALPHA));ll[axis].append(l)
                        counts['primitive_L_entries_completed']+=arg.size
                guard();save()
            self.B[label]=[sp.vstack(b,format='csr') for b in bb]
            if kind=='A': self.L[label]=[np.vstack(l) for l in ll]
        values=[]
        for first in range(0,len(points),256):
            guard();p=points[first:first+256]; counts['singular_entries_started']+=len(p)
            if kind=='K': value=np.maximum(p.sum(axis=1)**(-ALPHA)-RADIUS**(-ALPHA),0.)
            else: value=singular.A_psi_vec(ALPHA,RADIUS,p[:,0],p[:,1])
            values.append(finite('unmodified_singular',value))
            counts['singular_entries_completed']+=len(p);guard();save()
        self.singular=np.concatenate(values)

    def spline(self,matrix):
        self.guard();bu,bv=self.B['fine']
        if self.kind=='K':
            return finite('full_S',np.asarray(bv.multiply(bu@matrix).sum(axis=1)).ravel())
        lu,lv=self.L['fine']
        return finite('full_AS',np.sum(np.asarray(bv@matrix.T)*lu,axis=1)+
                      np.sum(np.asarray(bu@matrix)*lv,axis=1))

    def add(self,key,builder):
        if key in self.column_keys: raise Rejected('Direction formula column rebuilt.')
        self.guard();self.counts['direction_entries_started']+=len(self.points)
        if self.counts['direction_entries_started']>17280000:raise Rejected('Direction entry cap at reservation.')
        self.counts['direction_column_in_progress']={'pool':self.kind,'key':list(key),'entries':len(self.points)}
        self.save()
        values=builder();self.guard()
        self.columns.append(finite('direction_column',np.asarray(values).ravel()))
        self.column_keys.add(key);self.counts['direction_entries_completed']+=len(self.points)
        self.counts['direction_column_in_progress']=None;self.save()

    def orbit_column(self,label,ij):
        key=(label,*ij)
        if key in self.column_keys: raise Rejected('Repeated orbit formula build.')
        def build():
            bu,bv=self.B[label];shift=N if label=='fine' else 18
            out=np.zeros(len(self.points))
            for i,j in orbit(*ij):
                self.guard();bui=bu[:,i+shift].toarray().ravel();bvj=bv[:,j+shift].toarray().ravel()
                if self.kind=='K':out+=bui*bvj
                else:
                    lu,lv=self.L[label];out+=lu[:,i+shift]*bvj+bui*lv[:,j+shift]
            return out
        self.add(key,build)

    def matrix(self): return np.column_stack(self.columns)


def main():
    if (ROOT/'result.json').exists():
        raise SystemExit('Existing evidence: no repeat of the finite round05 discovery.')
    start=monotonic()
    counts={name:0 for name in ('LP_started','LP_completed','primitive_B_entries_started',
        'primitive_B_entries_completed','primitive_L_entries_started','primitive_L_entries_completed',
        'singular_entries_started','singular_entries_completed','direction_entries_started',
        'direction_entries_completed','baseline_entries_started','baseline_entries_completed',
        'residual_entries_started','residual_entries_completed','residual_batches_started',
        'residual_batches_completed','affine_matrix_rows_checked')}
    result={'status':'started','phase':'prepare','counts':counts,'solves':[],
            'candidate_exported':False,'limits':{'seconds_internal':SECONDS,'LP_calls':MAX_LP,
            'LP_soft_seconds':LP_SECONDS,'fine_classes':MAX_FINE,'variables':288,
            'per_pool_points':30000,'residual_batches':8,'residual_scalar_checks':480000,
            'direction_formula_entries':17280000,'MATLAB_existing_matrix_checks_separate':60000},
            'clock_scope':'Import separately recorded; outer wrapper must hard-kill entire process group at 300s.'}
    def save():
        result['elapsed_seconds']=monotonic()-start;dump(ROOT/'result.json',result)
    def guard():
        if monotonic()-start>=SECONDS: raise BudgetStop('internal_300s_deadline')
    def check_status(status,operation):
        if status!=highspy.HighsStatus.kOk: raise Rejected(operation+'_not_kOk')
    save()
    try:
        helper,helper_path=load_helper()
        source=ROOT/'input_baseline.json'
        if not source.exists():source=ROOT.parent/'round01/scaled_certificate.json'
        basis=helper.Basis(source)
        if Q(basis.spec['alpha'])!=ALPHA_Q or Q(basis.spec['R'])!=R_Q:
            raise Rejected('Fixed proven alpha/R changed.')
        c1=Q(basis.spec['c']); a1={tuple(map(int,k.split(','))):Q(v) for k,v in basis.spec['a'].items()}
        points_source=ROOT/'input_diagnostic_points.json'
        if not points_source.exists():points_source=ROOT.parent/'round04/points.json'
        refinement_path=ROOT/'refinement_identity.json'; mask=refinement(refinement_path)
        coarse_maps=[refined_direction(ij,mask) for ij in COARSE]
        # Check half-integrals under the actual signed-lattice expansion.
        for ij,direction in zip(COARSE,coarse_maps):
            expanded=sum((H_Q**2*len(orbit(*key))*v/2 for key,v in direction.items()),Q(0))
            if expanded!=COARSE_Q**2*len(orbit(*ij))/2:raise Rejected('Coarse D4 half mass mismatch.')
        keys,points,selected,targets=pools(points_source)
        result['pool_counts']={kind:{'validation':len(keys[kind]),'initial':len(selected[kind])} for kind in keys}
        dump(ROOT/'point_manifest.json',{kind:[[str(u),str(v)] for u,v in ps] for kind,ps in keys.items()})
        bspl,singular=helper.load_formulas();fd=helper.formula_directory()
        result['dependencies']={'reviewed_round02_sha256':sha256(helper_path.read_bytes()).hexdigest(),
            'baseline_sha256':sha256(source.read_bytes()).hexdigest(),
            'points_source_sha256':sha256(points_source.read_bytes()).hexdigest(),
            'refinement_identity_sha256':sha256(refinement_path.read_bytes()).hexdigest(),
            'formula_sha256':{n:sha256((fd/n).read_bytes()).hexdigest() for n in ('bspl_closed.py','sing2.py','sing3.py')}}
        result['support']={'old_actual_bound':'9/4','coarse_actual_bound':'5/2',
                           'declared_window_bound':'41/16','fine_activation_max_I':80,
                           'old_actual_support_may_expand':True}
        guard();save()
        geometry={}
        for kind in ('K','A'):
            result['phase']='cache_'+kind;save()
            geometry[kind]=Geometry(kind,points[kind],bspl,singular,guard,counts,save)
        baseline_matrix=full_matrix(a1)
        base={};spline={}
        for kind,g in geometry.items():
            counts['baseline_entries_started']+=len(points[kind]);guard()
            spline[kind]=g.spline(baseline_matrix)
            base[kind]=finite('baseline_'+kind,float(c1)*g.singular+spline[kind])
            counts['baseline_entries_completed']+=len(points[kind]);save()

        def residual(c,coefficients,label,known=None):
            counts['residual_batches_started']+=1
            if counts['residual_batches_started']>8:raise Rejected('Residual batch cap exceeded.')
            matrix=None if known is not None else full_matrix(coefficients)
            values={};violations={};summary={'label':label}
            for kind,g in geometry.items():
                guard();counts['residual_entries_started']+=len(points[kind]);save()
                values[kind]=known[kind] if known is not None else finite('expanded_'+kind,float(c)*g.singular+g.spline(matrix))
                violations[kind]=targets[kind]-values[kind]
                tolerance=K_TOL if kind=='K' else A_TOL
                bad=np.flatnonzero(violations[kind]>tolerance)
                summary[kind]={'minimum_value':float(values[kind].min()),
                    'maximum_violation':float(violations[kind].max()),'bad_indices':bad.tolist(),
                    'bad_count':len(bad),'tolerance':tolerance}
                counts['residual_entries_completed']+=len(points[kind]);guard();save()
            summary['passed']=all(summary[k]['bad_count']==0 for k in ('K','A'))
            counts['residual_batches_completed']+=1
            return values,violations,summary

        _,_,baseline_check=residual(c1,a1,'baseline',base)
        result['baseline_residual']=baseline_check;save()
        if not baseline_check['passed']:raise Rejected('Proven baseline fails unmodified numerical model; no rhs repair.')
        for kind,g in geometry.items():
            g.add(('global_c',),lambda g=g:float(c1)*g.singular)
            g.add(('global_S',),lambda kind=kind:spline[kind])
            for ij in COARSE:g.orbit_column('coarse',ij)
            save();guard()
        spline_mass=sum((H_Q**2*len(orbit(*ij))*v/2 for ij,v in a1.items()),Q(0))
        singular_half=ALPHA*RADIUS**(2-ALPHA)/(2-ALPHA)
        base_cost=float(c1)*singular_half+float(spline_mass)
        objective=[float(c1)*singular_half,float(spline_mass)]+[
            float(COARSE_Q**2*len(orbit(*ij))/2) for ij in COARSE]
        lower=[-.2,-.5]+[-.5]*190;upper=[.2,.5]+[.5]*190
        fine=[];direction_maps=list(coarse_maps)
        matrices={kind:g.matrix() for kind,g in geometry.items()}
        active=[(kind,index) for kind in ('K','A') for index in selected[kind]]
        active_set=set(active);scales=[]
        solver=highspy.Highs()
        for option,value in (('output_flag',False),('threads',1),('solver','simplex'),
            ('simplex_strategy',1),('presolve','off'),('primal_feasibility_tolerance',1e-10),
            ('dual_feasibility_tolerance',1e-10)):
            check_status(solver.setOptionValue(option,value),'option_'+option)
        check_status(solver.addVars(len(objective),np.asarray(lower),np.asarray(upper)),'initial_variables')
        check_status(solver.changeColsCost(len(objective),np.arange(len(objective),dtype=np.int32),np.asarray(objective)), 'cost')

        def raw_rows(row_keys):
            guard()
            matrix=sp.csr_matrix(np.vstack([matrices[k][i] for k,i in row_keys]))
            rhs=np.asarray([targets[k][i]-base[k][i] for k,i in row_keys])
            finite('raw_affine_rows',matrix.data);finite('raw_rhs',rhs)
            return matrix,rhs

        def append_rows(row_keys):
            raw,rhs=raw_rows(row_keys)
            scale=1./np.maximum(1.,abs(raw).max(axis=1).toarray().ravel())
            scaled=(sp.diags(scale)@raw).tocsr()
            check_status(solver.addRows(len(row_keys),rhs*scale,np.full(len(row_keys),highspy.kHighsInf),
                scaled.nnz,scaled.indptr[:-1].astype(np.int32),scaled.indices.astype(np.int32),scaled.data),'add_rows')
            scales.extend(scale.tolist());guard()

        append_rows(active)
        check_status(solver.setSolution(len(objective),np.arange(len(objective),dtype=np.int32),np.zeros(len(objective))), 'zero_primal_start')
        result['zero_primal_start_checked']=True;best=None
        last_python_x=np.empty(0);last_python_objective=np.empty(0)

        def expand(t):
            coeff={ij:(1.+t[1])*float(v) for ij,v in a1.items()}
            for value,direction in zip(t[2:],direction_maps):
                for ij,weight in direction.items():coeff[ij]=coeff.get(ij,0.)+float(value)*float(weight)
            return (1.+t[0])*float(c1),coeff

        def cost(c,coeff):
            return float(c)*singular_half+sum(float(H_Q**2*len(orbit(*ij))/2)*float(v) for ij,v in coeff.items())

        for iteration in range(MAX_LP):
            guard();result['phase']='LP_'+str(iteration+1)
            before_basis=bool(solver.getBasis().valid)
            if iteration and not before_basis:raise Rejected('Extended prior basis invalid; no cold retry.')
            soft=min(LP_SECONDS,SECONDS-(monotonic()-start))
            check_status(solver.setOptionValue('time_limit',soft),'time_limit')
            counts['LP_started']+=1
            record={'iteration':iteration+1,'status':'started','variables':len(objective),'rows':len(active),
                    'prior_extended_basis_valid':before_basis,'same_solver_instance':True,'time_limit_soft_seconds':soft}
            result['solves'].append(record);save();tic=monotonic()
            run_status=solver.run();counts['LP_completed']+=1
            model_status=solver.getModelStatus();solution=solver.getSolution()
            if solution.value_valid and np.all(np.isfinite(solution.col_value)):
                last_python_x=np.asarray(solution.col_value)
                last_python_objective=np.asarray([float(np.asarray(objective)@last_python_x)])
            else:last_python_x=np.empty(0);last_python_objective=np.empty(0)
            info=solver.getInfo();infeasibility=float(info.max_primal_infeasibility)
            record.update(status=solver.modelStatusToString(model_status),run_status=str(run_status),
                wall_seconds=monotonic()-tic,valid_primal=bool(solution.value_valid),
                max_scaled_primal_infeasibility=infeasibility if math.isfinite(infeasibility) else None)
            save();guard()
            if run_status!=highspy.HighsStatus.kOk or model_status!=highspy.HighsModelStatus.kOptimal or not solution.value_valid:
                result['termination']='LP_not_optimal_no_retry';break
            t=finite('direction_solution',np.asarray(solution.col_value))
            bound_error=max(0.,float(np.max(np.asarray(lower)-t)),float(np.max(t-np.asarray(upper))))
            record['variable_bound_error']=bound_error
            if bound_error>1e-10:raise Rejected('Variable bounds violated.')
            c,coeff=expand(t)
            values,violations,check=residual(c,coeff,'LP_'+str(iteration+1))
            record['full_expanded_residual']=check
            record['float_half_cost']=cost(c,coeff)
            record['affine_half_cost']=base_cost+float(np.asarray(objective)@t)
            record['cost_difference']=abs(record['float_half_cost']-record['affine_half_cost'])
            record['affine_full_max_difference']={k:float(np.max(abs(values[k]-(base[k]+matrices[k]@t)))) for k in ('K','A')}
            record['affine_matrix_rows_checked']=sum(len(p) for p in points.values())
            counts['affine_matrix_rows_checked']+=record['affine_matrix_rows_checked']
            active_bad=[(k,i) for k,i in active if violations[k][i]>(K_TOL if k=='K' else A_TOL)]
            record['active_raw_residual_bad_count']=len(active_bad)
            np.savez(ROOT/f'iteration_{iteration+1:02d}.npz',t=t,K=values['K'],A=values['A'])
            dump(ROOT/f'coefficients_{iteration+1:02d}.json',{'c':c,'a':{f'{i},{j}':v for (i,j),v in coeff.items()},'float_half_cost':record['float_half_cost']})
            save();guard()
            if active_bad:raise Rejected('Optimal solver violates active raw rows; no rescaling repair.')
            if check['passed']:
                best=(c,coeff,t.copy(),record['float_half_cost']);result['termination']='fixed_pool_feasible';break
            if iteration==MAX_LP-1:result['termination']='six_LP_cap_not_converged';break
            worst=[]
            for kind,cap,tol in (('K',256,K_TOL),('A',512,A_TOL)):
                omitted=[i for i in check[kind]['bad_indices'] if (kind,i) not in active_set]
                omitted.sort(key=lambda i:(-float(violations[kind][i]),*keys[kind][i]))
                worst.extend((kind,i) for i in omitted[:cap])
            if not worst:result['termination']='no_new_violated_rows';break
            # Deterministic activation, using fixed atom-free spline indices (not sorted coordinate IDs).
            witnesses=sorted(worst,key=lambda row:(-float(violations[row[0]][row[1]])/(K_TOL if row[0]=='K' else A_TOL),
                *keys[row[0]][row[1]],row[0]))
            proposed=[];existing=set(fine);limit=min(24,MAX_FINE-len(fine))
            for kind,index in (witnesses if limit else []):
                u,v=keys[kind][index];iu,iv=int(u//H_Q),int(v//H_Q)
                for i in range(iu-1,iu+3):
                    for j in range(iv-1,iv+3):
                        ij=canonical(i,j)
                        if ij[0]<=N and ij not in existing:
                            proposed.append(ij);existing.add(ij)
                            if len(proposed)>=limit:break
                    if len(proposed)>=limit:break
                if len(proposed)>=limit:break
            old_columns=len(objective)
            for ij in proposed:
                guard()
                for g in geometry.values():g.orbit_column('fine',ij)
                fine.append(ij);direction_maps.append({ij:Q(1)})
                objective.append(float(H_Q**2*len(orbit(*ij))/2));lower.append(-.25);upper.append(.25)
            if proposed:
                matrices={kind:g.matrix() for kind,g in geometry.items()}
                added=sp.csc_matrix(np.vstack([matrices[k][i,old_columns:] for k,i in active])*np.asarray(scales)[:,None])
                check_status(solver.addCols(len(proposed),np.asarray(objective[old_columns:]),
                    np.asarray(lower[old_columns:]),np.asarray(upper[old_columns:]),added.nnz,
                    added.indptr[:-1].astype(np.int32),added.indices.astype(np.int32),added.data),'extend_columns')
            append_rows(worst);active.extend(worst);active_set.update(worst)
            record['separation']={'added_K':sum(k=='K' for k,i in worst),'added_A':sum(k=='A' for k,i in worst),
                'activated_fine_classes':[list(ij) for ij in proposed],
                'basis_valid_after_extensions':bool(solver.getBasis().valid)}
            save();guard()

        result['fine_classes']=[list(ij) for ij in fine]
        result['direction_order']=['relative_c','relative_S']+[f'coarse_{i}_{j}' for i,j in COARSE]+[f'fine_{i}_{j}' for i,j in fine]
        # Export exactly the final existing finite model for a separately authorized single MATLAB solve.
        guard();raw,rhs=raw_rows(active)
        from scipy.io import savemat
        savemat(ROOT/'last_model.mat',{'A_ub':-raw,'b_ub':-rhs[:,None],'c':np.asarray(objective)[:,None],
            'lb':np.asarray(lower)[:,None],'ub':np.asarray(upper)[:,None],
            'python_x':last_python_x[:,None],'python_objective':last_python_objective,
            'base_cost':base_cost},do_compression=False)
        np.savez(ROOT/'finite_model_metadata.npz',objective=objective,lower=lower,upper=upper)
        dump(ROOT/'finite_model_rows.json',{'rows':[[k,i] for k,i in active],
            'row_count':len(active),'column_count':len(objective),'half_cost_constant':base_cost,
            'MATLAB_matrix_residual_checks_maximum':60000,
            'MATLAB_schema':'last_model.mat: A_ub*x<=b_ub; c is affine LP objective, python_objective excludes base_cost; empty python_x if no finite valid primal.'})
        save();guard()
        if best is not None and best[3]<base_cost-GAIN:
            c,coeff,t,value=best;den=2**50
            cq=Q(round(c*den),den);aq={ij:Q(round(v*den),den) for ij,v in coeff.items()}
            if cq<=0:raise Rejected('Rationalized singular coefficient nonpositive.')
            _,_,rational_check=residual(cq,aq,'unique_rationalized_candidate')
            result['rationalized_residual']=rational_check;save();guard()
            if not rational_check['passed']:raise Rejected('Rationalized candidate failed existing pool; no repair.')
            nz={ij:v for ij,v in aq.items() if v}
            actual_U=(max(i for i,j in nz)+2)*H_Q
            if actual_U<R_Q or actual_U>Q(41,16):raise Rejected('Actual exported support scan window invalid.')
            certificate={'alpha':str(ALPHA_Q),'R':str(R_Q),'h':str(H_Q),'N':N,'c':str(cq),
                         'a':{f'{i},{j}':str(v) for (i,j),v in sorted(nz.items())}}
            target=ROOT/'candidate.json'
            if target.exists():raise Rejected('Second rationalized export forbidden.')
            dump(target,certificate)
            result.update(candidate_exported=True,candidate_sha256=sha256(target.read_bytes()).hexdigest(),
                candidate_float_half_cost=cost(cq,aq),support_U=str(actual_U),
                rationalized_coefficient_denominator=str(den),exact_cost_and_full_checkers='Not run by this driver; the lead schedules the already planned exact cost and full-domain checks.')
        result['status']='finite_coarse_fine_discovery_closed'
    except BudgetStop as exc:result.update(status='interrupted_saved_checkpoint',reason=str(exc))
    except Exception as exc:result.update(status='failed_or_rejected',reason=f'{type(exc).__name__}: {exc}')
    assert counts['LP_started']<=6 and counts['LP_completed']<=counts['LP_started']
    assert counts['residual_batches_started']<=8 and counts['residual_entries_started']<=480000
    assert counts['direction_entries_started']<=17280000
    result.update(source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),module_import_seconds_separate=_IMPORT_SECONDS,
        environment={'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'highspy':highspy.Highs().version()},
        scope='Finite numerical model only. No automatic checker, strict cost proof, MATLAB solve, repair or restart.')
    save()
    print(json.dumps({'status':result['status'],'counts':counts,'seconds':result['elapsed_seconds'],
                      'candidate_exported':result['candidate_exported']},ensure_ascii=False),flush=True)


if __name__=='__main__':main()
