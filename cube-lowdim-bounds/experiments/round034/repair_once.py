#!/usr/bin/env python3
"""One exact 1023/1024 contraction, <=90s, no AS/Psi models or LP.

Reuse all 16272 rational local rows and one 51574-row cached float screen.
Only complete eligibility and a strict integer-root cost gain export a plain kernel.
"""
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
import os
from pathlib import Path
import sys
from time import monotonic
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'round14';SECONDS=90.;LAMBDA=Q(1023,1024)
BASELINE_SHA='21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
INPUTS={
 'build_result.json':'0262b7a6892c8325c42f5978f9c681a200541c9fc89676add1323a75f09c4ca2',
 'directions.json':'98c0e2fcc7f5d448aef33acc24f2b7fe22977f7a2361b20dc7c2ade926a0fd6e',
 'rows.json':'2affd7976754422037d2c1b396523e227b028bb82e67fffa8fb39885e034db84',
 'parameters.json':'68de308160353726ae7827c4320626ed8f44226317969168be9c1429d6f2d753',
 'solve_result.json':'9b292b2f5f40d2bad90ee342b0d872ea1a53dd3124ef0fa09e8f198c13854a43',
 'rejection_verification.json':'df7649a5faab5525ff7526894805449f1f6536ab20649d593dd57d76335926db'}
BASE_MANIFEST_SHA='0ec57b7775f2308a102fa221985dde67e7b6ee8e8945367ee9d21e6e42c8e29e'
LAST_MANIFEST_SHA='786fedf1c84ef759a4e1db400c9723db828ab663b1ec126d819a13c9e1c1cdd8'
POINTS_SHA='513c8d953695d9ddd1c7a097ee9af2ab624f79b875330d02e0d42ad2dcc70ebd'
B1_ARITHMETIC_SHA='1f7b1a0184101a93f5d64e463377ee4e86374660fee3e3c2de7676ff46237159'
COST_SOURCE_SHA='f5de232467e0060ffd3f0495b028cc6328372ef67c086a0a7693d5f38e762b67'
COLUMNS=(2,3,4,5,198,199,201,202,222,244,267,268)


class BudgetStop(Exception):pass


def dump(path,value):
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n');temp.replace(path)


def main():
    if any((ROOT/name).exists() for name in ('result.json','parameters.json','structured_kernel.json','exact_cost.json','candidate.json')):
        raise SystemExit('Existing repair evidence: no second contraction or retry.')
    start=monotonic()
    counts={'contractions':0,'baseline_row_arithmetic_checks':0,'exact_local_rows_checked':0,
        'cache_files_verified':0,'cached_point_entries_started':0,'cached_point_entries_completed':0,
        'cost_enclosures_started':0,'cost_enclosures_completed':0,'structure_checks':0,
        'new_model_calls':0,'new_point_function_calls':0,'LP_calls':0}
    result={'status':'started','candidate_exported':False,'candidate_sha256':None,'lambda':str(LAMBDA),'counts':counts,
        'limits':{'hard_seconds':90,'contractions':1,'exact_local_rows':16272,'cached_points':51574,
                  'integer_root_cost_enclosures':1,'new_models':0,'LP_calls':0},'inputs':[]}
    def guard():
        if monotonic()-start>=SECONDS:raise BudgetStop('internal_90s_deadline')
    def save():result['elapsed_seconds']=monotonic()-start;dump(ROOT/'result.json',result)
    def checked(path,wanted):
        guard();h=sha256()
        with path.open('rb') as stream:
            while True:
                guard();block=stream.read(4*1024*1024)
                if not block:break
                h.update(block)
        if h.hexdigest()!=wanted:raise ValueError('Pinned SHA mismatch: '+path.name)
        result['inputs'].append({'file':str(path.relative_to(ROOT.parent)),'sha256':wanted})
        return path
    def read(path):return json.loads(path.read_text())
    save()
    try:
        inputs={name:read(checked(OLD/name,wanted)) for name,wanted in INPUTS.items()}
        build,rd=inputs['build_result.json'],inputs['directions.json'];rows=inputs['rows.json']['rows']
        if not build['complete'] or not rd['complete'] or not inputs['rows.json']['complete']:
            raise ValueError('Incomplete old safe columns.')
        if build['directions_sha256']!=INPUTS['directions.json'] or build['rows_sha256']!=INPUTS['rows.json']:
            raise ValueError('Old build artifact binding mismatch.')
        source=inputs['solve_result.json'];rejection=inputs['rejection_verification.json']
        if source['candidate_exported'] or source['exact_bad_count']!=1 or source['exact_rows_checked']!=16272:
            raise ValueError('Frozen rejected parameter provenance changed.')
        if rejection['status']!='one_half_margin_rejection_verified_but_original_row_positive' or Q(rejection['original_safe_row'])<=0:
            raise ValueError('Old rejection was not the documented extra-margin failure.')
        directions=rd['directions']
        if len(directions)!=12 or tuple(d['old_column'] for d in directions)!=COLUMNS or len(rows)!=16272:
            raise ValueError('Frozen rows/direction columns changed.')
        oldq=[Q(v) for v in inputs['parameters.json']['q']]
        names=[d['name'] for d in directions]
        if len(oldq)!=26 or inputs['parameters.json']['direction_names']!=names or any(oldq[14+j]!=abs(oldq[2+j]) for j in range(12)):
            raise ValueError('Frozen q/z order invalid.')
        baseline_path=ROOT.parent/'round01/scaled_certificate.json'
        baseline=read(checked(baseline_path,BASELINE_SHA));c1=Q(baseline['c'])
        if (Q(baseline['alpha']),Q(baseline['R']),Q(baseline['h']),baseline['N'])!=(Q(47,40),Q(33,20),Q(1,32),80):
            raise ValueError('Fixed baseline analytic fields changed.')
        if rd['baseline_sha256']!=BASELINE_SHA or inputs['rows.json']['baseline_sha256']!=BASELINE_SHA:
            raise ValueError('Safe rows refer to another baseline.')
        bounds=[Q(d['bound']) for d in directions]
        def box_ok(q):
            return c1/2<=q[0]<=3*c1/2 and Q(1,2)<=q[1]<=Q(3,2) and all(
                abs(q[2+j])<=bounds[j] and q[14+j]==abs(q[2+j])<=bounds[j] for j in range(12))
        if not box_ok(oldq):raise ValueError('Old q outside frozen exact parameter box.')
        q=[(1-LAMBDA)*c1+LAMBDA*oldq[0],(1-LAMBDA)+LAMBDA*oldq[1]]+[LAMBDA*x for x in oldq[2:14]]
        q+=list(map(abs,q[2:14]));counts['contractions']=1
        if not box_ok(q) or any(q[14+j]!=LAMBDA*oldq[14+j] for j in range(12)):
            raise ValueError('Unique exact affine contraction failed box/abs identity.')
        dump(ROOT/'parameters.json',{'lambda':str(LAMBDA),'old_parameters_sha256':INPUTS['parameters.json'],
            'baseline_sha256':BASELINE_SHA,'q':[str(x) for x in q],'old_q':[str(x) for x in oldq],
            'direction_names':names,'rounding':'None after contraction; exact Fraction affine combination.'})
        bad=[];minimum=None;worst=None
        for index,row in enumerate(rows):
            guard();rho=Q(row['baseline_slack']);c_col=Q(row['c_column']);s_col=Q(row['s_column'])
            if c1*c_col+s_col!=rho or rho<=0:raise ValueError('Invalid reused positive baseline row.')
            counts['baseline_row_arithmetic_checks']+=1
            if len(row['t_midpoints'])!=12 or len(row['t_penalties'])!=12:raise ValueError('Incomplete signed row.')
            penalties=[Q(x) for x in row['t_penalties']]
            if any(x<0 for x in penalties):raise ValueError('Negative signed remainder penalty.')
            slack=q[0]*c_col+q[1]*s_col+sum((q[2+j]*Q(row['t_midpoints'][j])-q[14+j]*penalties[j] for j in range(12)),Q(0))-rho/2
            if slack<0:bad.append(index)
            if minimum is None or slack<minimum:
                minimum=slack;worst={'row_index':index,'leaf_node_id':row['leaf_node_id'],'bernstein_index':row['bernstein_index']}
            counts['exact_local_rows_checked']+=1
        result['exact_local_check']={'bad_count':len(bad),'bad_indices':bad,'minimum_half_margin_slack':str(minimum),'worst':worst}
        save()
        if bad:raise ValueError('Fixed contraction fails exact half-baseline gate; no second lambda.')

        # One numerical screen of the same saved rows; no B/L or full-kernel evaluation.
        import numpy as np
        base_dir=ROOT.parent/'round06';last_dir=ROOT.parent/'round07';cache=base_dir/'cache.local'
        prior=read(checked(base_dir/'result.json',BASE_MANIFEST_SHA))
        last_manifest=read(checked(last_dir/'local_artifact_manifest.json',LAST_MANIFEST_SHA))['files']
        point_doc=read(checked(base_dir/'point_manifest.json',POINTS_SHA))
        base_files={x['file']:x['sha256'] for x in prior['local_cache_files']}
        qf=np.array([float(x) for x in q[:14]])
        result['cached_point_checks']={};point_bad=0
        for kind,n,tolerance in (('K',26604,1e-9),('A',24970,1e-8)):
            for suffix in ('singular','baseline'):
                name=f'cache.local/{kind}_{suffix}.npy';checked(base_dir/name,base_files[name]);counts['cache_files_verified']+=1
            matrix_name=f'{kind}_directions_288.local.npy'
            checked(last_dir/matrix_name,last_manifest[matrix_name]['sha256']);counts['cache_files_verified']+=1
            matrix=np.load(last_dir/matrix_name,mmap_mode='r');singular=np.load(cache/f'{kind}_singular.npy',mmap_mode='r')
            oldbase=np.load(cache/f'{kind}_baseline.npy',mmap_mode='r');points=point_doc[kind]
            if matrix.shape!=(n,288) or singular.shape!=(n,) or oldbase.shape!=(n,) or len(points)!=n:
                raise ValueError('Fixed cached point shape invalid.')
            target=np.asarray([float(Q(u)+Q(v)<=1) for u,v in points]) if kind=='K' else np.zeros(n)
            slack=oldbase-target
            if not np.isfinite(matrix).all() or not np.isfinite(singular).all() or not np.isfinite(slack).all() or np.any(slack<0):
                raise ValueError('Nonfinite cache or negative baseline slack; no clamp.')
            guard();counts['cached_point_entries_started']+=n;save()
            small=np.column_stack([singular,matrix[:,1]]+[matrix[:,d['old_column']] for d in directions])
            values=small@qf
            if not np.isfinite(values).all():raise ValueError('Nonfinite cached finite-row screen.')
            original=target-values;strengthened=target+slack/2-values
            k=int(np.argmax(strengthened));original_bad=int(np.count_nonzero(original>tolerance))
            strengthened_bad=int(np.count_nonzero(strengthened>tolerance));point_bad+=original_bad+strengthened_bad
            result['cached_point_checks'][kind]={'original_bad':original_bad,'strengthened_bad':strengthened_bad,
                'original_max_violation':float(original.max()),'strengthened_max_violation':float(strengthened.max()),
                'tolerance':tolerance,'worst_strengthened_index':k,'worst_point':points[k],
                'scope':'Same cached float point model; no independent kernel/AS reevaluation.'}
            counts['cached_point_entries_completed']+=n;save();guard()
        if point_bad:raise ValueError('Unique cached point screen failed; no extra repair.')

        T={}
        for t,direction in zip(q[2:14],directions):
            guard()
            for key,value in direction['coefficients'].items():
                i,j=map(int,key.split(','))
                if not 0<=j<=i<=80 or i+j>48:raise ValueError('Frozen direction support invalid.')
                T[key]=T.get(key,Q(0))+t*Q(value)
        T={key:value for key,value in T.items() if value}
        coeff={key:q[1]*Q(value) for key,value in baseline['a'].items()}
        for key,value in T.items():coeff[key]=coeff.get(key,Q(0))+value
        coeff={key:value for key,value in coeff.items() if value}
        if q[0]<=0 or q[1]<0 or any(coeff.get(key,Q(0))!=q[1]*Q(baseline['a'].get(key,'0'))+T.get(key,Q(0)) for key in set(coeff)|set(baseline['a'])|set(T)):
            raise ValueError('Exact c/s/expanded structure failed.')
        actual_support=max((Q(sum(map(int,key.split(',')))+4,32) for key in T),default=Q(0))
        if actual_support>Q(3,4):raise ValueError('Selected-family actual inner support escaped 3/4.')
        structured={key:baseline[key] for key in ('alpha','R','h','N')}
        structured.update(beta=baseline['alpha'],c=str(q[0]),a={key:str(value) for key,value in sorted(coeff.items())},
            structure={'s':str(q[1]),'T':{key:str(value) for key,value in sorted(T.items())},'baseline_sha256':BASELINE_SHA})
        structured_path=ROOT/'structured_kernel.json';dump(structured_path,structured)
        cost_path=ROOT/'dependency_enclose_cost.py'
        if not cost_path.exists():cost_path=ROOT.parent/'round10/enclose_cost.py'
        checked(cost_path,COST_SOURCE_SHA)
        mod_spec=importlib.util.spec_from_file_location('round15_reviewed_exact_cost',cost_path)
        C=importlib.util.module_from_spec(mod_spec);sys.modules[mod_spec.name]=C;mod_spec.loader.exec_module(C)
        guard();structure_check=C.structure(structured_path,baseline_path);counts['structure_checks']=1
        counts['cost_enclosures_started']=1;save();exact=C.enclose(structured_path,digits=60)
        counts['cost_enclosures_completed']=1;guard()
        arithmetic_path=ROOT.parent/'round01/scaling_arithmetic.json'
        arithmetic=read(checked(arithmetic_path,B1_ARITHMETIC_SHA))
        if arithmetic['new_certificate_sha256']!=BASELINE_SHA:raise ValueError('B1 exact cost source invalid.')
        b1_lower=Q(arithmetic['new_cost_lower']);upper=Q(exact['upper'])
        exact.update(structure=structure_check,actual_T_support_bound=str(actual_support),lambda_value=str(LAMBDA),
            baseline_lower=str(b1_lower),strict_gain_lower=str(b1_lower-upper),improves_baseline=upper<b1_lower,
            below_user_target=upper<Q(17,5),user_target='17/5',cost_source_sha256=COST_SOURCE_SHA,
            execution_claim='35909/10000',below_execution_claim=upper<Q(35909,10000),
            scope='Exact cost and structure only; full-domain positivity still requires both independent checkers.')
        dump(ROOT/'exact_cost.json',exact);save()
        if not upper<b1_lower:raise ValueError('Unique exact cost enclosure does not improve B1.')
        if not upper<Q(35909,10000):raise ValueError('Unique exact cost exceeds fixed checker claim 3.5909.')
        plain={key:structured[key] for key in ('alpha','R','h','N','c','a')}
        if any(plain[key]!=structured[key] for key in plain):raise ValueError('Plain/structured kernel fields mismatch.')
        guard();dump(ROOT/'candidate.json',plain)
        candidate_sha=sha256((ROOT/'candidate.json').read_bytes()).hexdigest()
        exact['plain_candidate_sha256']=candidate_sha;dump(ROOT/'exact_cost.json',exact)
        result.update(status='unique_fixed_contraction_eligible',candidate_exported=True,candidate_sha256=candidate_sha,
            structured_kernel_sha256=sha256(structured_path.read_bytes()).hexdigest(),
            parameters_sha256=sha256((ROOT/'parameters.json').read_bytes()).hexdigest(),
            exact_cost_sha256=sha256((ROOT/'exact_cost.json').read_bytes()).hexdigest(),
            next='Lead runs the two planned full-domain checks; this is not yet a new C2 bound.')
    except BudgetStop as exc:result.update(status='interrupted_saved_checkpoint',reason=str(exc))
    except Exception as exc:result.update(status='failed_or_rejected',reason=f'{type(exc).__name__}: {exc}')
    assert counts['contractions']<=1 and counts['exact_local_rows_checked']<=16272
    assert counts['cached_point_entries_started']<=51574 and counts['cost_enclosures_started']<=1
    result['source_sha256']=sha256(Path(__file__).read_bytes()).hexdigest();save()
    print(json.dumps({'status':result['status'],'counts':counts,'candidate_exported':result['candidate_exported'],
        'seconds':result['elapsed_seconds']},ensure_ascii=False),flush=True)


if __name__=='__main__':main()
