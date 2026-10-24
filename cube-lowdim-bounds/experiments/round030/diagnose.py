#!/usr/bin/env python3
"""Three preselected parent-knot centers x two Arb formulas; <=6 calls/180s.

No import-time experiment. Source selection is fixed before any AK evaluation.
Only the lead-owned outer wrapper imposes hard kill during a long Arb call.
"""
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
import math
from pathlib import Path
import sys
from time import monotonic

ROOT=Path(__file__).resolve().parent
PRECISION,MAX_CALLS,SECONDS=192,6,180.
SOURCE_SHA='4a37443a3b5e9e894537e5b1d309f1318239461e1317c471ed9cc34afd7cf2b1'
KERNEL_SHA='e88f1a86ffe0dfb5e12140bd77e2defedc21249ef6514f45400fe5128112b0ad'
ADAPTER_SHA='028f2d38b0fbbd61c890b56d62fa870976b7ae7b56b567c97dbcb691de7b305b'
ORIGINAL_SHA='4bd7e823bdcb336f8e20c188d68d6d8c1e4aa22ec52353885059c24445317db5'
INDEPENDENT_SHA='c3d864b1368ad8150fd3366bf77f687a9a3cdd49b924a20d576aa678aed42ebf'
METHODS=('original_hlcert_adapter','independent_hl_core_point_methods')


class BudgetStop(Exception):pass


def dump(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    temporary.replace(path)


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module)
    return module


def main():
    result_path,trace_path=ROOT/'result.json',ROOT/'calls.jsonl'
    if result_path.exists() or trace_path.exists() or (ROOT/'points.json').exists():
        raise SystemExit('Existing diagnostic evidence: no repeated point evaluation.')
    start=monotonic();current=None
    counts={'started':0,'completed':0,'failed':0,'interrupted':0,'not_started':MAX_CALLS}
    result={'status':'started','precision_bits':PRECISION,'maximum_calls':MAX_CALLS,'counts':counts,
        'call_unit':'One complete scalar AK at one rational point with one formula implementation',
        'rows':[],'new_LP_calls':0,'all_domain_checker_calls':0,
        'independence_scope':'Two reviewed formula implementations sharing python-flint/Arb; no definition quadrature.'}
    dump(result_path,result);trace=trace_path.open('x')
    def guard():
        if monotonic()-start>=SECONDS:raise BudgetStop('180s_internal_deadline')
    def save():
        counts['not_started']=MAX_CALLS-counts['started']
        result['elapsed_seconds']=monotonic()-start;dump(result_path,result)
    def verify(path,wanted):
        guard();actual=sha256(path.read_bytes()).hexdigest()
        if actual!=wanted:raise ValueError('Frozen input SHA mismatch: '+path.name)
        return actual
    def emit(row):trace.write(json.dumps(row,ensure_ascii=False,allow_nan=False)+'\n');trace.flush()
    try:
        source=ROOT.parent/'round10_first.json';kernel_path=ROOT.parent/'round10/candidate_public_schema.json'
        verify(source,SOURCE_SHA);verify(kernel_path,KERNEL_SHA)
        report=json.loads(source.read_text());kernel=json.loads(kernel_path.read_text())
        if (Q(kernel['alpha']),Q(kernel['R']),Q(kernel['h']),kernel['N'])!=(Q(47,40),Q(33,20),Q(1,32),80):
            raise ValueError('Frozen alpha/R/h/N changed.')
        if 'beta' in kernel:raise ValueError('Only the reviewed alpha=beta public schema is allowed.')
        failures=report['A_fail'];ordered=[]
        if not failures:raise ValueError('No recorded A failure boxes.')
        for index,row in enumerate(failures):
            guard()
            if len(row)!=5 or not math.isfinite(float(row[4])):raise ValueError('Malformed/nonfinite failure entry.')
            bounds=tuple(Q(x) for x in row[:4]);lo=Q(str(row[4]))
            if not bounds[0]<bounds[1] or not bounds[2]<bounds[3]:raise ValueError('Empty failure box.')
            ordered.append((lo,bounds,index))
        ordered.sort(key=lambda item:(item[0],item[1],item[2]))
        h=Q(kernel['h']);radius=Q(kernel['R']);selected=[];parents=set();centers=set()
        for lower,box,index in ordered:
            guard();u=(box[0]+box[1])/2;v=(box[2]+box[3])/2
            parent=(int(u//h),int(v//h))
            if parent in parents:continue
            if (u,v) in centers:raise ValueError('Different parents unexpectedly share center.')
            # Validate the selected box; do not skip it and replace with another point.
            i,j=parent
            if not (i*h<=box[0]<box[1]<=(i+1)*h and j*h<=box[2]<box[3]<=(j+1)*h):
                raise ValueError('Selected failure box is not inside its recorded true knot cell.')
            if not (i*h<u<(i+1)*h and j*h<v<(j+1)*h and 0<v<=u and 0<u+v<radius):
                raise ValueError('Selected center is not in the reviewed smooth inside branch.')
            point={'id':f'A{len(selected)+1}','quantity':'AK','threshold':'0','u':str(u),'v':str(v),
                'source_index':index,'source_box':[str(x) for x in box],
                'source_lower_float_as_rational':str(lower),'source_lower_float':failures[index][4],
                'parent_knot_index':list(parent),'parent_knot_box':[str(i*h),str((i+1)*h),str(j*h),str((j+1)*h)],
                'selection':'Worst lower in each of first three distinct h-knot cells; exact rational center.',
                'away_axis_origin_truncation_knots':True}
            selected.append(point);parents.add(parent);centers.add((u,v))
            if len(selected)==3:break
        if len(selected)!=3:raise ValueError('Fewer than three different parent knot cells; no replacement schedule.')
        description={'status':'fixed_before_any_AK_call','source_sha256':SOURCE_SHA,'kernel_sha256':KERNEL_SHA,
            'selection':'Sort by lower, exact box endpoints, source index; take first three distinct floor(center/h) parents.',
            'parent_spacing':str(h),'points':selected,'maximum_calls':MAX_CALLS,'precision_bits':PRECISION}
        dump(ROOT/'points.json',description)
        result.update(source_failure_sha256=SOURCE_SHA,kernel_sha256=KERNEL_SHA,
            points_sha256=sha256((ROOT/'points.json').read_bytes()).hexdigest(),points=selected)
        save();guard()
        adapter_path=ROOT/'dependency_round04_diagnose.py'
        if not adapter_path.exists():adapter_path=ROOT.parent/'round04/diagnose.py'
        verify(adapter_path,ADAPTER_SHA)
        public=ROOT/'public_src'
        if not (public/'hlcert.py').is_file():public=ROOT.parent/'upstream-v4/checker'
        original_path,independent_path=public/'hlcert.py',public/'independent/hl_core.py'
        verify(original_path,ORIGINAL_SHA);verify(independent_path,INDEPENDENT_SHA)
        helper=load('round11_reviewed_round04_adapter',adapter_path)
        original=load('round11_original_hlcert',original_path)
        independent=load('round11_independent_hl_core',independent_path)
        from flint import arb,ctx
        import flint
        ctx.prec=PRECISION  # after both import-time defaults; before either kernel's constants
        adapters={METHODS[0]:helper.OriginalAdapter(original,kernel_path),
                  METHODS[1]:helper.IndependentAdapter(independent,kernel_path)}
        result.update(implementation_source_sha256={'round04_adapter':ADAPTER_SHA,'hlcert.py':ORIGINAL_SHA,
            'hl_core.py':INDEPENDENT_SHA},flint_version=flint.__version__,python=sys.version)
        save();guard();values={}
        for point in selected:
            u,v=Q(point['u']),Q(point['v'])
            for method in METHODS:
                guard()
                if counts['started']>=MAX_CALLS:raise BudgetStop('six_call_cap')
                current={'call':counts['started']+1,'point_id':point['id'],'u':str(u),'v':str(v),
                    'quantity':'AK','threshold':'0','kernel_sha256':KERNEL_SHA,'implementation':method,
                    'source_box':point['source_box'],'source_lower_float':point['source_lower_float'],
                    'source_index':point['source_index'],'source_failure_sha256':SOURCE_SHA,
                    'parent_knot_index':point['parent_knot_index'],'status':'started','started_seconds':monotonic()-start}
                counts['started']+=1;result['rows'].append(current);emit(current);save()
                try:
                    value,parts=adapters[method].AK(u,v)
                    guard()
                    if not value.is_finite() or not all(part.is_finite() for part in parts.values()):
                        raise ValueError('Nonfinite point interval; no sign conclusion.')
                    ball=value.str(70);recorded=arb(ball)
                    if recorded.upper()<0:classification='strict_negative_AK_at_this_point'
                    elif recorded.lower()>=0:classification='nonnegative_AK_at_this_point_only'
                    else:classification='inconclusive_interval_crosses_zero'
                    current.update(status='completed',ball=ball,
                        lower_endpoint_ball=recorded.lower().str(70),upper_endpoint_ball=recorded.upper().str(70),
                        classification=classification,parts={name:part.str(70) for name,part in parts.items()},
                        strict_negative_neighborhood_justified=classification=='strict_negative_AK_at_this_point',
                        completed_seconds=monotonic()-start)
                    values[point['id'],method]=recorded;counts['completed']+=1;emit(current);save();current=None
                except BudgetStop as exc:
                    counts['interrupted']+=1;current.update(status='interrupted',reason=str(exc));emit(current);save();current=None;raise
                except Exception as exc:
                    counts['failed']+=1;current.update(status='failed_unknown',reason=f'{type(exc).__name__}: {exc}')
                    emit(current);save();current=None;raise
        comparisons=[]
        for point in selected:
            a,b=(values[point['id'],method] for method in METHODS);difference=a-b
            comparisons.append({'point_id':point['id'],
                'intervals_overlap':not (a.upper()<b.lower() or b.upper()<a.lower()),
                'difference_contains_zero':not (difference.lower()>0 or difference.upper()<0),
                'difference_ball':difference.str(70)})
        result.update(status='six_fixed_Arb_point_calls_closed',cross_implementation_comparisons=comparisons,
            continuity_argument='Selected centers have positive distance from axes, origin, r=R and knots. Apsi is locally analytic; AS is a finite sum of continuous B and |x/h-m|^(3-alpha). A strict upper<0 therefore implies an open negative neighborhood of positive area, violating a.e. AK>=0 through the reviewed away-origin distribution bridge. No neighborhood size computed.',
            scope='Point facts only; nonnegative selected centers do not prove their boxes or the full domain nonnegative.')
    except BudgetStop as exc:result.update(status='interrupted_saved_checkpoint',reason=str(exc))
    except Exception as exc:result.update(status='failed_preparation_or_diagnostic',reason=f'{type(exc).__name__}: {exc}')
    finally:trace.close()
    assert counts['started']<=MAX_CALLS
    assert counts['started']==counts['completed']+counts['failed']+counts['interrupted']
    result.update(source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        trace_sha256=sha256(trace_path.read_bytes()).hexdigest())
    save();print(json.dumps({'status':result['status'],'counts':counts,'seconds':result['elapsed_seconds']},ensure_ascii=False),flush=True)


if __name__=='__main__':main()
