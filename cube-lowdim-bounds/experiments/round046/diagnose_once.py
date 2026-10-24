#!/usr/bin/env python3
"""Two frozen inner points x two reviewed Arb AK paths, <=4 calls/60s.

Static preparation only until the lead launches its owned hard-timeout wrapper.
No import-time experiment, kernel search, box model, or whole-domain checker.
"""
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
import math
from pathlib import Path
import sys
from time import monotonic

ROOT = Path(__file__).resolve().parent
PRECISION, MAX_CALLS, SECONDS = 192, 4, 60.0
PINS = {
    '../round25/candidate.json': '7081cfef6fb47e5d47c95d3740c7a6ba48d48bb6197928a84d8cb6340c4166d5',
    '../round26/first_report.json': 'cde8d4c4547a706d420df92a0c776643aaa86d82a71a7144cf9a1d584b506b77',
    '../round26/second_report.json': '755071dc5049e885f239e2eacdd83d6247ab4796d55889de7c25fb60c3972c76',
    '../round04/diagnose.py': '028f2d38b0fbbd61c890b56d62fa870976b7ae7b56b567c97dbcb691de7b305b',
    '../round17/fixed_checker/hlcert.py': '4c347b7782120f64e20df07a44b51453feff8eb2e3ed9ab4abfec68f2776fa05',
    '../upstream-v4/checker/independent/hl_core.py': 'c3d864b1368ad8150fd3366bf77f687a9a3cdd49b924a20d576aa678aed42ebf',
    '../upstream-v4/checker/independent/hl_check.py': '0a5ee85e9b6c250494aaf88abce55bb28a4ba5fc46075c531eeefaa0cc2888f8',
}
METHODS = ('corrected_first_exact_B_adapter', 'independent_core_point_methods')
POINTS = ((Q(23497, 65536), Q(10819, 65536)),
          (Q(47105, 131072), Q(21905, 131072)))


class BudgetStop(Exception):
    pass


def dump(path, value):
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
    temp.replace(path)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def exact_endpoint(value):
    """Only a finite exact endpoint Arb is passed, never an approximate midpoint."""
    if not value.is_finite() or not value.is_exact():
        raise ValueError('Endpoint is not exact and finite.')
    mantissa, exponent = value.man_exp()
    m, e = int(mantissa), int(exponent)
    return Q(m * 2**e) if e >= 0 else Q(m, 2**(-e))


def enclosure(value):
    if not value.is_finite():
        raise ValueError('Nonfinite ball; no sign conclusion.')
    lo, hi = exact_endpoint(value.lower()), exact_endpoint(value.upper())
    if lo > hi:
        raise ValueError('Reversed outward endpoints.')
    return {'ball': value.str(70), 'lower': str(lo), 'upper': str(hi)}


def frozen_points(first, second):
    failures = first['A_fail']
    if not failures:
        raise ValueError('No first A_fail list.')
    for row in failures:
        if len(row) != 5 or not math.isfinite(float(row[-1])):
            raise ValueError('Malformed/nonfinite saved failure row.')
    index, first_box = min(enumerate(failures), key=lambda pair: (Q(str(pair[1][-1])), pair[0]))
    bounds1 = tuple(Q(x) for x in first_box[:4])
    if bounds1 != (Q(2937,8192), Q(11749,32768), Q(5409,32768), Q(2705,16384)):
        raise ValueError('Frozen first worst box changed.')
    parents = [row for row in second['c4']['results'] if not row['ok']]
    if len(parents) != 1 or (parents[0]['p'], parents[0]['q']) != (11,5):
        raise ValueError('Second failed parent changed.')
    saved = parents[0]['fails'][0]
    bounds2 = tuple(Q(x) for x in saved[:4])
    if bounds2 != (Q(1,2), Q(1025,2048), Q(89,256), Q(713,2048)):
        raise ValueError('Frozen second first saved box changed.')
    centers = (((bounds1[0]+bounds1[1])/2, (bounds1[2]+bounds1[3])/2),
               ((Q(11)+(bounds2[0]+bounds2[1])/2)/32,
                (Q(5)+(bounds2[2]+bounds2[3])/2)/32))
    if centers != POINTS:
        raise ValueError('Frozen rational centers do not match saved reports.')
    points = []
    for i, (u,v) in enumerate(centers):
        if not (0 < v <= u and u+v < Q(33,20) and int(u*32)==11 and int(v*32)==5):
            raise ValueError('Point is not in the reviewed inner branch/parent.')
        points.append({'id': 'P'+str(i+1), 'u': str(u), 'v': str(v), 'r': str(u+v),
            'quantity': 'AK', 'threshold': '0', 'parent': [11,5],
            'source': 'first_worst_A_fail' if i==0 else 'second_first_saved_fail',
            'source_index': index if i==0 else 0,
            'source_coordinate_system': 'physical' if i==0 else 'local_h_knot',
            'source_box': [str(x) for x in (bounds1 if i==0 else bounds2)],
            'source_sufficient_lower_float': first_box[-1] if i==0 else saved[-1]})
    return points


def main():
    if any((ROOT/name).exists() for name in ('result.json', 'points.json', 'calls.jsonl')):
        raise SystemExit('Existing point evidence: no repeated evaluation.')
    start = monotonic()
    counts = {'started':0, 'completed':0, 'failed':0, 'interrupted':0, 'not_started':MAX_CALLS}
    result = {'status':'started', 'counts':counts, 'rows':[], 'precision_bits':PRECISION,
        'maximum_calls':MAX_CALLS, 'budget_seconds':SECONDS,
        'call_unit':'One complete scalar AK at one frozen rational point with one implementation',
        'new_K_calls':0, 'new_box_models':0, 'LP_calls':0, 'full_checker_calls':0,
        'independence_scope':'Two reviewed point formula paths sharing Arb and spline analytic identities; no independent definition quadrature.',
        'source_sha256_before':digest(Path(__file__))}
    out, trace_path = ROOT/'result.json', ROOT/'calls.jsonl'
    dump(out,result)
    trace = trace_path.open('x')
    def guard():
        if monotonic()-start >= SECONDS:
            raise BudgetStop('internal_60s_deadline')
    def save():
        counts['not_started']=MAX_CALLS-counts['started']
        result['elapsed_seconds']=monotonic()-start
        dump(out,result)
    def emit(row):
        trace.write(json.dumps(row,ensure_ascii=False,allow_nan=False)+'\n'); trace.flush()
    try:
        before = {}
        for relative,wanted in PINS.items():
            guard(); actual=digest(ROOT/relative)
            if actual != wanted:
                raise ValueError('Pinned dependency changed: '+relative)
            before[relative]=actual
        result['input_source_sha256_before']=before; save()
        kernel_path = ROOT/'../round25/candidate.json'
        kernel=json.loads(kernel_path.read_text())
        if set(kernel) != {'alpha','R','h','N','c','a'} or (
            Q(kernel['alpha']),Q(kernel['R']),Q(kernel['h']),kernel['N']) != (Q(47,40),Q(33,20),Q(1,32),80):
            raise ValueError('Frozen public schema or analytic parameters changed.')
        points=frozen_points(json.loads((ROOT/'../round26/first_report.json').read_text()),
                             json.loads((ROOT/'../round26/second_report.json').read_text()))
        dump(ROOT/'points.json',{'status':'fixed_before_any_point_call','points':points,
             'maximum_calls':MAX_CALLS,'precision_bits':PRECISION,'source_pins':PINS})
        result.update(points=points,points_sha256=digest(ROOT/'points.json'));save();guard()
        helper=load('round27_reviewed_adapter',ROOT/'../round04/diagnose.py')
        first=load('round27_fixed_hlcert',ROOT/'../round17/fixed_checker/hlcert.py')
        second=load('round27_independent_core',ROOT/'../upstream-v4/checker/independent/hl_core.py')
        from flint import ctx
        import flint
        ctx.prec=PRECISION  # AFTER module defaults, BEFORE kernel constants/construction.
        adapters={METHODS[0]:helper.OriginalAdapter(first,kernel_path),
                  METHODS[1]:helper.IndependentAdapter(second,kernel_path)}
        result.update(flint_version=flint.__version__,python=sys.version,
            inside_method='First exact B/D with C_alpha*r^(-2alpha)+4I+2I+2I; second AS_point/APsi_point. Both guarded 0<v<=u,r<R; no exterior integral branch.')
        save();guard();values={}
        for point in points:
            u,v=Q(point['u']),Q(point['v'])
            for method in METHODS:
                guard()
                if counts['started'] >= MAX_CALLS:
                    raise BudgetStop('four_scalar_cap')
                row={'call':counts['started']+1,'point_id':point['id'],'u':str(u),'v':str(v),
                     'quantity':'AK','threshold':'0','method':method,'status':'started',
                     'started_seconds':monotonic()-start}
                counts['started']+=1;result['rows'].append(row);emit(row);save()
                try:
                    value,parts=adapters[method].AK(u,v)
                    guard()
                    packed=enclosure(value)
                    packed_parts={name:enclosure(part) for name,part in parts.items()}
                    lo,hi=Q(packed['lower']),Q(packed['upper'])
                    verdict='strict_negative_AK' if hi<0 else ('strict_positive_AK' if lo>0 else 'inconclusive_crosses_or_touches_zero')
                    row.update(status='completed',**packed,parts=packed_parts,classification=verdict,
                        strict_negative_neighborhood_justified=hi<0,completed_seconds=monotonic()-start)
                    values[point['id'],method]=value
                    counts['completed']+=1;emit(row);save()
                except BudgetStop as exc:
                    counts['interrupted']+=1;row.update(status='interrupted',reason=str(exc));emit(row);save();raise
                except Exception as exc:
                    counts['failed']+=1;row.update(status='failed_unknown',reason=f'{type(exc).__name__}: {exc}');emit(row);save();raise
        comparisons=[]
        for point in points:
            a,b=(values[point['id'],method] for method in METHODS)
            comparisons.append({'point_id':point['id'],
                'intervals_overlap':not(a.upper()<b.lower() or b.upper()<a.lower()),
                'difference_contains_zero':(a-b).contains(0), 'difference':enclosure(a-b)})
        result.update(status='four_fixed_AK_calls_closed',cross_implementation_comparisons=comparisons,
            continuity_argument='Both points are away from axes/origin/r=R. APsi is locally continuous and finite AS consists of continuous cubic B and positive powers |x/h-m|^(3-alpha), including at knots. Strict AK upper<0 therefore yields an open negative neighborhood and excludes a.e./distribution nonnegativity via the reviewed bridge. No neighborhood radius computed.',
            scope='Two point facts only. Positive points do not certify boxes, parent, or full kernel.')
    except BudgetStop as exc:
        result.update(status='interrupted_saved_checkpoint',reason=str(exc))
    except Exception as exc:
        result.update(status='failed_preparation_or_diagnostic',reason=f'{type(exc).__name__}: {exc}')
    finally:
        trace.close()
        try:
            after={relative:digest(ROOT/relative) for relative in PINS}
            result['input_source_sha256_after']=after
            result['inputs_sources_unchanged']=after==result.get('input_source_sha256_before')
            result['source_sha256_after']=digest(Path(__file__))
            result['source_unchanged']=result['source_sha256_after']==result['source_sha256_before']
            if not result['inputs_sources_unchanged'] or not result['source_unchanged']:
                result['status']='provenance_changed_invalid'
        except Exception as exc:
            result.update(status='provenance_recheck_failed',provenance_error=str(exc))
        result['trace_sha256']=digest(trace_path)
        save()
    assert counts['started']<=MAX_CALLS
    assert counts['started']==counts['completed']+counts['failed']+counts['interrupted']
    print(json.dumps({'status':result['status'],'counts':counts,'seconds':result['elapsed_seconds']}),flush=True)
    return 0 if result['status']=='four_fixed_AK_calls_closed' else 1


if __name__=='__main__':
    raise SystemExit(main())
