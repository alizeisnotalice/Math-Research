#!/usr/bin/env python3
"""48 scalar point checks at most; no optimization, quadrature or box rerun."""
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys
from time import monotonic

ROOT = Path(__file__).resolve().parent
PRECISION = 256
MAX_CALLS = 48
SECONDS = 180.


class BudgetStop(Exception):
    pass


def dump(path, value):
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def exact_B(x):
    """First adapter: centred B evaluated in exact rational arithmetic."""
    x = abs(Q(x))
    if x <= 1:
        return (4-6*x*x+3*x*x*x)/6
    if x <= 2:
        return (2-x)**3/6
    return Q(0)


class OriginalAdapter:
    def __init__(self, module, path):
        self.module = module
        self.kernel = module.Kern(json.loads(path.read_text()))

    def K(self, u, v):
        M, K = self.module, self.kernel
        spline = M.arb(0)
        for comp in K.comps:
            h = Q(str(comp.h_q))
            bu = {i: exact_B(u/h-i) for i in range(-comp.N, comp.N+1)}
            bv = {j: exact_B(v/h-j) for j in range(-comp.N, comp.N+1)}
            bu = {i: b for i, b in bu.items() if b}
            bv = {j: b for j, b in bv.items() if b}
            value = M.fmpq(0)
            for i, x in bu.items():
                for j, y in bv.items():
                    value += comp.A[i+comp.N][j+comp.N]*M.Q(x)*M.Q(y)
            spline += M.arb(value)
        r, radius = u+v, Q(str(K.R_q))
        if not 0 < r < radius:
            raise ValueError('The frozen diagnostic only uses the smooth inside branch.')
        singular = K.c*(M.apow(M.arb(M.Q(r)), -K.al)-M.apow(K.R, -K.al))
        return singular+spline, {'c_Psi': singular, 'S': spline}

    def AK(self, u, v):
        M, K = self.module, self.kernel
        A1 = M.arb(0); A2 = M.arb(0)
        for comp in K.comps:
            h = Q(str(comp.h_q)); s, t = u/h, v/h
            bu = {i: exact_B(s-i) for i in range(-comp.N, comp.N+1)}
            bv = {j: exact_B(t-j) for j in range(-comp.N, comp.N+1)}
            bu = {i: M.Q(b) for i, b in bu.items() if b}
            bv = {j: M.Q(b) for j, b in bv.items() if b}
            one = M.arb(0); two = M.arb(0)
            for m, row in comp.D.items():
                Dv = sum((row[j+comp.N]*b for j, b in bv.items()), M.fmpq(0))
                Du = sum((row[i+comp.N]*b for i, b in bu.items()), M.fmpq(0))
                # Exact zero gets 0 for positive exponent 3-alpha, never log(0).
                if Dv and s != m:
                    one += M.apow(M.arb(M.Q(abs(s-m))), comp.p)*M.arb(Dv)
                if Du and t != m:
                    two += M.apow(M.arb(M.Q(abs(t-m))), comp.p)*M.arb(Du)
            A1 += comp.sc*one; A2 += comp.sc*two
        r, radius = u+v, Q(str(K.R_q))
        if not 0 < r < radius:
            raise ValueError('The frozen diagnostic only uses the smooth inside branch.')
        ua, va, ra = M.arb(M.Q(u)), M.arb(M.Q(v)), M.arb(M.Q(r))
        Apsi = K.Calpha*M.apow(ra, -2*K.al)+4*K.I0(K.R-ra)+2*K.I0(K.R+ua-va)+2*K.I0(K.R+va-ua)
        singular = K.c*Apsi
        return singular+A1+A2, {'c_A_Psi': singular, 'A1S': A1, 'A2S': A2, 'AS': A1+A2}


class IndependentAdapter:
    def __init__(self, module, path):
        self.module = module
        self.kernel = module.Spec(path)

    def K(self, u, v):
        M, K = self.module, self.kernel
        singular = K.cA*K.Psi(M.A(u+v))
        spline = M.A(K.S_point(u, v))
        return singular+spline, {'c_Psi': singular, 'S': spline}

    def AK(self, u, v):
        K = self.kernel
        one, two = K.AS_point(u, v)
        singular = K.cA*K.APsi_point(u, v)
        return singular+one+two, {'c_A_Psi': singular, 'A1S': one, 'A2S': two, 'AS': one+two}


def main():
    result_path, trace_path = ROOT/'result.json', ROOT/'calls.jsonl'
    if result_path.exists() or trace_path.exists():
        raise SystemExit('Existing diagnostic evidence: no repeated point evaluation.')
    start = monotonic()
    counts = {'started': 0, 'completed': 0, 'failed': 0, 'interrupted': 0, 'not_started': MAX_CALLS}
    result = {'status': 'started', 'counts': counts, 'rows': [], 'precision_bits': PRECISION,
              'maximum_calls': MAX_CALLS, 'call_unit': 'one scalar K or AK at one point/kernel/implementation',
              'independence_scope': 'Two formula implementations sharing python-flint/Arb; no independent definition quadrature.'}
    dump(result_path, result)
    trace = trace_path.open('x')

    def save():
        counts['not_started'] = MAX_CALLS-counts['started']
        result['elapsed_seconds'] = monotonic()-start
        dump(result_path, result)

    def guard():
        if monotonic()-start >= SECONDS:
            raise BudgetStop('180s_internal_deadline')

    try:
        description = json.loads((ROOT/'points.json').read_text())
        points = description['points']
        if len(points) != 12 or description['maximum_calls'] != MAX_CALLS:
            raise ValueError('Frozen point schedule changed.')
        for p in points:
            u, v = Q(p['u']), Q(p['v'])
            if not 0 < v <= u or not 0 < u+v < Q(33, 20):
                raise ValueError('Point on axis/origin/outside fixed smooth branch.')
            if p['quantity'] == 'K' and u+v > 1:
                raise ValueError('K<1 test point outside D.')
        directory = ROOT/'public_src'
        if not (directory/'hlcert.py').is_file():
            directory = ROOT.parent/'upstream-v4/checker'
        original_path, independent_path = directory/'hlcert.py', directory/'independent/hl_core.py'
        first = load('round04_original_hlcert', original_path)
        second = load('round04_independent_hl_core', independent_path)
        from flint import arb, fmpq, ctx
        import flint
        ctx.prec = PRECISION  # Set after BOTH modules' import-time defaults, before constants.
        adapters, kernels = {}, {}
        for item in description['kernels']:
            kernel_path = ROOT/'kernels'/(item['id']+'.json')
            if not kernel_path.is_file(): kernel_path = ROOT.parent/item['file']
            actual_hash = sha256(kernel_path.read_bytes()).hexdigest()
            if actual_hash != item['sha256']:
                raise ValueError('Kernel source hash changed.')
            kernels[item['id']] = actual_hash
            adapters[item['id']] = {
                'original_hlcert_adapter': OriginalAdapter(first, kernel_path),
                'independent_hl_core_point_methods': IndependentAdapter(second, kernel_path)}
            guard()
        result.update(points_sha256=sha256((ROOT/'points.json').read_bytes()).hexdigest(),
                      kernels_sha256=kernels,
                      implementation_source_sha256={'hlcert.py': sha256(original_path.read_bytes()).hexdigest(),
                                                   'hl_core.py': sha256(independent_path.read_bytes()).hexdigest()},
                      flint_version=flint.__version__, python=sys.version)
        save()
        values = {}
        for point in points:
            u, v = Q(point['u']), Q(point['v'])
            threshold_q = Q(point['threshold'])
            threshold = arb(fmpq(threshold_q.numerator, threshold_q.denominator))
            for kernel in description['kernels']:
                for implementation in description['implementations']:
                    guard()
                    if counts['started'] >= MAX_CALLS: raise BudgetStop('48_call_limit')
                    row = {'call': counts['started']+1, 'point_id': point['id'], 'u': point['u'], 'v': point['v'],
                           'quantity': point['quantity'], 'threshold': point['threshold'],
                           'kernel': kernel['id'], 'implementation': implementation, 'status': 'started',
                           'started_seconds': monotonic()-start}
                    counts['started'] += 1; result['rows'].append(row)
                    trace.write(json.dumps(row)+'\n'); trace.flush(); save()
                    try:
                        adapter = adapters[kernel['id']][implementation]
                        value, parts = getattr(adapter, point['quantity'])(u, v)
                        guard()
                        if not value.is_finite() or not all(part.is_finite() for part in parts.values()):
                            raise ValueError('Nonfinite point enclosure: no sign conclusion.')
                        text = value.str(70)
                        recorded = arb(text)  # Classify the outward displayed enclosure saved in the artifact.
                        if recorded.upper() < threshold:
                            verdict = 'strict_condition_violation_at_this_point'
                        elif recorded.lower() >= threshold:
                            verdict = 'condition_certified_at_this_point_only'
                        else:
                            verdict = 'inconclusive_crosses_threshold'
                        row.update(status='completed', ball=text, classification=verdict,
                                   parts={name: part.str(70) for name, part in parts.items()},
                                   completed_seconds=monotonic()-start)
                        values[(point['id'], kernel['id'], implementation)] = recorded
                        counts['completed'] += 1
                    except BudgetStop as exc:
                        counts['interrupted'] += 1
                        row.update(status='interrupted', reason=str(exc)); save()
                        trace.write(json.dumps(row)+'\n'); trace.flush()
                        raise
                    except Exception as exc:
                        counts['failed'] += 1
                        row.update(status='failed_unknown', reason=f'{type(exc).__name__}: {exc}')
                    trace.write(json.dumps(row)+'\n'); trace.flush(); save()
            print(json.dumps({'point_completed': point['id'], 'counts': counts}), flush=True)
        comparisons = []
        for point in points:
            for kernel in description['kernels']:
                keys = [(point['id'], kernel['id'], name) for name in description['implementations']]
                if not all(key in values for key in keys):
                    comparisons.append({'point_id': point['id'], 'kernel': kernel['id'], 'status': 'missing_completed_value'})
                    continue
                a, b = (values[key] for key in keys)
                difference = a-b
                comparisons.append({'point_id': point['id'], 'kernel': kernel['id'],
                                    'intervals_overlap': not (a.upper() < b.lower() or b.upper() < a.lower()),
                                    'difference_contains_zero': not (difference.lower() > 0 or difference.upper() < 0),
                                    'difference_ball': difference.str(70)})
        result.update(status='finite_point_diagnostic_closed', cross_implementation_comparisons=comparisons,
                      scope='Local point facts only; no all-box positivity proof, new weak-type bound, LP, or kernel repair.')
    except BudgetStop as exc:
        result.update(status='interrupted_saved_checkpoint', reason=str(exc))
    except Exception as exc:
        result.update(status='failed_preparation_or_diagnostic', reason=f'{type(exc).__name__}: {exc}')
    finally:
        trace.close()
    assert counts['started'] <= MAX_CALLS
    assert counts['started'] == counts['completed']+counts['failed']+counts['interrupted']
    result.update(source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                  trace_sha256=sha256(trace_path.read_bytes()).hexdigest())
    save()
    print(json.dumps({'status': result['status'], 'counts': counts, 'seconds': result['elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    main()
