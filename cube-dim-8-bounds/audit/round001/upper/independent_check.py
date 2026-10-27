#!/usr/bin/env python3
"""Independent n=8 scalar proof; all mathematical certification is Fraction-only.

No author code is imported or executed. The integral polynomial is built by
seven convolutions, log(t) is bounded directly (no range reduction), and the
target is certified by squaring positive rational inequalities, without sqrt.
"""
from fractions import Fraction as Q
from pathlib import Path
import datetime, hashlib, json, resource, sys, time

D = Path(__file__).resolve().parent
DEFAULT = Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube8-terra-scalar/cube-dim-8-bounds/upper/round001/scalar')
S = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def direct_log_bounds(x, count=400):
    # log(x)=2 atanh(z), z=(x-1)/(x+1). For 0<=z<1 the omitted
    # positive tail is <=2 z^(2N+1)/((2N+1)(1-z^2)).
    assert x > 1
    z = (x-1)/(x+1)
    power, subtotal = z, Q(0)
    for k in range(count):
        subtotal += power/Q(2*k+1)
        power *= z*z
    lo = 2*subtotal
    hi = lo+2*power/(Q(2*count+1)*(1-z*z))
    return lo, hi

def multiply_linear(coeffs, q):
    out = [Q(0)]*(len(coeffs)+1)
    for j, value in enumerate(coeffs):
        out[j] += q*value
        out[j+1] -= value/2
    return out

def scalar_coefficients(t):
    # x=1/rho converts I into int_(1/t)^1 (q-x/2)^7 dx/x.
    # Build that polynomial without binomial coefficients or author coeff().
    q = 1+Q(1,2)/t
    poly = [Q(1)]
    for _ in range(7):
        poly = multiply_linear(poly, q)
    expected = [Q(1), Q(-7,2), Q(21,4), Q(-35,8),
                Q(35,16), Q(-21,32), Q(7,64), Q(-1,128)]
    assert poly == [expected[j]*q**(7-j) for j in range(8)]
    nonlog_integral = sum((poly[j]*(1-t**(-j))/j for j in range(1,8)), Q(0))
    a = 1-((1+1/t)/2)**8-4*nonlog_integral
    b = 8-4*poly[0]
    return a, b

def square_root_enclosure(x, steps=220):
    # Bisection using rational comparisons only. This encloses rather than
    # evaluates sqrt, and its invariants are rechecked by exact squares.
    assert x > 0
    lo, hi = Q(0), Q(1)
    while hi*hi <= x:
        hi *= 2
    for _ in range(steps):
        mid = (lo+hi)/2
        if mid*mid <= x:
            lo = mid
        else:
            hi = mid
    assert lo*lo <= x < hi*hi
    return lo, hi

def enclose(t):
    ll, lh = direct_log_bounds(t)
    a,b = scalar_coefficients(t)
    blo,bhi = sorted([a+b*ll,a+b*lh])
    assert blo>0
    sl = square_root_enclosure(blo)[0]
    sh = square_root_enclosure(bhi)[1]
    base = (1+1/t)**4
    flo,fhi = (base+sl)**2,(base+sh)**2
    return dict(t=t, B_lower=blo, B_upper=bhi, F_lower=flo,
                F_upper=fhi, base=base, log_width=lh-ll,
                log_coefficient=b)

def decimal_outward(x, places, upper=False):
    scale = 10**places
    k = x.numerator*scale//x.denominator
    if upper and Q(k,scale)<x:
        k += 1
    return f'{k//scale}.{k%scale:0{places}d}'

def main():
    registration = (S/'registration.md').read_text()
    assert 't=6,7,8,9,10 in that order' in registration
    source_result = json.loads((S/'result.json').read_text())
    trace = json.loads((S/'trace.json').read_text())
    source_receipt = json.loads((S/'actualreceipt.json').read_text())
    manifest = json.loads((S/'manifest.json').read_text())
    for name,digest in manifest['files'].items():
        assert sha(S/name)==digest, name
    for name,digest in source_receipt['source_sha256'].items():
        assert sha(S/name)==digest, name
    assert source_receipt['status']=='completed' and source_receipt['exit_code']==0
    assert source_receipt['objective_calls']==len(trace)==19
    assert source_result['optimizer_nfev']==14<=80
    assert [Q.from_float(float.fromhex(r['t_hex'])) for r in trace[:5]]==list(map(Q,[6,7,8,9,10]))
    assert all(Q(2)<=Q.from_float(float.fromhex(r['t_hex']))<=Q(32) for r in trace[5:])
    best = min(trace,key=lambda row:float.fromhex(row['F_hex']))
    assert best==source_result['best_observed']
    exact_float = Q.from_float(float.fromhex(best['t_hex']))
    # Python Fraction round is the exact nearest-integer, ties-to-even rule.
    center = round(exact_float*2**20)
    assert center==7674058==source_result['grid_center_integer']
    assert source_result['grid_denominator']==2**20
    # This is a check of saved callbacks, not a new optimizer/search run.
    # The author objective uses ordinary float elementary operations; allow
    # 1e-11 roundoff while independently enclosing the exact scalar at each
    # recorded binary-rational input. This is provenance consistency only.
    for callback in trace:
        callback_t=Q.from_float(float.fromhex(callback['t_hex']))
        callback_f=Q.from_float(float.fromhex(callback['F_hex']))
        true=enclose(callback_t)
        assert true['F_lower']-Q(1,10**11)<callback_f<true['F_upper']+Q(1,10**11)
    ts = [Q(center+j,2**20) for j in (-1,0,1)]
    rows = [enclose(t) for t in ts]
    winner = rows[1]
    assert winner['t']==Q(3837029,524288)
    for j,row in enumerate(source_result['rows']):
        assert Q(row['t'])==rows[j]['t']
        assert row['numerator_on_grid']==center+j-1
        for quantity in ('B','F'):
            ownlo,ownhi = rows[j][quantity+'_lower'],rows[j][quantity+'_upper']
            savedlo,savedhi = Q(row[quantity]['lower']),Q(row[quantity]['upper'])
            # The independently computed rational enclosure contains each
            # author Arb interval. This verifies consistency, not Arb internals.
            assert ownlo<=savedlo<=savedhi<=ownhi
    assert all(winner['F_upper']<other['F_lower'] for other in (rows[0],rows[2]))
    assert source_result['winner']==source_result['rows'][1]
    U = Q(56627260639,2000000000)
    assert Q(source_result['strict_rational_upper'])==U
    saved_u=Q(source_result['winner']['F']['upper'])
    assert U==Q(saved_u.numerator*10**10//saved_u.denominator+1,10**10)
    # Independent target comparison. Since base and B are positive,
    # F=base^2+B+2 base sqrt(B). It is enough to show
    # d=U-base^2-Bhi>0 and d^2>4 base^2 Bhi. No sqrt is used.
    base=winner['base']
    # Compress the proof object outward to an exact 50-place decimal rational.
    # The direct series fractions themselves can have >4300 decimal digits.
    certificate_blo=Q(decimal_outward(winner['B_lower'],50))
    bhi=Q(decimal_outward(winner['B_upper'],50,True))
    assert certificate_blo<=winner['B_lower']<=winner['B_upper']<=bhi
    d=U-base*base-bhi
    residual=d*d-4*base*base*bhi
    assert d>0 and residual>0 and winner['F_upper']<U
    result = {
        'status':'PASS',
        'scope':'independent finite scalar certificate; operator chain reviewed separately in proofreview.md',
        'winner':str(winner['t']), 'strict_upper':str(U),
        'direct_atanh_terms':400,'rational_sqrt_bisection_steps':220,
        'polynomial_method':'seven exact convolutions after x=1/rho',
        'integrated_coefficients':['-7/2','21/8','-35/24','35/64','-21/160','7/384','-1/896'],
        'B_interval_40dp':[decimal_outward(winner['B_lower'],40),decimal_outward(winner['B_upper'],40,True)],
        'F_interval_40dp':[decimal_outward(winner['F_lower'],40),decimal_outward(winner['F_upper'],40,True)],
        'positive_squared_residual_lower_30dp':decimal_outward(residual,30),
        'target_gap_lower_30dp':decimal_outward(U-winner['F_upper'],30),
        'source_manifest_hashes_verified':True,
        'source_arb_intervals_contained_in_independent_enclosures':True,
        'registered_seed_order_trace_count_and_grid_verified':True,
        'saved_callback_values_consistent_to_1e_minus_11':True,
        'strict_three_candidate_winner_verified':True,
        'source_hashes':{p.name:sha(p) for p in sorted(S.iterdir()) if p.is_file()},
        'exact_certificate':{'B_lower':str(certificate_blo),'B_upper':str(bhi),
                             'target_squared_residual':str(residual)},
    }
    (D/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('exact_certificate','source_hashes')},indent=2))

if __name__=='__main__':
    started=datetime.datetime.now(datetime.timezone.utc).isoformat()
    wall,cpu=time.perf_counter(),time.process_time()
    code,error=0,None
    try:
        main()
    except BaseException as exc:
        code,error=1,repr(exc)
        raise
    finally:
        receipt={'status':'PASS' if code==0 else 'FAIL','exit_code':code,'error':error,
                 'start_utc':started,'finish_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                 'wall_seconds':time.perf_counter()-wall,'cpu_seconds':time.process_time()-cpu,
                 'peak_rss_native_units':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                 'python':sys.version,'executable':sys.executable,
                 'source_directory':str(S),'checker_sha256':sha(Path(__file__)),
                 'mathematical_dependencies':'Python standard library; Fraction only'}
        (D/'actualreceipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
