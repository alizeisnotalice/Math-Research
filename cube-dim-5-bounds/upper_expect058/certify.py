#!/usr/bin/env python3
"""One registered exact-arithmetic batch; no samples or quadrature reruns."""
from pathlib import Path
from fractions import Fraction as F
import hashlib, json, math, os, time

HERE = Path(__file__).resolve().parent
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def put(p, obj):
    tmp = p.with_suffix(p.suffix + '.tmp')
    tmp.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n')
    tmp.replace(p)
def s(x): return str(x)
def dec(x, digits=45):
    sign = '-' if x < 0 else ''
    x = abs(x); n = x.numerator * 10**digits // x.denominator
    return sign + str(n // 10**digits) + '.' + str(n % 10**digits).zfill(digits)

def main():
    start = time.monotonic()
    if (HERE/'results.json').exists():
        print('Existing results left unchanged; use verify.py.'); return
    reg = json.loads((HERE/'registration.json').read_text())
    assert sha(__file__) == reg['code_sha256']['certify.py']
    for p, h in reg['sources'].items(): assert sha(p) == h, p
    budget = json.loads(Path(reg['budget_source']).read_text())
    pi_lo, pi_hi = map(F, budget['pi_machin80_interval'])
    assert F(3) < pi_lo < pi_hi < F(22,7)

    # All Gamma values are exact rational multiples of sqrt(pi).
    g8 = math.prod(F(2*j-1,2) for j in range(1,9))
    g10 = math.prod(F(2*j-1,2) for j in range(1,11))
    factorial7 = math.factorial(7)
    ccoef = 16*g8/factorial7
    assert g8 == F(2027025,256)
    assert g10 == F(654729075,1024)
    assert ccoef == F(6435,256)
    # B0^5 = bcoef*pi^4, with L_iso=1/(16384*sqrt(pi)).
    bcoef = ccoef**5 * factorial7 / (32*g10*16384)
    assert bcoef == 2*g8**5/(factorial7**4*g10)
    gatecoef = F(1,10000)**5 * 60**4 * 32**9 / 3125

    power_bounds=[]; iterations=0; scale=10**80
    def root_bracket(label, x, n):
        nonlocal iterations
        target = x.numerator * scale**n // x.denominator
        low, high = 0, 1 << ((target.bit_length()+n-1)//n)
        while high-low>1:
            mid=(low+high)//2
            if mid**n <= target: low=mid
            else: high=mid
            iterations += 1
            if iterations > 4096 or time.monotonic()-start >= 120:
                put(HERE/'checkpoint.json',{'status':'STOP_registered_limit',
                    'pid':os.getpid(),'elapsed_seconds':time.monotonic()-start,
                    'integer_root_iterations':iterations,'completed_brackets':power_bounds,
                    'pending_label':label,'pending_input':s(x),'degree':n,
                    'binary_low':str(low),'binary_high':str(high)})
                raise RuntimeError('Registered finite limit reached; checkpoint saved.')
        lo, hi = F(low,scale), F(low+1,scale)
        assert lo**n <= x < hi**n
        power_bounds.append({'label':label,'input':s(x),'degree':n,
                             'lower':s(lo),'upper':s(hi)})
        return lo,hi

    blo = root_bracket('B0_lower_pi',bcoef*pi_lo**4,5)[0]
    bhi = root_bracket('B0_upper_pi',bcoef*pi_hi**4,5)[1]
    glo = root_bracket('Jgate_lower_pi',gatecoef*pi_lo**7,5)[0]
    ghi = root_bracket('Jgate_upper_pi',gatecoef*pi_hi**7,5)[1]
    ulo = root_bracket('universal_upper_expression_lower_pi',3**5*bcoef**2*pi_lo**8,10)[0]
    uhi = root_bracket('universal_upper_expression_upper_pi',3**5*bcoef**2*pi_hi**8,10)[1]
    gap = blo-ghi
    assert gap > 0
    assert bcoef > gatecoef*pi_hi**3  # fifth-power comparison without roots
    oldlo, oldhi = map(F, budget['J_flag_full_source_interval'])
    assert oldhi < ulo  # new Jensen crude cap is weaker than locked flag cap

    result = {
        'status':'PASS_finite_global_Jensen_mechanism_excluded_not_actual_expectation_excluded',
        'proof_scope':'One complete Dirichlet probability cell; any positive normalized group/time average',
        'C_coefficient_of_pi':s(ccoef),
        'B0_fifth_power_coefficient_of_pi4':s(bcoef),
        'Jgate_fifth_power_coefficient_of_pi7':s(gatecoef),
        'B0_interval':[s(blo),s(bhi)],
        'Jgate_interval':[s(glo),s(ghi)],
        'global_Jensen_cap_expression_interval':[s(ulo),s(uhi)],
        'B0_lower_minus_Jgate_upper':s(gap),
        'fifth_power_gap_coefficient_lower':s(bcoef-gatecoef*pi_hi**3),
        'locked_Jflag_upper':s(oldhi),
        'display':{'B0':[dec(blo),dec(bhi)],'Jgate':[dec(glo),dec(ghi)],
                   'global_Jensen_cap':[dec(ulo),dec(uhi)],'certified_gap_lower':dec(gap)},
        'power_bounds':power_bounds,
        'new_samples':0,'new_fixed_shapes':0,'new_graph_states':0,'quadrature_evaluations':0,
        'new_I_or_C5_bound':False,
        'runtime':{'pid':os.getpid(),'seconds':time.monotonic()-start,
                   'integer_root_iterations':iterations,'arithmetic':'Fraction and Python int only'},
        'sources':reg['sources'], 'code_sha256':sha(__file__)
    }
    for p,h in reg['sources'].items(): assert sha(p)==h, p
    put(HERE/'results.json',result)
    print(json.dumps({k:result[k] for k in ('status','display','runtime')},ensure_ascii=False))

if __name__=='__main__': main()
