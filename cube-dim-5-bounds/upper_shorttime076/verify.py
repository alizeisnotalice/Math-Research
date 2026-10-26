#!/usr/bin/env python3
"""Independent normalization/root path; read saved inputs, never invoke producer."""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import hashlib
import json
import os
import signal
import time

HERE = Path(__file__).resolve().parent
D = 2**320


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def alarm(_signum, _frame):
    raise TimeoutError('180-second complete verification limit')


def dyadic_out(lo, hi):
    low = F((lo*D).numerator//(lo*D).denominator, D)
    high = F((hi*D).numerator//(hi*D).denominator+1, D)
    assert low <= lo < hi < high
    return low, high


def newton_root(q):
    n = q.numerator*D**5//q.denominator
    r = 1 << ((n.bit_length()+4)//5)
    while True:
        nxt = (4*r+n//r**4)//5
        if nxt >= r:
            break
        r = nxt
    while r**5 > n:
        r -= 1
    while (r+1)**5 <= n:
        r += 1
    low, high = F(r,D), F(r+1,D)
    assert low**5 <= q < high**5
    return low, high


def interval(row, keys=('lower_dyadic','upper_dyadic')):
    lo,hi = [F(row[k]) for k in keys]
    assert 0 < lo < hi
    return lo,hi


def main():
    start = time.monotonic()
    signal.signal(signal.SIGALRM,alarm)
    signal.alarm(180)
    reg = json.loads((HERE/'registration.json').read_text())
    data = json.loads((HERE/'results.json').read_text())
    assert data['registration_sha256'] == sha(HERE/'registration.json')
    assert data['source_sha256'] == sha(HERE/'certify.py')
    for row in reg['sources']:
        assert sha(row['path']) == row['sha256']
    source = json.loads(Path(reg['sources'][0]['path']).read_text())['constants']['constants']
    p3l,p3u = interval(source['P3'])
    p4l,p4u = interval(source['P4'])
    sl,su = interval(source['S'])
    pi_saved = json.loads(Path(reg['sources'][3]['path']).read_text())['pi_machin80_interval']
    pil,piu = dyadic_out(F(pi_saved[0]),F(pi_saved[1]))
    j = isqrt(5*D*D)
    root5l,root5u = F(j,D),F(j+1,D)
    assert root5l**2 < 5 < root5u**2
    # This path derives F5 from pi and sqrt5 rather than using saved F5.
    fl,fu = 16/(piu*root5u),16/(pil*root5l)
    c3l,c3u = fl*p4l,fu*p4u
    c1l,c1u = 7680*sl*p3l/(piu**2*root5u),7680*su*p3u/(pil**2*root5l)
    old = json.loads(Path(reg['sources'][4]['path']).read_text())
    oldl,oldu = interval(old['E_eps_Lq_norm'],('lower','upper'))
    assert len(data['cases']) == 4
    checks=[]
    for index,row in enumerate(data['cases']):
        eps = F(reg['epsilons'][index])
        assert eps == F(row['epsilon'])
        lo5 = F(3125,7776)*eps**6*(c3l+c1l*eps**2)/piu**4
        hi5 = F(3125,7776)*eps**6*(c3u+c1u*eps**2)/pil**4
        lo,_ = newton_root(lo5)
        _,hi = newton_root(hi5)
        savedL5,savedU5=interval(row['E_fifth_power'],('lower','upper'))
        savedLo,savedHi=interval(row['E_norm'],('lower','upper'))
        assert savedLo**5 <= savedL5 < savedU5 < savedHi**5
        assert max(lo5,savedL5) < min(hi5,savedU5)
        cap=F(reg['preselected_strict_E_caps'][index])
        multiplier=F(reg['preselected_old_E_multipliers'][index])
        assert hi < cap
        assert hi < multiplier*oldl
        below=hi<oldl
        above=lo>oldu
        assert below != above
        expected='strictly_below_old_E' if below else 'strictly_above_old_E'
        assert expected == row['comparison_to_old_E']
        checks.append({'epsilon':str(eps), 'E_lower':str(lo),'E_upper':str(hi),
                       'E_fifth_power_lower':str(lo5),'E_fifth_power_upper':str(hi5),
                       'fixed_strict_cap':str(cap),'cap_passed':True,
                       'old_multiplier':str(multiplier),'multiplier_comparison_passed':True,
                       'comparison_to_old_E':expected,
                       'ratio_new_old_lower':str(lo/oldu),'ratio_new_old_upper':str(hi/oldl)})
    out={'id':'upper_shorttime076_independent_verification',
         'status':'four_caps_and_old_budget_comparisons_verified','pid':os.getpid(),
         'source_sha256':sha(__file__),'input_results_sha256':sha(HERE/'results.json'),
         'registration_sha256':sha(HERE/'registration.json'),
         'independent_root_denominator':str(D),'sqrt5_interval':[str(root5l),str(root5u)],
         'F5_independent_interval':[str(fl),str(fu)],'checks':checks,
         'new_Arb_Gamma_Green_calls':0,'old_E_formula_recomputations':0,
         'seconds':time.monotonic()-start}
    signal.alarm(0)
    (HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['status','pid','seconds']}))


if __name__ == '__main__':
    main()
