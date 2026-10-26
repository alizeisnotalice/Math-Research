"""PREPARED, not authorized until primary reviews the registered code SHA.

One eta=1/10, epsilon=1/512. Read exact saved076 C1/C3/pi endpoints only.
No H/log/gradient/Gamma/Arb, source checker import, or old baseline rerun.
"""
import argparse
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import os
import time

HERE = Path(__file__).resolve().parent


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def interval(obj):
    lo,hi = F(obj['lower']),F(obj['upper'])
    assert 0 < lo <= hi
    return lo,hi


def pair(lo,hi):
    assert 0 < lo <= hi
    return {'lower':str(lo),'upper':str(hi)}


def fifth_roots(lo,hi,bits=256):
    """Exact integer fifth-root bisection, then verify both directed endpoints."""
    D = 1 << bits
    def floor_root(x):
        low,high = 0,D
        while F(high,D)**5 <= x: high *= 2
        while high-low > 1:
            mid = (low+high)//2
            if mid**5*x.denominator <= x.numerator*D**5: low = mid
            else: high = mid
        assert low**5*x.denominator <= x.numerator*D**5 < (low+1)**5*x.denominator
        return low
    lower,upper = F(floor_root(lo),D),F(floor_root(hi)+1,D)
    assert lower**5 <= lo <= hi < upper**5
    return pair(lower,upper)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--reviewed-registration-sha',required=True)
    args = ap.parse_args()
    rp = HERE/'condition_registration.json'
    assert sha(rp) == args.reviewed_registration_sha
    reg = json.loads(rp.read_bytes())
    assert reg['scope'] == 'one_eta1over10_epsilon1over512_finite_new_condition_A_and_R_bounds'
    assert sha(__file__) == reg['code_sha256']
    for row in reg['source_locks']: assert sha(row['path']) == row['sha256']
    out = HERE/'condition_results.json'
    assert not out.exists(),'Preserve completed finite result; no rerun'
    started = time.monotonic()
    data = json.loads(Path(reg['saved076_results']).read_bytes())
    C1L,C1U = interval(data['derived_constants']['C1'])
    C3L,C3U = interval(data['derived_constants']['C3'])
    piL,piU = interval(data['derived_constants']['pi_outward_384_dyadic'])
    eta,eps = F(1,10),F(1,512)
    d = 5-9*eta+5*eta*eta
    assert d == F(83,20) and eta*eta < F(4,5)
    nearL,nearU = C3L/d,C3U/d
    far_epsL,far_epsU = C1L*eps**2/(30*eta**2),C1U*eps**2/(30*eta**2)
    gap = nearL-far_epsU
    assert gap > 0,'Sufficient condition not certified; do not claim new A'
    ceL,ceU = nearL+C1L*eps**2,nearU+C1U*eps**2
    D = 1 << 256
    Cbar = F((ceU.numerator*D+ceU.denominator-1)//ceU.denominator,D)
    assert Cbar >= ceU
    # New076-style bound; no old E expression is recomputed.
    e076L = F(5,6)**5*eps**6*ceL/piU**4
    e076U = F(5,6)**5*eps**6*ceU/piL**4
    # Consume133's proved exact interface, no133 toy/Beta/Gamma work repeated.
    e133L = F(625,524288)*ceL*eps**6
    e133U = F(625,524288)*ceU*eps**6
    result = {'task':'upper_facepeak134','status':'single_condition_and_new_A_R_bounds_certified',
      'pid':os.getpid(),'registration_sha256':sha(rp),'code_sha256':sha(__file__),
      'eta':str(eta),'epsilon':str(eps),'d_eta':str(d),
      'near_tminus3_coefficient':pair(nearL,nearU),
      'far_epsilon_scaled_tminus3_coefficient':pair(far_epsL,far_epsU),
      'strict_condition_gap_lower':str(gap),
      'condition_covers_all_t':'0<t<=1/512 because far/near comparison multiplies t^2',
      'A_formula':'(20/83)*C3*t^-3+C1*t^-1',
      'C_epsilon_new':pair(ceL,ceU),'Cbar_upper_256dyadic':str(Cbar),
      'new076_norm_fifth_power':pair(e076L,e076U),
      'new076_norm':fifth_roots(e076L,e076U),
      'new133_norm_fifth_power':pair(e133L,e133U),
      'new133_norm':fifth_roots(e133L,e133U),
      'new133_bound_with_Cbar_fifth_power':str(F(625,524288)*Cbar*eps**6),
      'new_H_gradient_log_Gamma_Arb_calls':0,'old_baseline_runs':0,
      'new_I_g_or_C5_bound':False,'matlab_calls':0,'seconds':time.monotonic()-started}
    with out.open('x') as f: json.dump(result,f,ensure_ascii=False,sort_keys=True,indent=2); f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k in ['task','status','pid','eta','epsilon','d_eta','seconds','new_H_gradient_log_Gamma_Arb_calls','old_baseline_runs']},ensure_ascii=False))


if __name__ == '__main__': main()
