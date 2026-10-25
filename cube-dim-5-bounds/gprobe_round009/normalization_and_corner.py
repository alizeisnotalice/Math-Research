"""Exact normalization audit and m=5 corner optimization; no cone-moment rerun."""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import hashlib, json, math, sys, time
import sympy as sp
from flint import arb, ctx

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
sys.path.insert(0, str(PROJECT / 'numerics'))
from body_fan_geometry import base_terms, sorted_chamber_certificate

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def polyval(coeff, x):
    answer = F(0)
    for c in reversed(coeff):
        answer = answer*x + c
    return answer

def interval_poly(coeff, lo, hi):
    lower = upper = F(0)
    for c in reversed(coeff):
        products = (lower*lo, lower*hi, upper*lo, upper*hi)
        lower, upper = min(products)+c, max(products)+c
    return lower, upper

def variations(values):
    signs = [1 if x > 0 else -1 for x in values if x]
    return sum(a != b for a,b in zip(signs, signs[1:]))

def arctan_bounds(q, n=60):
    # n even: alternating partial sum lies below arctan(1/q).
    assert n % 2 == 0
    lo = sum((F((-1)**k, (2*k+1)*q**(2*k+1)) for k in range(n)), F(0))
    return lo, lo + F(1, (2*n+1)*q**(2*n+1))

def display(x):
    return format(float(x), '.17g')

def compact_bounds(lo, hi, digits):
    scale = 10**digits
    return F((lo*scale).__floor__(),scale),F((hi*scale).__ceil__(),scale)

def main():
    began = time.monotonic()
    ctx.prec = 256
    input_path = PROJECT/'gprobe_round008/high_cos_tail.json'
    old_path = PROJECT/'gprobe_round004/moment_m1.json'
    old = json.loads(old_path.read_text())
    data = json.loads(input_path.read_text())
    moment = [F(v) for v in data['rows'][0]['corner']['raw_h_even_0_to_10']]
    geometry = sorted_chamber_certificate()
    raw0 = 3840*sum((F(det, math.prod(rates)) for _,det,rates in base_terms()), F(0))
    assert raw0 == moment[0]
    assert raw0/32 == F(old['M0_times_pi5'])
    assert 2*moment[1]/5 == F(old['coordinate_raw_integral'])
    # P(y)/(32*pi^5) with y=T^-2, and P'(y) = y*Q(y)/2.
    qcoeff = [(-1)**k * moment[k] for k in range(6)]
    pcoeff = [F(0),F(0)] + [qcoeff[k]/(4+2*k) for k in range(6)]
    assert all((k+1)*pcoeff[k+1] == (qcoeff[k-1]/2 if k >= 1 else 0) for k in range(7))
    y = sp.Symbol('y')
    qpoly = sp.Poly(sum(sp.Rational(v.numerator,v.denominator)*y**k for k,v in enumerate(qcoeff)),y)
    sturm = sp.sturm(qpoly.as_expr(), y)
    # SymPy first takes the monic squarefree part; certify it has not
    # changed the root set or multiplicities before accepting its chain.
    assert sp.gcd(qpoly,qpoly.diff()).degree() == 0
    assert sp.Poly(sturm[0],y) == qpoly.monic()
    assert sp.Poly(sturm[1],y) == qpoly.monic().diff()
    for j in range(2,len(sturm)):
        assert sp.Poly(sturm[j],y) == -sp.rem(sp.Poly(sturm[j-2],y),sp.Poly(sturm[j-1],y))
    assert sp.Poly(sturm[-1],y).degree() == 0
    sturm_coeff = [[F(str(v)) for v in reversed(sp.Poly(p,y).all_coeffs())] for p in sturm]
    v0 = variations([c[0] for c in sturm_coeff])
    vinf = variations([c[-1] for c in sturm_coeff])
    assert v0-vinf == 1
    lo,hi = F(17),F(18)
    assert polyval(qcoeff,lo) > 0 > polyval(qcoeff,hi)
    for _ in range(32):
        mid = (lo+hi)/2
        if polyval(qcoeff,mid) > 0: lo=mid
        else: hi=mid
    assert variations([polyval(c,lo) for c in sturm_coeff])-variations([polyval(c,hi) for c in sturm_coeff]) == 1
    # Machin's identity supplies exact rational pi enclosures.
    a,b = arctan_bounds(5),arctan_bounds(239)
    pi_lo,pi_hi = compact_bounds(16*a[0]-4*b[1],16*a[1]-4*b[0],40)
    assert F(314159,100000) < pi_lo < pi_hi < F(314160,100000)
    corner_p = polyval(pcoeff,F(16))
    old_corner_p = F(data['rows'][0]['corner']['values']['5']['1/4']['exact_over_pi5'].split(' / ')[0])
    assert corner_p == old_corner_p
    corner_interval = compact_bounds(corner_p/(32*pi_hi**5),corner_p/(32*pi_lo**5),24)
    pl,pu = interval_poly(pcoeff,lo,hi)
    assert pl > 0
    maximum_interval = compact_bounds(pl/(32*pi_hi**5), pu/(32*pi_lo**5),24)
    target = F(11,100000)
    assert maximum_interval[1] < target
    def arb_fraction(f): return arb(f.numerator)/f.denominator
    out = {
        'status':'exact_normalization_correction_and_corner_obstruction_for_all_T',
        'scope':'Only the odd m=5 cosine-tail minorant is excluded as a global g>=11/100000 certificate. No counterexample to G>=11/100000 is claimed.',
        'm':5, 'fixed_T':'1/4', 'target_g':str(target),
        'normalization':{
            'correct':'1/(32*pi^5) = (2*pi)^-5',
            'round008_display':'1/pi^5; missing factor 1/32',
            'cones':26880, 'base_cones':7, 'hyperoctahedral_group_order':3840,
            'raw0_direct_7_cone_sum_times_3840':str(raw0),
            'raw0_div32_matches_prior_p1_zero_times_pi5':str(raw0/32),
            'raw1_corner_times2_div5_matches_prior_coordinate_second_moment':str(2*moment[1]/5),
            'geometry_certificate':geometry,
        },
        'stored_corner_moments':[str(v) for v in moment],
        'polynomial_variable':'y=T^-2',
        'P_ascending_coefficients':[str(v) for v in pcoeff],
        'Q_ascending_coefficients':[str(v) for v in qcoeff],
        'derivative_identity':'P_prime(y)=y*Q(y)/2',
        'corner_fixed_T':{'correct_exact_over_pi5':str(corner_p/32),
                          'rational_enclosure':[str(v) for v in corner_interval],
                          'enclosure_display':[display(v) for v in corner_interval],
                          'arb_ball':str(arb_fraction(corner_p)/(32*arb.pi()**5))},
        'Sturm':{'sequence_ascending_coefficients':[[str(v) for v in c] for c in sturm_coeff],
                 'variations_at_zero':v0,'variations_at_positive_infinity':vinf,
                 'positive_roots':v0-vinf,'bisections':32,
                 'root_y_interval':[str(lo),str(hi)],
                 'T_squared_interval':[str(1/hi),str(1/lo)],
                 'root_y_interval_display':[display(lo),display(hi)]},
        'corner_supremum_all_T':{'rational_enclosure':[str(v) for v in maximum_interval],
                                  'enclosure_display':[display(v) for v in maximum_interval],
                                  'strictly_below_target':True},
        'pi_enclosure':{'method':'Machin 16 atan(1/5)-4 atan(1/239), 60-term alternating rational bounds',
                        'lower':str(pi_lo),'upper':str(pi_hi)},
        'budget':{'full_5d_Bernstein_tests':0,'branch_boxes':0,'Sturm_bisections':32,
                  'recomputed_26880_cone_moments':False,'unresolved_cover':None,
                  'reason_no_5d_cover':'The requested lower bound on this minorant is false at the corner, for every T>0.'},
        'seconds':time.monotonic()-began,
        'inputs_sha256':{str(p.relative_to(PROJECT)):sha(p) for p in (input_path,old_path,PROJECT/'gprobe_round008/high_cos_tail.py',PROJECT/'numerics/body_fan_geometry.py')},
        'source_sha256':sha(Path(__file__)),
        'runtime':{'python':sys.executable,'sympy':sp.__version__},
    }
    (ROOT/'normalization_and_corner.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('status','corner_fixed_T','corner_supremum_all_T','seconds','source_sha256')},indent=2))

if __name__ == '__main__': main()
