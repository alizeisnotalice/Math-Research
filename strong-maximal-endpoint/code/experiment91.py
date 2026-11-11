"""r091 exact phase-stripe tests and scalar price diagnostics.

Exact Fraction calculations cover the two-dimensional phase integral and box
clipping. mpmath rows are numerical scale checks only; they do not certify a
physical witness theorem or a finite-N onset.
"""
import os
for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_name] = '1'

import json
import math
import multiprocessing as mp
import time
from fractions import Fraction as Q
from pathlib import Path
import sys

import mpmath as m

HERE = Path(__file__).resolve().parent
T_FALLBACK = Path('/Volumes/data/CodexProjects/T')
CANDIDATES = (HERE/'code', HERE,
              T_FALLBACK/'Math-Research-strong-maximal/strong-maximal-endpoint/code')
CODE = next((p for p in CANDIDATES if (p/'verify85.py').is_file() and
             (p/'sliver79.py').is_file()), None)
if CODE is None:
    raise ImportError('could not locate verify85.py and sliver79.py in local code or T fallback')
sys.path.insert(0, str(CODE))
from verify85 import constants85
from sliver79 import logi

STRIPE_P = (Q(1, 4), Q(1, 3), Q(1, 2), Q(3, 4))
STRIPE_M = (8, 16, 32, 64)
UNIT_M = (6, 12, 24, 48)
SCALAR_K = (64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048, 4096)
PRECISIONS = (80, 160)
BOXES = (
    (Q(0), Q(1), Q(0), Q(1)),
    (Q(1), Q(2), Q(0), Q(1, 2)),
)
RECEIVER = (Q(3, 2), Q(101, 200))
RECTS = (
    ('R1', (Q(1), Q(2), Q(0), Q(51, 100))),
    ('R2', (Q(1), Q(2), Q(0), Q(13, 25))),
    ('R3', (Q(11, 10), Q(17, 10), Q(1, 7), Q(13, 25))),
    ('R4', (Q(7, 6), Q(19, 11), Q(2, 9), Q(53, 100))),
)
UNIT_BOXES = ((Q(0), Q(1), Q(0), Q(1)),)
UNIT_RECT = (Q(1, 8), Q(12, 11), Q(1, 7), Q(5, 7))
UNIT_RECEIVER = (Q(12, 11), Q(1, 2))


def _area(rect):
    return (rect[1]-rect[0])*(rect[3]-rect[2])


def _intersection(a, b):
    x0, x1 = max(a[0], b[0]), min(a[1], b[1])
    y0, y1 = max(a[2], b[2]), min(a[3], b[3])
    if x0 >= x1 or y0 >= y1:
        return None
    return x0, x1, y0, y1


def _positive_part_square(z):
    return max(Q(0), z)**2/2


def phase_low_area(rect, eps, p):
    """Exact area where frac((x+y)/eps) lies in [0,p), by four corners."""
    x0, x1, y0, y1 = rect
    max_t = x1+y1
    last = max_t//eps + 1

    def primitive(t):
        return sum((_positive_part_square(t-j*eps) -
                    _positive_part_square(t-(j+p)*eps)
                    for j in range(int(last)+1)), Q(0))

    return (primitive(x1+y1)-primitive(x1+y0)-
            primitive(x0+y1)+primitive(x0+y0))


def _clipped_low_area(rect, boxes, eps, p):
    total = Q(0)
    for box in boxes:
        part = _intersection(rect, box)
        if part is not None:
            total += phase_low_area(part, eps, p)
    return total


def _phi2_high_interval():
    log_iv = logi(Q(5, 4))
    return Q(5, 4)*(1+log_iv[0]), Q(5, 4)*(1+log_iv[1])


def _rectangle_price_interval(low_density, base_density):
    phi_lo, phi_hi = _phi2_high_interval()
    return (low_density+(base_density-low_density)*phi_lo,
            low_density+(base_density-low_density)*phi_hi)


def _global_price_interval(p, total_area):
    phi_lo, phi_hi = _phi2_high_interval()
    return (total_area*(p+(1-p)*phi_lo),
            total_area*(p+(1-p)*phi_hi))


def _mp_price(p, total_area, local_low, local_base, wp):
    with m.workdps(wp):
        pp = m.mpf(p.numerator)/p.denominator
        phi_high = m.mpf(5)/4*(1+m.log(m.mpf(5)/4))
        global_price = total_area*(pp+(1-pp)*phi_high)
        local_price = local_low+(local_base-local_low)*phi_high
        return {
            'working_precision': wp,
            'global_price': m.nstr(global_price, 65),
            'local_rectangle_price': m.nstr(local_price, 65),
        }


def _rect_row(rect_id, rect, boxes, p, eps, receiver, K0, must_witness, total_area):
    area = _area(rect)
    base_area = sum((_area(part) for box in boxes
                     if (part := _intersection(rect, box)) is not None), Q(0))
    low_area = _clipped_low_area(rect, boxes, eps, p)
    base = base_area/area
    low = low_area/area
    discrepancy = abs(low-p*base)
    distances = []
    for box in boxes:
        dx = max(Q(0), box[0]-receiver[0], receiver[0]-box[1])
        dy = max(Q(0), box[2]-receiver[1], receiver[1]-box[3])
        distances.append(max(dx, dy))
    bx = min(distances)
    if bx <= 0:
        raise ArithmeticError(f'receiver is not strictly outside support for {rect_id}')
    proof_bound = 2*K0*eps/bx
    mean = Q(5, 4)*base-Q(1, 4)*low
    strict_margin = mean-1
    if not (rect[0] < receiver[0] < rect[1] and rect[2] < receiver[1] < rect[3]):
        raise ArithmeticError(f'receiver is not interior to {rect_id}')
    witness = strict_margin > 0
    if must_witness and not witness:
        raise ArithmeticError(f'expected actual strict witness failed for {rect_id}, p={p}, eps={eps}')
    price_iv = _rectangle_price_interval(low, base)
    return {
        'id': rect_id,
        'rect': [str(v) for v in rect],
        'receiver': [str(v) for v in receiver],
        'boxes': [[str(v) for v in b] for b in boxes],
        'b_x': str(bx),
        'base_density': str(base),
        'low_density': str(low),
        'density_error': str(discrepancy),
        'proof_bound_2K0eps_over_bx': str(proof_bound),
        'actual_mean_1p25base_minus_0p25low': str(mean),
        'strict_margin': str(strict_margin),
        'actual_witness': witness,
        'rectangle_price_interval': [str(price_iv[0]), str(price_iv[1])],
        'precision_rows': [_mp_price(p, total_area, low, base, wp) for wp in PRECISIONS],
    }


def stripe_job(arg):
    p_s, m_value = arg
    p = Q(p_s)
    eps = Q(1, m_value)
    total_area = sum((_area(b) for b in BOXES), Q(0))
    shares = []
    for i, box in enumerate(BOXES):
        low = phase_low_area(box, eps, p)
        share = low/_area(box)
        if share != p:
            raise ArithmeticError(f'global phase share mismatch p={p}, m={m_value}, box={i}')
        shares.append({'box_index': i, 'box': [str(v) for v in box],
                       'low_share': str(share), 'box_area': str(_area(box))})
    rect_rows = []
    for rid, rect in RECTS:
        rect_rows.append(_rect_row(rid, rect, BOXES, p, eps, RECEIVER, 2,
                                   rid in ('R1', 'R2'), total_area))
    phi_iv = logi(Q(5, 4))
    price_interval = _global_price_interval(p, total_area)
    price_density_iv = _global_price_interval(p, Q(1))
    return {
        'p': str(p), 'm': m_value, 'epsilon': str(eps),
        'K0': 2, 'boxes': [[str(v) for v in b] for b in BOXES],
        'total_support_area': str(total_area), 'global_low_shares': shares,
        'phi2_high_log_interval': [str(v) for v in phi_iv],
        'global_price_interval': [str(v) for v in price_interval],
        'global_price_density_interval': [str(v) for v in price_density_iv],
        'global_price_precision_rows': [_mp_price(p, total_area, p, Q(1), wp)
                                        for wp in PRECISIONS],
        'rectangles': rect_rows,
        'pid': os.getpid(),
    }


def unit_audit_job(m_value):
    p = Q(1, 3)
    eps = Q(1, m_value)
    total_area = Q(1)
    rect_area = _area(UNIT_RECT)
    base_area = _area(_intersection(UNIT_RECT, UNIT_BOXES[0]))
    base = base_area/rect_area
    low_area = _clipped_low_area(UNIT_RECT, UNIT_BOXES, eps, p)
    low = low_area/rect_area
    err = abs(low-p*base)
    mean = Q(5, 4)*base-Q(1, 4)*low
    bx = Q(1, 11)
    bound = 2*eps/bx
    expected_errors = {
        6: Q(11, 48960), 12: Q(11, 257040),
        24: Q(0), 48: Q(0),
    }
    expected_means = {
        6: Q(41393, 39168), 12: Q(217327, 205632),
        24: Q(539, 510), 48: Q(539, 510),
    }
    if err != expected_errors[m_value] or mean != expected_means[m_value]:
        raise ArithmeticError(f'audit91 exact unit result mismatch at m={m_value}: {err}, {mean}')
    if phase_low_area(UNIT_BOXES[0], eps, p) != p:
        raise ArithmeticError(f'unit global share mismatch at m={m_value}')
    price_iv = _rectangle_price_interval(low, base)
    return {
        'p': str(p), 'm': m_value, 'epsilon': str(eps), 'K0': 1,
        'boxes': [[str(v) for v in b] for b in UNIT_BOXES],
        'total_support_area': '1', 'rectangles': [{
            'id': 'unit_audit_R', 'rect': [str(v) for v in UNIT_RECT],
            'receiver': [str(v) for v in UNIT_RECEIVER], 'b_x': str(bx),
            'base_density': str(base), 'low_density': str(low),
            'density_error': str(err), 'proof_bound_2K0eps_over_bx': str(bound),
            'actual_mean_1p25base_minus_0p25low': str(mean),
            'strict_margin': str(mean-1), 'actual_witness': mean > 1,
            'audit_expected_density_error': str(expected_errors[m_value]),
            'audit_expected_mean': str(expected_means[m_value]),
            'rectangle_price_interval': [str(price_iv[0]), str(price_iv[1])],
            'precision_rows': [_mp_price(p, total_area, low, base, wp) for wp in PRECISIONS],
        }],
        'pid': os.getpid(),
    }


def _mp_constant(mu_interval):
    start = sum((m.mpf(Q(q).numerator)/Q(q).denominator for q in mu_interval), m.mpf(0))/2
    return m.findroot(lambda x: x-1-m.exp(-x), start)


def scalar_job(arg):
    k, mu_interval = arg
    rows = []
    guard = math.ceil(k*math.log10(2)/2)+80
    for wp in PRECISIONS:
        with m.workdps(wp+guard):
            N_exact = 1 << k
            N = m.mpf(N_exact)
            L = m.mpf(k)*m.log(2)
            mu = _mp_constant(mu_interval)
            sigma = m.sqrt(mu)
            z0 = m.sqrt(L)
            logtail = lambda x: m.log(m.erfc(x/m.sqrt(2))/2)
            z = m.findroot(lambda x: logtail(x)-m.log(sigma*x)+L/2,
                           (z0*m.mpf('.8'), z0*m.mpf('1.1')))
            t0_sq = L/2-4*m.log(L)
            t1_sq = L-8*m.log(L)
            if not (t0_sq > 0 and t1_sq > 0):
                raise ArithmeticError(f'nonpositive tail scale at k={k}')
            t0, t1 = m.sqrt(t0_sq), m.sqrt(t1_sq)
            shift = 1/m.sqrt(L)
            p0, p1 = m.erfc(t0/m.sqrt(2))/2, m.erfc((t1+shift)/m.sqrt(2))/2
            h = L*m.exp(-L/4)
            error1 = L**m.mpf('2.5')*m.exp(-L/2)
            log_p0_over_h = m.log(p0)-m.log(h)
            log_p1_over_error = m.log(p1)-m.log(error1)

            sqrtN = m.sqrt(N)
            d = sigma*sqrtN*z
            tau = sigma*z/sqrtN
            log_phi1 = m.log1p(tau)+(N-1)*m.log1p(m.log1p(tau))
            log_phi2 = m.log1p(2*tau)+(N-1)*m.log1p(m.log1p(2*tau))
            log_phi1_minus_d = log_phi1-d
            log_phi2_minus_2d = log_phi2-2*d
            log_tau = m.log(tau)
            fee_minus_d = log_tau+log_phi2_minus_2d-log_phi1_minus_d
            rows.append({
                'working_precision': wp,
                'guarded_working_precision': m.mp.dps,
                'N': str(N_exact), 'L': m.nstr(L, 65),
                'mu': m.nstr(mu, 65), 'sigma': m.nstr(sigma, 65),
                'z': m.nstr(z, 65), 't0': m.nstr(t0, 65), 't1': m.nstr(t1, 65),
                'shift': m.nstr(shift, 65), 'p0': m.nstr(p0, 65), 'p1': m.nstr(p1, 65),
                'h': m.nstr(h, 65), 'error1': m.nstr(error1, 65),
                'log_p0_over_h': m.nstr(log_p0_over_h, 65),
                'log_p1_over_error': m.nstr(log_p1_over_error, 65),
                't0_over_z': m.nstr(t0/z, 65), 't1_over_z': m.nstr(t1/z, 65),
                'tau': m.nstr(tau, 65), 'd': m.nstr(d, 65), 'log_tau': m.nstr(log_tau, 65),
                'logPhi1_minus_d': m.nstr(log_phi1_minus_d, 65),
                'logPhi2_minus_2d': m.nstr(log_phi2_minus_2d, 65),
                'fee_minus_d': m.nstr(fee_minus_d, 65),
                'fee_minus_d_over_L': m.nstr(fee_minus_d/L, 65),
            })
    return {'k': k, 'N': str(1 << k), 'precision_rows': rows, 'pid': os.getpid()}


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('output_directory')
    args = parser.parse_args()
    outdir = Path(args.output_directory)
    outdir.mkdir(parents=True, exist_ok=True)
    constants = constants85()
    mu_interval = constants['mu_interval']
    start = time.monotonic()
    stripe_jobs = [(str(p), mval) for p in STRIPE_P for mval in STRIPE_M]
    with mp.get_context('spawn').Pool(12) as pool:
        stripes = pool.map(stripe_job, stripe_jobs)
        unit_rows = pool.map(unit_audit_job, UNIT_M)
        scalar_rows = pool.map(scalar_job, [(k, mu_interval) for k in SCALAR_K])
    all_worker_pids = {r['pid'] for r in stripes+unit_rows+scalar_rows}
    actual_workers = len(all_worker_pids)
    data = {
        'pass_all': True,
        'version': m.__version__,
        'configured_workers': 12,
        'actual_workers': actual_workers,
        'worker_pids': sorted(all_worker_pids),
        'working_precision': list(PRECISIONS),
        'scope': 'Exact phase-stripe arithmetic and scalar lower-bound diagnostics only; no finite-N physical onset and no claim that the atom is the entire expected price.',
        'stripe_inputs': {
            'boxes': [[str(v) for v in b] for b in BOXES], 'K0': 2,
            'p_values': [str(p) for p in STRIPE_P], 'm_values': list(STRIPE_M),
            'receiver': [str(v) for v in RECEIVER],
            'rectangle_ids': [rid for rid, _ in RECTS],
            'rectangles': [{'id': rid, 'rect': [str(v) for v in rect]} for rid, rect in RECTS],
        },
        'counts': {'stripe_cases': len(stripes), 'stripe_rectangle_tests': sum(len(x['rectangles']) for x in stripes),
                   'unit_audit_cases': len(unit_rows), 'scalar_dimensions': len(scalar_rows),
                   'scalar_precision_rows': sum(len(x['precision_rows']) for x in scalar_rows)},
        'constants': {'mu_interval': mu_interval},
        'seconds': time.monotonic()-start,
        'stripe_rows': stripes,
        'unit_audit_rows': unit_rows,
        'scalar_rows': scalar_rows,
    }
    (outdir/'experiments91.json').write_text(json.dumps(data, indent=2)+'\n')
    print(json.dumps({k: v for k, v in data.items()
                      if k not in ('stripe_rows', 'unit_audit_rows', 'scalar_rows')}, indent=2))


if __name__ == '__main__':
    main()
