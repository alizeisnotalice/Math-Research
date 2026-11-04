"""High-precision numerical evidence, not an interval certificate.
Independent area uses the exact B/C radial widths, with detected crossings split.
The finite four-scale MATLAB implementation is the independent geometric check.
"""
import argparse
import concurrent.futures
import json
import platform
from pathlib import Path
import mpmath as mp


def qfun(mu):
    a = 1 + 2 * mu
    if a <= 0:
        return -mp.mpf(1)/6 - 2*mu
    if a >= mp.mpf('0.5'):
        return mp.mpf(7)/6 + 2*mu
    return mp.mpf(5)/6-a+8*a**3/3


def roots_split(fun, top):
    points = sorted(set([top*j/96 for j in range(1, 97)]
                        + [top*mp.mpf(10)**(-j) for j in range(1,13)]))
    roots = []
    left, fl = points[0], fun(points[0])
    for right in points[1:]:
        fr = fun(right)
        if fl*fr < 0:
            lo, hi, flo = left, right, fl
            for _ in range(230):
                mid = (lo+hi)/2
                fm = fun(mid)
                if fm*flo > 0:
                    lo, flo = mid, fm
                else:
                    hi = mid
            roots.append((lo+hi)/2)
        left, fl = right, fr
    return [mp.mpf(0), *roots, top]


def run_case(case):
    mode, pstr, tstr = case
    mp.mp.dps = 80
    p, t = mp.mpf(pstr), mp.mpf(tstr)
    beta = 1+p*t if mode == 'outer' else 1-t/2+p*t*t
    lam = (beta-1)/t
    alpha = 1-t
    diagonal = max(-mp.expm1(mp.log(alpha)/beta),
                   mp.expm1(-mp.log(alpha)/(2-beta)))/2
    top = diagonal/t
    def wb(u):
        return -2*t*u-mp.expm1(mp.log1p(-t)+(1-beta)*mp.log1p(-2*t*u))
    def wc(u):
        return -mp.expm1(mp.log1p(-t)+(2-beta)*mp.log1p(2*t*u))
    split = roots_split(lambda u: wc(u)-wb(u), top)
    area = 8*t*sum(mp.quad(lambda u: max(wb(u), wc(u), 0), [lo, hi])
                      for lo, hi in zip(split, split[1:]))
    m = max(-lam, 1+lam)
    # The affine absolute-value integral is split exactly at its possible zero.
    a, b = 1+2*lam-2*lam*t, 2*lam*t
    cuts = [mp.mpf(0), mp.mpf('0.5')]
    if b and 0 < -a/b < mp.mpf('0.5'):
        cuts.insert(1, -a/b)
    integral = sum(mp.quad(lambda u: u*abs(a+b*u), [lo, hi])
                   for lo, hi in zip(cuts, cuts[1:]))
    composite = 2*t*t+t**3*(1+8*integral)+2*m*m*t**4
    asym = (2*t*t+max(-2*p, 2+2*p)*t**3 if mode == 'outer'
            else 2*t*t+t**3+qfun(p)*t**4)
    return {'mode': mode, 'parameter': pstr, 't': tstr,
            'area': mp.nstr(area, 55),
            'composite_error_over_t5': float((area-composite)/t**5),
            'window_error_scaled': float((area-asym)/t**(4 if mode=='outer' else 5)),
            'interior_switches_detected': len(split)-2}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=8)
    args = parser.parse_args()
    ts = ['0.02','0.01','0.003','0.001','0.0001']
    params = {'outer':['-2','-1','-0.75','-0.5','-0.25','0','0.5','1','2'],
              'inner':['-1','-0.6','-0.5','-0.45','-0.4','-0.375','-0.3','-0.25','-0.1','0','1']}
    cases = [(mode,p,t) for mode, ps in params.items() for p in ps for t in ts]
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        rows = list(pool.map(run_case, cases))
    summary = {mode: {'cases': sum(r['mode']==mode for r in rows),
                'max_abs_window_scaled_error': max(abs(r['window_error_scaled']) for r in rows if r['mode']==mode)}
               for mode in params}
    result = {'kind':'non-certified-high-precision-numerical-evidence',
              'python':platform.python_version(),'mpmath':mp.__version__,
              'decimal_digits':80,'workers':args.workers,'random_seed':None,
              'crossing_scan':'96 uniform nodes plus 12 geometric nodes near zero','bisection_steps':230,
              'coverage_note':'Finite deterministic parameter grid; crossing scan is not a root-count certificate.',
              'max_abs_composite_error_over_t5':max(abs(r['composite_error_over_t5']) for r in rows),
              'summary':summary,'rows':rows}
    assert all(abs(r['composite_error_over_t5'])<100 for r in rows)
    out=Path(__file__).resolve().parents[1]/'results'/'critical_square_python.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
    matlab_path = out.with_name('critical_square_matlab.json')
    if matlab_path.exists():
        differences = []
        for mode, param, t, area, *_ in json.loads(matlab_path.read_text())['rows']:
            ref = next(r for r in rows if r['mode']==('outer' if mode==1 else 'inner')
                       and float(r['parameter'])==param and float(r['t'])==t)
            differences.append(abs(area/float(ref['area'])-1))
        comparison = {'kind':'cross-tool-numerical-comparison',
                      'matched_cases':len(differences),
                      'max_relative_area_difference':max(differences)}
        out.with_name('critical_square_crosscheck.json').write_text(json.dumps(comparison,indent=2)+'\n')
        print(json.dumps(comparison))

if __name__=='__main__':
    main()
