"""Non-certified high-precision checks of polynomial envelopes and side windows.
The rigorous error bound is analytic and stated in the manuscript.
"""
import concurrent.futures
import json
import platform
import os
import resource
import time
from pathlib import Path
import mpmath as mp


def integrate_envelope(functions, top):
    functions = [*functions, lambda u: mp.mpf(0)]
    nodes = sorted(set([top*j/64 for j in range(1,65)]
                       + [top*mp.mpf(10)**(-j) for j in range(1,13)]))
    values = [[f(u) for f in functions] for u in nodes]
    roots = [mp.mpf(0), top]
    for i in range(len(functions)):
        for j in range(i):
            differences = [v[i]-v[j] for v in values]
            if max(abs(d) for d in differences) < mp.mpf('1e-65'):
                continue
            for k in range(len(nodes)-1):
                fl, fr = differences[k:k+2]
                if fl == 0:
                    roots.append(nodes[k])
                if fr == 0:
                    roots.append(nodes[k+1])
                if fl*fr < 0:
                    lo, hi = nodes[k], nodes[k+1]
                    for _ in range(220):
                        mid = (lo+hi)/2
                        fm = functions[i](mid)-functions[j](mid)
                        if fl*fm > 0:
                            lo, fl = mid, fm
                        else:
                            hi = mid
                    roots.append((lo+hi)/2)
    roots.sort()
    clean = [roots[0]]
    for r in roots[1:]:
        if r-clean[-1] > mp.mpf('1e-55'):
            clean.append(r)
    total = mp.mpf(0)
    for lo, hi in zip(clean, clean[1:]):
        winner = max(range(len(functions)),key=lambda i: functions[i]((lo+hi)/2))
        total += mp.quad(functions[winner], [lo,hi])
    return total


def exact_functions(beta,t):
    ra = -mp.expm1(mp.log1p(-t)/beta)/2
    rd = mp.expm1(-mp.log1p(-t)/(2-beta))/2
    return [lambda u:ra/t-u,
            lambda u:(-t*u-mp.expm1(mp.log1p(-t)+(1-beta)*mp.log1p(-2*t*u)))/t-u,
            lambda u:(t*u-mp.expm1(mp.log1p(-t)+(2-beta)*mp.log1p(2*t*u)))/t-u,
            lambda u:rd/t-u]


def polynomials(beta,t,n):
    coeff=[[] for _ in range(4)]
    for j in range(1,n+1):
        coeff[0].append((-1)**(j+1)*mp.binomial(1/beta,j)/2*t**(j-1))
        coeff[3].append((-1)**j*mp.binomial(-1/(2-beta),j)/2*t**(j-1))
        coeff[1].append((-mp.binomial(1-beta,j)*(-2)**j*t**(j-1),
                         mp.binomial(1-beta,j-1)*(-2)**(j-1)*t**(j-1)))
        coeff[2].append((-mp.binomial(2-beta,j)*2**j*t**(j-1),
                         mp.binomial(2-beta,j-1)*2**(j-1)*t**(j-1)))
    return [lambda u:sum(coeff[0])-u,
            lambda u:sum(a*u**j+b*u**(j-1) for j,(a,b) in enumerate(coeff[1],1))-2*u,
            lambda u:sum(a*u**j+b*u**(j-1) for j,(a,b) in enumerate(coeff[2],1)),
            lambda u:sum(coeff[3])-u]


def gfun(q):
    if q<=0:
        return mp.mpf(0)
    if q>=1:
        return 2*q-mp.mpf(2)/3
    return 2*q*q-2*q**3/3


def run(job):
    mp.mp.dps=80
    mode,pstr,tstr=job
    p,t=mp.mpf(pstr),mp.mpf(tstr)
    beta=p if mode=='uniform' else (mp.mpf('0.5') if mode=='low' else mp.mpf('1.5'))+p*t
    top=mp.mpf(5) if mode=='uniform' else mp.mpf(2)
    area=8*t*t*integrate_envelope(exact_functions(beta,t),top)
    if beta == 1:
        assert abs(area-2*t*t/(1-t)) < mp.mpf('1e-55')
    if mode!='uniform':
        f=-4-16*p+gfun(1+4*p) if mode=='low' else 12+16*p+gfun(-3-4*p)
        return {'mode':mode,'parameter':pstr,'t':tstr,'area':mp.nstr(area,55),
                'error_over_t4':float((area-4*t*t-f*t**3)/t**4)}
    rows=[]
    for n in range(1,7):
        approx=8*t*t*integrate_envelope(polynomials(beta,t,n),top)
        if n == 1:
            c = min(beta, 2-beta)
            tangent = 1/c**2 if c < mp.mpf('0.5') else 2/c
            assert abs(approx-tangent*t*t) < mp.mpf('1e-55')
        scaled=abs(area-approx)/t**(n+2)
        # J=[0.3,1.8], eta=0.2, K=B=5, exactly as in the manuscript.
        dk=lambda k:10**k*mp.rf(5,k)*2**(5+k)
        en=(dk(n+1)+(n+1)*dk(n))/mp.factorial(n+1)
        ratio=scaled/(40*en)
        assert ratio<1
        rows.append({'N':n,'scaled_error':float(scaled),'analytic_bound_ratio':float(ratio)})
    return {'mode':mode,'beta':pstr,'t':tstr,'checks':rows}


def main():
    started=time.perf_counter()
    load_before=os.getloadavg()
    betas=['0.3','0.49','0.5','0.51','0.9','1','1.1','1.49','1.5','1.51','1.8']
    jobs=[('uniform',b,t) for b in betas for t in ['0.02','0.003']]
    ps=['-1.25','-1','-0.875','-0.75','-0.5','-0.25','-0.125','0','0.5','1']
    jobs.append(('uniform','1.75','0.0001'))
    jobs += [(mode,p,t) for mode in ['low','high'] for p in ps for t in ['0.02','0.01','0.003','0.0001']]
    with concurrent.futures.ProcessPoolExecutor(max_workers=8) as pool:
        rows=list(pool.map(run,jobs))
    order_summary=[]
    for n in range(1,7):
        checks=[r['checks'][n-1] for r in rows if r['mode']=='uniform']
        order_summary.append({'N':n,'integral_cases':len(checks),
                 'max_error_over_t_Nplus2':max(c['scaled_error'] for c in checks),
                 'max_ratio_to_analytic_bound':max(c['analytic_bound_ratio'] for c in checks)})
    sharp_case=next(r for r in rows if r['mode']=='uniform' and r['beta']=='1.75')
    sharp_rows=[]
    for entry in sharp_case['checks']:
        n=entry['N']
        predicted=float(8*mp.rf(4,n+1)/mp.factorial(n+1))
        sharp_rows.append({'N':n,'predicted_limit':predicted,
                           'observed_scaled_error_t_0p0001':entry['scaled_error']})
    side=[r for r in rows if r['mode']!='uniform']
    result={'kind':'non-certified-high-precision-numerical-evidence','python':platform.python_version(),
            'mpmath':mp.__version__,'decimal_digits':80,'workers':8,'random_seed':None,
            'compute':{'logical_cpus':os.cpu_count(),'elapsed_seconds':time.perf_counter()-started,
                       'child_cpu_seconds':resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime,
                       'largest_child_peak_rss_bytes':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
                       'load_before':load_before,'load_after':os.getloadavg(),
                       'rss_units':'bytes on macOS; maximum single-child peak, not sum'},
            'root_scan':'64 uniform plus 12 geometric nodes, 220 bisection steps; not a root-count certificate',
            'uniform_beta_grid':betas,'uniform_t_grid':['0.02','0.003'],
            'uniform_J':['0.3','1.8'],'uniform_K':5,'uniform_orders':order_summary,
            'sharpness_probe_beta':'1.75','sharpness_probe':sharp_rows,
            'side_rows':side}
    path=Path(__file__).resolve().parents[1]/'results'/'uniform_square_python.json'
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='side_rows'},indent=2))
    print('Side cases:',len(side),'max scaled error:',max(abs(r['error_over_t4']) for r in side))
    matlab=path.with_name('uniform_square_matlab.json')
    if matlab.exists():
        errors=[]
        for mode,p,t,area,*_ in json.loads(matlab.read_text())['rows']:
            ref=next(r for r in side if r['mode']==('low' if mode==1 else 'high')
                     and float(r['parameter'])==p and float(r['t'])==t)
            errors.append(abs(area/float(ref['area'])-1))
        compare={'kind':'cross-tool-numerical-comparison','cases':len(errors),'max_relative_area_difference':max(errors)}
        path.with_name('uniform_square_crosscheck.json').write_text(json.dumps(compare,indent=2)+'\n')
        print(compare)

if __name__=='__main__':
    main()
