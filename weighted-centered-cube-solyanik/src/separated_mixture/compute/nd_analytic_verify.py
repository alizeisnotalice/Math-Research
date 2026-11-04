#!/usr/bin/env python3
"""Independent closed-form audit of saved nD Gauss-24 residual masses.

The saved experiment has its other site at (D,0,...), D=1.3 and r<=1/8.
Thus ||x-(D,0,...)||_infty=D-x_1 on every measured cube.  No quadrature is
needed for that residual.  This script deliberately reads the existing summary
and does not rerun the 12/24/48 Gauss node evaluations.
"""
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parent

def power_difference(upper, lower, exponent):
    # upper^p-lower^p, stable also for the negative exponents used here.
    return lower**exponent * math.expm1(exponent*math.log1p((upper-lower)/lower))

def other_mass(n, c, D, sigma, radius):
    exponent=sigma-n+1
    return c*(2*radius)**(n-1)*power_difference(D+radius,D-radius,exponent)/exponent

def main():
    source=json.loads((ROOT/'nd_blowup_summary.json').read_text())
    checks=[]
    for item in source['results']:
        n=item['n']; tau=item['tau']
        for row in item['rows']:
            for name,box in row['boxes'].items():
                r=box['radius']; assert 0<r<1.3/2
                pure=.8*(2**n)*n*r**tau/tau
                assert math.isclose(pure,box['pure_exact'],rel_tol=2e-14,abs_tol=1e-14)
                analytic=pure + .35*(2*r)**n + other_mass(n,.55,1.3,.75,r)
                error=abs(analytic-box['mixed_gauss24'])
                assert math.isfinite(error) and error<=1e-12*max(1.,abs(analytic))
                checks.append({'n':n,'tau':tau,'k':row['k'],'box':name,'radius':r,
                               'analytic_mixed':analytic,'saved_gauss24':box['mixed_gauss24'],
                               'absolute_error':error})
    assert len(checks)==140
    out={'scope':'closed-form audit of Gauss24 nonsingular residual only',
         'formula':'c(2r)^(n-1)*((D+r)^(sigma-n+1)-(D-r)^(sigma-n+1))/(sigma-n+1)',
         'checks':checks,'count':len(checks),'analytic_max_abs_error':max(q['absolute_error'] for q in checks),
         'note':'No Gauss nodes rerun; this audits the existing saved values.'}
    (ROOT/'nd_analytic_audit.json').write_text(json.dumps(out,indent=2,allow_nan=False))
    print(f"checks={out['count']} analytic_max_abs_error={out['analytic_max_abs_error']:.3e}")

if __name__=='__main__': main()
