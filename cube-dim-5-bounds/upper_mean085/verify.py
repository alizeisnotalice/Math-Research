#!/usr/bin/env python3
"""Independent exact moment check by the uniquely largest active coordinate."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import os
import time
import sympy as sp

HERE=Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rat(x):
    x=F(x)
    return sp.Rational(x.numerator,x.denominator)


def conditional_moment(bounds,anchor):
    intervals=[(F(lo)-F(c),F(hi)-F(c)) for (lo,hi),c in zip(bounds,anchor)]
    active=[j for j,(lo,hi) in enumerate(intervals) if hi>lo]
    fixed=[abs(lo) for lo,hi in intervals if lo==hi]
    threshold=max(fixed,default=F(0))
    # Atom: all active absolute coordinates are <= the largest fixed one.
    atom_probability=F(1)
    for j in active:
        lo,hi=intervals[j]
        atom_probability*=max(F(0),min(hi,threshold)-max(lo,-threshold))/(hi-lo)
    answer=rat(threshold**2*atom_probability)
    x=sp.Symbol('x',real=True)
    knots={F(0),threshold}|{abs(v) for span in intervals for v in span}
    for i in active:
        ilo,ihi=intervals[i]
        cuts=sorted({ilo,ihi}|{sign*v for v in knots for sign in (-1,1) if ilo<sign*v<ihi})
        for left,right in zip(cuts,cuts[1:]):
            mid=(left+right)/2
            radius=abs(mid)
            if radius<=threshold:
                continue
            r=x if mid>0 else -x
            integrand=x*x/rat(ihi-ilo)
            for j in active:
                if i==j:
                    continue
                lo,hi=intervals[j]
                overlap=min(hi,radius)-max(lo,-radius)
                if overlap<=0:
                    integrand=sp.Integer(0)
                    break
                top=rat(hi) if hi<radius else r
                bottom=rat(lo) if lo>-radius else -r
                integrand*= (top-bottom)/rat(hi-lo)
            answer+=sp.integrate(sp.expand(integrand),(x,rat(left),rat(right)))
    answer=sp.cancel(answer)
    assert answer.is_Rational
    return F(int(answer.p),int(answer.q))


def main():
    start=time.monotonic()
    reg=json.loads((HERE/'registration.json').read_text())
    result=json.loads((HERE/'results.json').read_text())
    assert result['registration_sha256']==sha(HERE/'registration.json')
    assert result['source_sha256']==sha(HERE/'certify.py')
    assert result['moment_source_sha256']==sha(HERE/'moment.py')
    for row in reg['sources']:
        assert sha(row['path'])==row['sha256']
    rows=[]
    for case,saved in zip(reg['cases'],result['cases']):
        assert case['id']==saved['id']
        m2=conditional_moment(case['bounds'],case['anchor_tail'])
        assert m2==F(saved['E_max_abs_squared'])
        mu=sp.Matrix([rat((F(lo)+F(hi))/2-F(c)) for (lo,hi),c in zip(case['bounds'],case['anchor_tail'])])
        covariance=sp.diag(*[rat((F(hi)-F(lo))**2/12) for lo,hi in case['bounds']])
        Q=sp.eye(4)/4-sp.ones(4)/36
        eta=sp.cancel(sp.trace(Q*covariance)+(mu.T*Q*mu)[0])
        eta=F(int(eta.p),int(eta.q))
        assert eta==F(saved['eta_old_Q'])
        assert F(25,48)*m2<=F(5,2)*eta
        assert F(25,48)*m2==F(saved['face_25_over_24_coefficient'])
        rows.append({'id':case['id'],'M2_independent':str(m2),
                     'eta_independent':str(eta),'matches':True})
    assert len(rows)==reg['case_count']==16
    out={'id':'upper_mean085_independent_exact_verification','status':'all_16_conditional_max_moments_match',
         'pid':os.getpid(),'source_sha256':sha(__file__),'input_results_sha256':sha(HERE/'results.json'),
         'sympy_version':sp.__version__,'method':'Condition on unique largest active coordinate, add fixed-coordinate atom; exact symbolic polynomial integration. No production moment import.',
         'cases':rows,'seconds':time.monotonic()-start,'new_H_gradient_Green_calls':0,
         'new_Arb_transcendental_calls':0,'no_wall_time_cap':True}
    (HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['status','pid','seconds']}))


if __name__=='__main__':main()
