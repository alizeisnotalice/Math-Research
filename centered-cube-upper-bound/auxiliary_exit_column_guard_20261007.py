"""Exact arithmetic diagnostics for a gate-free auxiliary exit column.
This does not construct actual FIRST inputs. No randomness or fitted order.
"""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json

def log2_interval(m):
    lo = sum((F(2, (2*k+1)*3**(2*k+1)) for k in range(m)), F(0))
    return lo, lo + F(2, (2*m+1)*3**(2*m+1))*F(9,8)

lo, hi = log2_interval(24)
lo2, hi2 = log2_interval(40)
assert lo < lo2 < hi2 < hi
rows=[]
for stage, dims in enumerate(((576,1024,2304),(4096,6400,9216),(16384,25600,36864)),1):
    for n in dims:
        # sigma=2/n; survival <=(1-1/(2n))^n, all exact.
        den=(2*n)**n
        num=den-(2*n-1)**n
        assert 3*num > den
        rt=isqrt(n)
        assert rt*rt==n
        lower=(1+n*lo)/3
        upper=1+n*hi
        # Independent cone integration: core 1 plus n*int_a^b dR/R.
        split=F(1,3)+F(n,3)*lo
        assert lower==split
        rows.append(dict(stage=stage,n=n,sigma=str(F(2,n)),
                         exact_exit_lower_gt_one_third=True,
                         column_lower=float(lower),column_upper=float(upper),
                         column_lower_over_sqrt_n=float(lower/rt),
                         log2_interval_width=float(hi-lo)))
out=dict(scope='gate-free kernel only; not actual FIRST or geom closure',
         random_seed=None,three_stages=rows,
         log2_interval=[str(lo),str(hi)],all_checks_passed=True)
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(rows=len(rows),all_checks_passed=True,
                     last_lower_over_sqrt_n=rows[-1]['column_lower_over_sqrt_n'])))
