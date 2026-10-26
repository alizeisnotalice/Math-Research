#!/usr/bin/env python3
"""Independent integer-only rescaled-root check of saved new ratios."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import os
import time

HERE=Path(__file__).resolve().parent
D=2**224


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def root(q):
    n=q.numerator*D*D//q.denominator
    lo,hi=0,1<<((n.bit_length()+1)//2)
    while hi-lo>1:
        mid=(lo+hi)//2
        if mid*mid<=n:lo=mid
        else:hi=mid
    low,high=F(lo,D),F(lo+1,D)
    assert low*low<=q<high*high
    return low,high


def main():
    start=time.monotonic();reg=json.loads((HERE/'registration.json').read_text())
    result=json.loads((HERE/'results.json').read_text())
    for key,name in [('registration_sha256','registration.json'),('source_sha256','certify.py'),('ratio_source_sha256','ratio.py')]:
        assert result[key]==sha(HERE/name)
    for src in reg['sources']:assert sha(src['path'])==src['sha256']
    assert len(result['lifts'])==len(reg['lift_cases'])==12
    for raw,lift in zip(reg['lift_cases'],result['lifts']):
        u=list(map(F,raw));radius=max(map(abs,u));a=list(map(F,lift['coefficients']));V=lift['columns']
        assert len(a)==len(V)==5 and all(v[0]==1 and all(abs(x)==1 for x in v) for v in V)
        assert [sum((a[i]*V[i][j] for i in range(5)),F(0)) for j in range(5)]==[F(0)]+u
        assert sum(a)==0 and sum(map(abs,a))==radius and a[0]==-radius/2
        assert all(v>=0 for v in a[1:])
        sorted_abs=sorted(map(abs,u),reverse=True)+[F(0)]
        assert a[1:]==[(sorted_abs[k]-sorted_abs[k+1])/2 for k in range(4)]
        assert lift['order']==sorted(range(4),key=lambda j:(-abs(u[j]),j))
        assert lift['signs']==[1 if v>=0 else -1 for v in u]
    checks=[]
    for row,old in zip(result['ratios'],reg['saved_matched_comparisons']):
        x=F(row['x']);lo1,hi1=root(16+x*x);lo2,hi2=root(256+x*x)
        low=(lo1+x)**2*(lo2+x)**8/16**9
        high=(hi1+x)**2*(hi2+x)**8/16**9
        arblo,arbhi=F(row['R_new']['lower']),F(row['R_new']['upper'])
        apiL,apiU=map(F,row['integer_api_interval'])
        assert apiL<=low<=arblo<arbhi<=high<=apiU
        assert high<F(old['old_asinh_lower'])
        assert F(row['log_R_new']['upper'])<x
        if old['old_exp_lower'] is not None:
            assert high<F(old['old_exp_lower'])
        assert F(old['source_cell_h'])/F(old['source_time_start'])==x
        checks.append({'x':str(x),'independent_lower':str(low),'independent_upper':str(high),
                       'encloses_saved_Arb_ball':True,'contained_in_integer_API_interval':True,
                       'strictly_below_saved_old_asinh':True,
                       'strictly_below_saved_exp':old['old_exp_lower'] is not None,
                       'saved_log_upper_below_x':True})
    assert len(checks)==6
    assert F(result['ratios'][4]['integer_api_interval'][1])<843
    out={'id':'upper_shift091_independent_integer_verification','status':'12_lifts_and_6_ratios_verified',
         'pid':os.getpid(),'source_sha256':sha(__file__),'input_results_sha256':sha(HERE/'results.json'),
         'rescaled_formula':'(sqrt(16+x^2)+x)^2*(sqrt(256+x^2)+x)^8/16^9',
         'root_denominator':str(D),'checks':checks,'seconds':time.monotonic()-start,
         'new_Arb_log_exp_asinh_H_gradient_Green_calls':0,'old_baseline_recomputations':0,
         'wall_time_cap':None}
    (HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['status','pid','seconds']}))


if __name__=='__main__':main()
