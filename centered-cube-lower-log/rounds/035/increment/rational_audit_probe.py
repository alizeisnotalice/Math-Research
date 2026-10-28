"""Exact rational autocorrelation/increment/one-sided derivative checks."""
from fractions import Fraction as F
from bisect import bisect_left,bisect_right
from itertools import product
from pathlib import Path
import json

def atoms(u):
    k=len(u)
    return [(j+x)/k for j,x in enumerate(u)]
def primitive(a,x,left=False):
    x=x%1
    return F((bisect_left if left else bisect_right)(a,x))-len(a)*x
def corr(a,t):
    k=len(a); t=t%1
    cuts=sorted({F(0),F(1),*a,*[(x-t)%1 for x in a]})
    value=F(0)
    for l,r in zip(cuts,cuts[1:]):
        mid=(l+r)/2;ell=r-l
        value+=ell*primitive(a,mid)*primitive(a,mid+t)+F(k*k,12)*ell**3
    return value
def deriv_right(a,t):
    mean=F(1,2)-sum((len(a)*x-j for j,x in enumerate(a)),F(0))/len(a)
    return sum((primitive(a,x-t,left=True)for x in a),F(0))-len(a)*mean

configs=[[F(1,10),F(9,10),F(1,10)],
         [F(99,100),F(1,100),F(3,7),F(5,9),F(1,3)],
         [F(1,2)]*5]
local_checks=0; derivative_checks=0; global_checks=0; max_change=F(0)
for u in configs:
    a=atoms(u);k=len(a)
    gaps=[(a[(j+1)%k]-a[j])%1 for j in range(k)]
    for tau in [F(1,100),F(1,10),F(1,2),F(1)]:
        t=tau/k
        candidate=tau/2-tau*tau/2+sum((max(F(0),t-d)for d in gaps),F(0))
        assert corr(a,0)-corr(a,t)==candidate
        local_checks+=1
    differences={(x-y)%1 for x in a for y in a}
    for t in [F(0),F(1,17),F(2,7),F(5,8),F(16,17)]:
        ordered=[(xx-yy)%1 for xx in a for yy in a if xx!=yy]
        global_candidate=F(k)*t/2-F(k*k)*t*t/2+sum((max(F(0),t-z)for z in ordered),F(0))
        assert corr(a,0)-corr(a,t)==global_candidate
        global_checks+=1
        distance=min(((z-t)%1 for z in differences if (z-t)%1>0),default=F(1))
        h=min(distance/3,F(1,100000))
        exact_slope=(corr(a,t+h)-corr(a,t))/h-F(k*k)*h/2
        assert exact_slope==deriv_right(a,t)
        derivative_checks+=1

# Bounded-difference enumeration, including deterministic jump coincidences.
us=[F(1,6),F(1,2),F(5,6)]
for u in product(us,repeat=3):
    a=atoms(u)
    for l in range(3):
        for new in us:
            v=list(u);v[l]=new;b=atoms(v)
            for t in [F(0),F(1,13),F(1,3),F(1,2),F(11,13)]:
                change=abs(deriv_right(a,t)-deriv_right(b,t))
                max_change=max(max_change,change)
                assert change<=5

# Exact integration of the single gap-overlap moments, for tau<=1/2.
gap_moments=[]
for tau in [F(1,100),F(1,10),F(1,2)]:
    mu=tau**3/6
    second=tau**4/12
    var_sum_coefficient=second-3*mu*mu
    assert var_sum_coefficient==tau**4*(1-tau*tau)/12
    gap_moments.append(dict(tau=str(tau),mean=str(mu),second=str(second),
                            variance_S_times_K=str(var_sum_coefficient)))
out=dict(scope='Exact rational small configurations only; no MC.',
         local_increment_checks=local_checks,global_ordered_pair_checks=global_checks,
         one_sided_derivative_checks=derivative_checks,
         max_derivative_single_jitter_change_enumerated=str(max_change),
         gap_moments=gap_moments,
         endpoint_note='Derivative is right-sided at atom differences; Rprime(0+)=-K/2, while 1- has K/2.')
# Exact polynomial integrals: R(0) is an average of U^2-U+1/3.
poly=[F(1,3),F(-1),F(1)]
mean=sum((a/F(i+1)for i,a in enumerate(poly)),F(0))
square=[sum((poly[j]*poly[i-j]for j in range(len(poly)) if 0<=i-j<len(poly)),F(0))for i in range(5)]
second=sum((a/F(i+1)for i,a in enumerate(square)),F(0))
assert mean==F(1,6) and second-mean*mean==F(1,180)
out['R0_mean']=str(mean)
out['R0_variance_times_K']=str(second-mean*mean)
Path(__file__).with_name('rational_audit_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
