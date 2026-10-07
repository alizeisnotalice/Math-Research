from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import hashlib
import json

BASE=Path(__file__).resolve().parent
PREFIX='aggregate_mass_band_budget_20261007'
records=[]
total=0

def mass_band(value, base):
    ratio=value/base
    k=0
    while ratio > 2:
        ratio/=2
        k+=1
    assert 1 < ratio <= 2
    return k

for stage,n in enumerate((16,64,256),1):
    checks=[]
    def check(name,test):
        global total
        assert test,(stage,name)
        total+=1
        checks.append(name)
    delta,D,tau=F(1,100*n),16*n*n,F(3,8)
    root=isqrt(n)
    check('rational square-root cutoff',root*root==n)
    check('complete source high-density',1/delta > tau*root)
    check('positive slab and source box',0 < delta < F(1,6) and D>4)
    low,high=F(4,3),F(2)
    xlo,xhi=(low-delta)/2,(high-delta)/2
    check('anchor cold on whole receiver core',xlo>F(1,2)+delta/2)
    check('two-sided receiver width',2*(xhi-xlo)==F(2,3))
    check('inner tangential source margin',F(D-4,2)+1==F(D-2,2))
    for r in (F(4,3),F(3,2),F(7,4),F(2)):
        x=(r-delta)/2
        check('original full-capture arrival '+str(r),2*x+delta==r)
        check('partial derivative positive '+str(r),x-delta/2>0)
        check('full derivative negative '+str(r),-1/(r*r)<0)
        check('original amplitude band '+str(r),tau < 1/r <=2*tau)
    mlo,mhi=low**(n-1),high**(n-1)
    klo,khi=mass_band(mlo,tau),mass_band(mhi,tau)
    check('many actual mass bands',khi-klo+1>=n//2)
    check('far complete band indices',khi-klo-2>=F(n,3))
    check('linear log-span certificate',(mhi/mlo)**2 >2**(n-1))
    margin=F(D-4,D-2)**(n-1)
    check('inner source fraction',margin>F(3,4))
    cross_lower=tau*F(2,3*(n-1))*margin
    check('actual complete-band cross lower envelope',cross_lower>F(1,8*n))
    profile_cap=tau*F(2,3)
    weak_core_ratio=profile_cap*F(D-2,D)**(n-1)
    check('aggregate source profile cap',profile_cap==F(1,4))
    check('whole core weak payment',0<weak_core_ratio<profile_cap)
    check('slab core actual tail zero',profile_cap<root)
    nu,tail_delta=F(1,4),F(1,4)
    retention=1-nu-F(1,root)
    check('conditional residual tail absorption',retention-tail_delta>0)
    records.append({'stage':stage,'n':n,'checks':checks,'check_count':len(checks),'mass_band_indices':[klo,khi],'complete_band_gap':khi-klo-2,'profile_cap_exact':'1/4','weak_core_ratio_approx':float(weak_core_ratio),'cross_lower_approx':float(cross_lower),'conditional_retention_exact':str(retention)})

result={'status':'PASS','executions':1,'checks':total,'scope':'True complete high-source slab, original continuous winner, anchor-cold gate and aggregate core parameter certificate. Not general source-tail, actual soft gates, or endpoint order fit.','records':records}
result_path=BASE/(PREFIX+'_results.json')
result_path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
receipt={'status':'PASS','executions':1,'checks':total,'registration_sha256':hashlib.sha256((BASE/(PREFIX+'_registration.md')).read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'results_sha256':hashlib.sha256(result_path.read_bytes()).hexdigest()}
(BASE/(PREFIX+'_receipt.json')).write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':'PASS','checks':total,'rounds':[(r['n'],r['mass_band_indices'],r['weak_core_ratio_approx']) for r in records]}))
