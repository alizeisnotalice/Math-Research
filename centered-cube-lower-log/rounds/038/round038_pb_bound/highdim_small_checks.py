"""Read-only high-dimensional implementation checks; no full cache, no RNG/MC."""
from pathlib import Path
from fractions import Fraction as F
import importlib.util,sys,json,math,hashlib
import numpy as np

HERE=Path(__file__).resolve().parent
OTHER=HERE.parent/'round038_highdim'
sys.dont_write_bytecode=True
sys.path.insert(0,str(OTHER))
spec=importlib.util.spec_from_file_location('audit_source',OTHER/'source_tables.py')
source=importlib.util.module_from_spec(spec);spec.loader.exec_module(source)
spec=importlib.util.spec_from_file_location('audit_driver',OTHER/'probe_highdim.py')
driver=importlib.util.module_from_spec(spec);spec.loader.exec_module(driver)
src=source.Source()
delta,h=F(1,src.w*src.K*src.v),F(1,src.Q)
support=(delta+h)/2
def cdf(u):
    a,b=max(delta,h),min(delta,h);d,e=(a+b)/2,(a-b)/2
    if u<=-d:return F(0)
    if u<-e:return (u+d)**2/(2*a*b)
    if u<=e:return F(1,2)+u/a
    if u<d:return 1-(d-u)**2/(2*a*b)
    return F(1)
cases=0;max_count_error=0.
for table in range(4):
    for plus in [True,False]:
        indices=range(src.cut)if plus else range(src.cut,src.K)
        base=[]
        for j in indices:
            label=(int(src.labels[j])+(table if j<src.cut else -table))%4
            base.append((F(j)+F(2*label+1,8))/src.K)
        atoms=[a+s for s in [-1,0,1,2] for a in base]
        for L in [F(10,9),F(7,6),F(11,9)]:
            ell=L-1
            centers=[]
            for atom in [base[0],base[-1]]:
                for side in [-1,1]:
                    for shift in [-support,0,support]:
                        center=atom+shift+side*ell/2
                        center=(center-F(1,2))%1+F(1,2)
                        centers.append(center)
            expected=[sum((cdf(c+ell/2-a)-cdf(c-ell/2-a)for a in atoms),F(0))for c in centers]
            actual,transition=src._peak_bin_mean_and_transition_count(
                np.array(list(map(float,centers))),float(ell),table,plus)
            err=max(abs(float(e)-v)for e,v in zip(expected,actual))
            assert err<1e-9 and transition.max()<=2
            max_count_error=max(max_count_error,err);cases+=len(centers)

ratio_cases=0
for a,b in [(F(3,10),F(7,10)),(F(1,2),F(1,2)),(F(2,5),F(1,2))]:
    for e in [F(1,100),F(1,20)]:
        p=a/(a+b)
        for da in [-e,0,e]:
            for db in [-e,0,e]:
                true=(a+da)/(a+b+da+db)
                assert abs(true-p)<=e/(a+b)
                assert abs(true*(1-true)-p*(1-p))<=e/(a+b)
                ratio_cases+=1
assert F(2803,5000)**2*F(29,50)**2/F(21,100)<1
assert F(49,151)>F(3,10)
corners=0
for cut in [0,3,9]:
    for mu in [1.,4.,10.]:
        for var in [1.,3.]:
            e=.02
            low,up=driver.cdf_enclosure(cut,np.array([mu]),np.array([var]),np.array([e]),0.)
            for mm in [mu-e,mu,mu+e]:
                for vv in [var-e,var,var+e]:
                    value=.5*math.erfc(-(cut+.5-mm)/math.sqrt(2*vv))
                    assert low[0]<=value+1e-15 and value<=up[0]+1e-15
                    corners+=1

reg=json.loads((OTHER/'HIGH_DIM_REGISTRATION.json').read_text())
plan=json.loads((OTHER/'HIGH_DIM_PLAN.json').read_text())
result=json.loads((OTHER/'HIGH_DIM_RESULTS.json').read_text())
for key,name in [('script_SHA256','probe_highdim.py'),('source_SHA256','source_tables.py'),('builder_SHA256','build_cache.py')]:
    assert plan[key]==hashlib.sha256((OTHER/name).read_bytes()).hexdigest()
assert result['registration_SHA256']==hashlib.sha256((OTHER/'HIGH_DIM_REGISTRATION.json').read_bytes()).hexdigest()
assert len(set(reg['cuts']))==3 and 0<=reg['cuts'][0]<reg['cuts'][-1]<src.n
out={'scope':'No cache build, no RNG, no new high-dimensional samples',
     'peak_endpoint_Fraction_cases':cases,'max_peak_count_FP_error':max_count_error,
     'normalized_class_ratio_Fraction_cases':ratio_cases,
     'CDF_rectangle_corner_interior_cases':corners,'BE_constant_rational_check':True,
     'p_range_from_eta49to51_rational_check':True,'frozen_hashes_match':True,
     'cuts_distinct_inside_support':reg['cuts'],'N':reg['N'],'J':reg['J'],
     'known_FP_certification_limit':'Floating count errors are observed only; no interval certification.',
     'inputs_SHA256':{name:hashlib.sha256((OTHER/name).read_bytes()).hexdigest()for name in [
       'source_tables.py','build_cache.py','probe_highdim.py','HIGH_DIM_REGISTRATION.json','HIGH_DIM_RESULTS.json']}}
target=HERE/'highdim_small_checks.json'
if target.exists():raise RuntimeError('Refusing replacement')
target.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
