"""Fixed 4x8 generic points plus tau0 four calibrations; one mp80 group."""
from pathlib import Path
from fractions import Fraction as F
import runpy,json,time,math,hashlib,warnings
import numpy as np
import mpmath as mp
from scipy.integrate import quad,IntegrationWarning
OUT=Path(__file__).resolve().parent
NS=runpy.run_path(str(OUT/'cb11_frame_green.py'))
POINTS=[(1,F(2,7),F(1,11)),(1,F(3,11),F(2,13)),(1,F(4,9),F(1,17)),(1,F(5,11),F(3,19)),(1,F(2,3),F(1,13)),(1,F(7,10),F(2,23)),(1,F(9,11),F(4,25)),(1,F(11,13),F(5,29))]
TAUS=[F(1,10),F(1,3),F(1),F(3)]

def run():
    start=time.monotonic();rows=[];allcases=[(tau,x) for tau in TAUS for x in POINTS]+[(F(0),x) for x in POINTS[:4]]
    # Algebraic nondegeneracy before evaluating the preregistered fixed cases.
    for tau,x in allcases:
        for cone,rates,_ in NS['_prepared'](tau):
            b=[sum(F(v)*xx for v,xx in zip(r,x))/k for r,k in zip(cone,rates)]
            assert 0 not in b and len(set(b))==3,(tau,x,cone,b)
    oldpath=Path('/Users/zhengzhihao/Desktop/T/Math-Research-finite-dim/cube-finite-dim-bounds/upper/continuous_bridge/cb9_frame_green.py')
    assert hashlib.sha256(oldpath.read_bytes()).hexdigest()=='d23501eeaf417a736730a6e61398e27b9364c231d35ec6ee8bf573e357058bf5'
    old=runpy.run_path(str(oldpath))
    with mp.workdps(80):
        for tau,x in allcases:
            value=NS['green'](x,tau,method='float');high=NS['green'](x,tau,method='mp',dps=80);rel=float(abs(mp.mpf(value)-high)/high)
            assert rel<1e-8,(tau,x,value,high,rel)
            t=3/7;direct=0j
            for cone in NS['CONES']:
                rates=[float(F(a)+tau*b) for a,b in map(NS['rate_coeff'],cone)];ell=[sum(float(v)*float(xx) for v,xx in zip(r,x)) for r in cone]
                direct+=1/math.prod(complex(t*k,-b) for k,b in zip(rates,ell))
            direct/=8*math.pi**3;pair=NS['heat_kernel'](t,x,tau);prel=abs(direct.real-pair)/pair;assert prel<1e-10 and abs(direct.imag)<1e-12
            with warnings.catch_warnings(record=True) as ww:
                warnings.simplefilter('always',IntegrationWarning)
                a,e1=quad(lambda tt:NS['heat_kernel'](tt,x,tau),0,1,epsabs=2e-11,epsrel=2e-10,limit=250)
                b,e2=quad(lambda tt:NS['heat_kernel'](tt,x,tau),1,np.inf,epsabs=2e-11,epsrel=2e-10,limit=250)
            qrel=abs((a+b)-float(high))/float(high);assert qrel<2e-8,(tau,x,a+b,high,qrel)
            rec={'tau_exact':str(tau),'x_exact':list(map(str,x)),'G_float':value,'G_mp80':mp.nstr(high,82),'float_relative_difference':rel,'G_direct_time_integral':a+b,'time_quadrature_error_estimate':e1+e2,'time_relative_difference':qrel,'time_warnings':[str(w.message) for w in ww],'p_full96_complex':direct.real,'p_antipodal48':pair,'p_relative_difference':prel}
            if tau==0:
                reference=old['green'](x);rec['tau0_cb9_positive_integral']=reference;rec['tau0_relative_difference']=abs(reference-float(high))/float(high);assert rec['tau0_relative_difference']<1e-10
            rows.append(rec)
    res={'status':'PASS fixed 32 positive-tau cases + 4 tau0 calibrations','positive_taus':list(map(str,TAUS)),'points_per_positive_tau':8,'tau0_points':4,'only_high_precision_group_dps':80,'seconds':time.monotonic()-start,'records':rows,'module_sha256':hashlib.sha256((OUT/'cb11_frame_green.py').read_bytes()).hexdigest(),'scope':'representation/normalization calibration only; no cost search, no certified cube g/J, no C_n improvement'}
    (OUT/'cb11_results.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps({k:v for k,v in res.items() if k!='records'},indent=2));print('max float/mp rel',max(r['float_relative_difference'] for r in rows));print('max time rel',max(r['time_relative_difference'] for r in rows));print('warnings',sum(bool(r['time_warnings']) for r in rows))
if __name__=='__main__':run()
