"""Actual one-dimensional source/query diagnostics; no n-length sample or Gaussian model."""
from pathlib import Path
import importlib.util,sys,json,math,time,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'round038_highdim/source_tables.py'
sys.dont_write_bytecode=True
spec=importlib.util.spec_from_file_location('source38',SOURCE)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
def stats(plus,minus,weights):
    a=plus+minus;p=plus/a;z=np.log(a)
    mass=float(np.dot(weights,a));q=weights*a/mass
    mz=float(q@z);mp=float(q@p)
    vz=float(q@(z-mz)**2);vp=float(q@(p-mp)**2)
    cov=float(q@((z-mz)*(p-mp)))
    return {'query_tilt_mass':mass,'mean_log_total':mz,'variance_log_total':vz,
            'mean_plus_posterior':mp,'variance_plus_posterior':vp,
            'covariance_log_total_posterior':cov,
            'correlation_log_total_posterior':cov/math.sqrt(vz*vp),
            'p_min':float(p.min()),'p_max':float(p.max())}
def point_values(src,x,L,b):
    plus=np.empty(len(x));minus=np.empty(len(x))
    for first in range(0,len(x),256):
        xx=x[first:first+256]
        plus[first:first+256]=src.point_class_kernel(xx,L,b,True)
        minus[first:first+256]=src.point_class_kernel(xx,L,b,False)
    return plus,minus
def main():
    target=HERE/'correlation_results.json'
    if target.exists():raise RuntimeError('Refusing replacement')
    started=time.time();src=mod.Source();L=7/6
    reg={'source_SHA256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
         'L':'7/6','tables':[0,1,2,3],'methods':['bin_mean_channel','point_midpoint','point_Gauss2_per_bin'],
         'n':src.n,'Q':src.Q,'scope':'True finite source formulas; deterministic 1D quadrature; no MC/CLT claim; ordinary FP'}
    (HERE/'correlation_registration.json').write_text(json.dumps(reg,indent=2)+'\n')
    rows=[]
    for b in reg['tables']:
        evaluated=src.evaluate(L,b)
        methods={'bin_mean_channel':stats(evaluated['mean_plus'],evaluated['mean_minus'],np.full(src.Q,1/src.Q))}
        plus,minus=point_values(src,src.bins_center,L,b)
        methods['point_midpoint']=stats(plus,minus,np.full(src.Q,1/src.Q))
        x=(src.bins_center[:,None]+src.binwidth*np.array([-1,1])[None,:]/(2*math.sqrt(3))).ravel()
        plus,minus=point_values(src,x,L,b)
        methods['point_Gauss2_per_bin']=stats(plus,minus,np.full(len(x),1/len(x)))
        for method,d in methods.items():
            d.update({'table':b,'method':method,'source_eta':src.pchi,
                      'E_p_minus_source_eta':d['mean_plus_posterior']-src.pchi,
                      'n_times_mean_log_total':src.n*d['mean_log_total'],
                      'n_times_variance_log_total':src.n*d['variance_log_total'],
                      'n_times_variance_plus_posterior':src.n*d['variance_plus_posterior']})
            rows.append(d)
        print('table',b,'corr',methods['point_Gauss2_per_bin']['correlation_log_total_posterior'],flush=True)
    out={'registration':reg,'rows':rows,'seconds':time.time()-started,
         'qualification':'All quadrature/FP un-certified. Bin-channel and continuous point quadratures explicitly distinct.'}
    target.write_text(json.dumps(out,indent=2)+'\n')
    print('seconds',out['seconds'],flush=True)
if __name__=='__main__':main()
