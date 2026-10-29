"""Exact tiny mixed-tilt checks and deterministic 1D finite-source scans, no tail MC."""
from pathlib import Path
from fractions import Fraction as F
import sys,importlib.util,json,math,itertools,time,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent
sys.dont_write_bytecode=True
def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
old=load(HERE.parent/'round038_highdim/source_tables.py','oldsource')
new=load(HERE.parent/'round040_source/source_tables.py','newsource')
def stats(plus,minus):
    a=plus+minus;p=plus/a;z=np.log(a);q=a/a.sum()
    eta=float(q@p);mu=float(q@z);vz=float(q@(z-mu)**2);vp=float(q@(p-eta)**2)
    cov=float(q@((z-mu)*(p-eta)))
    return {'eta':eta,'mean_log_total':mu,'variance_log_total':vz,'variance_p':vp,
            'covariance':cov,'correlation':cov/math.sqrt(vz*vp),
            'signal_fraction':vp/(eta*(1-eta))}
def pb(ps):
    out=[F(1)]
    for p in ps:
        nxt=[F(0)]*(len(out)+1)
        for k,v in enumerate(out):nxt[k]+=v*(1-p);nxt[k+1]+=v*p
        out=nxt
    return out
def main():
    target=HERE/'mixed_tilt_results.json'
    if target.exists():raise RuntimeError('Refusing replacement')
    started=time.time();src=old.Source();L=7/6;response=[]
    x=(src.bins_center[:,None]+src.binwidth*np.array([-1,1])[None,:]/(2*math.sqrt(3))).ravel()
    for b in range(4):
        zz=[];pp=[];aa=[]
        for first in range(0,len(x),256):
            plus=src.point_class_kernel(x[first:first+256],L,b,True)
            minus=src.point_class_kernel(x[first:first+256],L,b,False)
            a=plus+minus;aa.extend(a);pp.extend(plus/a);zz.extend(np.log(a))
        a=np.array(aa);p=np.array(pp);z=np.array(zz);eta=src.pchi
        for t in [-2.,-1.,-.5,0.,.5,1.,2.,4.]:
            weight=a*np.exp(t*z);norm=float(weight.mean()/a.mean());et=float(weight@p/weight.sum())
            cov=float((weight/weight.sum())@((z-float(weight@z/weight.sum()))*(p-et)))
            response.append({'b':b,'t':t,'eta_t':et,'eta_t_minus_eta':et-eta,
                             'standardized_count_shift':math.sqrt(src.n)*(et-eta)/math.sqrt(eta*(1-eta)),
                             'per_coordinate_score_MGF':norm,'eta_t_derivative_covariance':cov})
    scan=[]
    for R in [4,8,16,32]:
        source=new.Source(R=R)
        if R==4:assert np.array_equal(source.atoms,src.atoms)
        for length in [10/9,7/6,11/9]:
            for b in range(R):
                d=source.evaluate(length,b)
                row=stats(d['mean_plus'],d['mean_minus'])
                row.update({'R':R,'b':b,'L':length,'method':'bin-channel only',
                            'n_variance_log_total':source.n*row['variance_log_total']})
                scan.append(row)
        print('new-source R',R,'done',flush=True)
    # Exact finite identity E[w(X_1)...w(X_n) PB(I|X)] = m^n Bin(n,eta_w)(I).
    ps=[F(1,3),F(2,3)];scores=[F(3,4),F(5,4)];n=3;exact=[]
    for t in [-1,0,1,2]:
        m=sum(a**t for a in scores)/2;et=sum(a**t*p for a,p in zip(scores,ps))/(2*m)
        lhs=F(0)
        for obs in itertools.product([0,1],repeat=n):
            post=pb([ps[o]for o in obs]);w=math.prod(scores[o]**t for o in obs)
            lhs+=w*sum(post[:2],F(0))/2**n
        rhs=m**n*sum(pb([et]*n)[:2],F(0));assert lhs==rhs
        exact.append({'t':t,'eta_t':str(et),'LHS':str(lhs),'RHS':str(rhs),'equal':True})
    # Zero covariance does not imply mixed-tilt independence.
    probs=[F(9,20),F(3,5),F(9,20)];scores=[F(4,5),F(1),F(5,4)]
    eta=sum(probs)/3;assert eta==F(1,2)
    # log scores are (-h,0,h), symmetric endpoints p are equal: covariance exactly zero.
    eta1=sum(a*p for a,p in zip(scores,probs))/sum(scores)
    assert eta1!=eta
    out={'continuous_response':response,'new_source_bin_scan':scan,'Fraction_mixed_identity':exact,
         'zero_covariance_counterexample':{'p':list(map(str,probs)),'positive_score_weights':list(map(str,scores)),
                                          'eta':str(eta),'eta_at_t1':str(eta1),'log_score_covariance_exact':0},
         'source_SHA256':{name:hashlib.sha256((HERE.parent/path).read_bytes()).hexdigest()for name,path in [
             ('38','round038_highdim/source_tables.py'),('40','round040_source/source_tables.py')]},
         'seconds':time.time()-started,
         'scope':'Deterministic true-source 1D and exact abstract channel; no Gaussian/CLT substitution or tail MC'}
    target.write_text(json.dumps(out,indent=2)+'\n');print('seconds',out['seconds'],flush=True)
if __name__=='__main__':main()
