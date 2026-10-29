"""Frozen deterministic Gaussian-channel coupling probe, not a cube model."""
from pathlib import Path
from statistics import NormalDist
import hashlib,json,math,sys,time
import numpy as np

HERE=Path(__file__).resolve().parent
ND=NormalDist()
def cdf(x):
    a=np.asarray(x)
    return np.fromiter((.5*math.erfc(-float(t)/math.sqrt(2))for t in a.flat),float,a.size).reshape(a.shape)
def interval(lo,hi):
    # Exact real Gaussian formula with a stable upper-tail branch.
    lo,hi=np.broadcast_arrays(lo,hi)
    return np.where(lo>=0,cdf(-lo)-cdf(-hi),cdf(hi)-cdf(lo))
def lse(a,axis=0):
    m=np.max(a,axis=axis,keepdims=True)
    return np.squeeze(m+np.log(np.exp(a-m).sum(axis=axis,keepdims=True)),axis=axis)
def profile(logg,y,dy):
    logw=-.5*y*y-.5*math.log(2*math.pi)+math.log(dy)
    order=np.argsort(logg)[::-1]
    logtail=np.logaddexp.accumulate(logw[order])
    score=logg[order]+logtail
    i=int(np.argmax(score))
    return {'coefficient':float(math.exp(score[i])), 'log_threshold':float(logg[order[i]]),
            'boundary_y':float(y[order[i]]),
            'output_mass_quadrature':float(math.exp(float(lse(logg+logw))))}
def model(tau,c,R,dy):
    s=math.sqrt(1-tau*tau)
    cuts=np.array([-math.inf]+[ND.inv_cdf(j/R)for j in range(1,R)]+[math.inf])
    med=np.array([ND.inv_cdf((j+.5)/R)for j in range(R)])
    shifts=c-tau*med
    y=np.arange(-12,c+12+dy/2,dy)
    total=np.zeros_like(y)
    priorlog=[]
    # RN(gate) = RN(component) * conditional segment probability.
    for b,t in enumerate(shifts):
        z=y-t
        pb=interval((cuts[b]-tau*z)/s,(cuts[b+1]-tau*z)/s)
        total+=np.exp(t*y-t*t/2)*pb
        priorlog.append(t*y-t*t/2-math.log(R))
    assert total.min()>0
    gate=profile(np.log(total),y,dy)
    prior=profile(lse(np.array(priorlog)),y,dy)
    limit_gate=profile(-math.log(s)+y*y/2-(y-c)**2/(2*s*s),y,dy)
    sp=math.sqrt(1+tau*tau)
    limit_prior=profile(-math.log(sp)+y*y/2-(y-c)**2/(2*sp*sp),y,dy)
    # Direct Gaussian input-output mixture density equals the RN formula.
    iy=[int(len(y)*q)for q in [.25,.5,.75]]
    check=[]
    nodes,weights=np.polynomial.legendre.leggauss(64)
    for index in iy:
        yy=float(y[index]);direct=0.
        for b,t in enumerate(shifts):
            # Integrate source prior quantile u exactly in probability space.
            u=(b+(nodes+1)/2)/R
            xx=np.array([ND.inv_cdf(float(v))for v in u])
            density=np.exp(-.5*((yy-t-tau*xx)/s)**2)/(s*math.sqrt(2*math.pi))
            direct+=float(weights@density)/(2*R)
        rn=total[index]*math.exp(-yy*yy/2)/math.sqrt(2*math.pi)
        check.append({'y':yy,'RN_density':rn,'source_quadrature_density':direct,
                      'absolute_difference':abs(rn-direct)})
    return {'tau':tau,'c':c,'R':R,'dy':dy,'gate':gate,'prior':prior,
            'ratio_of_empirical_fixed_scale_weak_functionals':gate['coefficient']/prior['coefficient'],
            'continuum_gate':limit_gate,'continuum_prior':limit_prior,
            'continuum_ratio':limit_gate['coefficient']/limit_prior['coefficient'],
            'source_quadrature_checks':check}
def main():
    out=HERE/'results.json'
    if out.exists():raise RuntimeError('Refusing replacement')
    plan={'tau':[.15,1/3],'c':[4.,8.],'R':[2,4,8,16,32,64], 'dy':.002,
          'confirmation':[{'tau':1/3,'c':8.,'R':64,'dy':.001}],
          'script_SHA256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'scope':'Deterministic Gaussian channel prototype; no true cube averages or C_n inference.'}
    (HERE/'registration.json').write_text(json.dumps(plan,indent=2)+'\n')
    started=time.time()
    rows=[model(t,c,r,plan['dy'])for t in plan['tau']for c in plan['c']for r in plan['R']]
    confirm=[model(**case)for case in plan['confirmation']]
    data={'plan_SHA256':hashlib.sha256((HERE/'registration.json').read_bytes()).hexdigest(),
          'numpy':np.__version__,'python':sys.version,'cases':rows,'confirmation':confirm,
          'seconds':time.time()-started,'interval_certified':False}
    out.write_text(json.dumps(data,indent=2)+'\n')
    for row in rows:
        print(row['tau'],row['c'],row['R'],row['gate']['coefficient'],row['prior']['coefficient'],row['continuum_ratio'])
    print('seconds',data['seconds'])
if __name__=='__main__':main()
