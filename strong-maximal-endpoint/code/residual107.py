"""One-step Gaussian and second-Hermite comparison using exact ex-Gaussian convolution."""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse,json,math
import mpmath as m


def evaluate(spec):
    j,scale,K=spec
    values=[]
    for wp in (70,140):
        m.mp.dps=wp+int(j*1.3)+45
        t=m.mpf(2)**j;mu=1+m.lambertw(1/m.e);sig=m.sqrt(mu);a=(9-mu)/24
        u=m.mpf(scale)*m.sqrt(m.log(t));b=u*sig*m.sqrt(t);sdev=sig*m.sqrt(t-1)
        c=mu+mu*u/(2*sig*m.sqrt(t));prob=1/(1+m.exp(c))
        def cdf(x):return m.erfc(-x/m.sqrt(2))/2
        def density(x,v):return m.exp(-(x/v)**2/2)/(m.sqrt(2*m.pi)*v)
        def exp_integral(x):
            f=density(x,sdev);df=-x*f/sdev**2
            conv=m.exp(-x+sdev**2/2)*cdf(x/sdev-sdev)
            return (cdf(x/sdev) if u<=0 else -cdf(-x/sdev))-conv-a*sdev**2*(df-f+conv)
        lhs=(1-prob)*exp_integral(b+mu)+prob*exp_integral(b+mu-c)
        f=density(b,sig*m.sqrt(t));rhs=(cdf(u) if u<=0 else -cdf(-u))-a*(u*u-1)*f-K*m.log(t/(t-1))*f
        value=t*(lhs-rhs)/f;values.append(value)
    error=abs(values[1]-values[0]);assert error<m.mpf('1e-45'),(spec,error)
    return dict(log2_t=j,u_scale=scale,K=K,scaled_residual=m.nstr(values[1],35),limiting_value=m.nstr(K-(mu-1)/8,35),precision_difference=m.nstr(error,8))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('output',type=Path);a=ap.parse_args();cases=[(j,z,k) for j in (12,16,24,32,48,64,96,128,192,256) for z in ('-2','-1','-0.5','0','0.5','1','2') for k in (0,1)]
    with ProcessPoolExecutor(12) as pool:rows=list(pool.map(evaluate,cases))
    a.output.write_text(json.dumps(dict(cases=len(rows),pass_all=True,scope='Finite matched-action residual without projected initial correction. Analytic all-action and all-state estimates are in the proof.',rows=rows),indent=2)+'\n')
