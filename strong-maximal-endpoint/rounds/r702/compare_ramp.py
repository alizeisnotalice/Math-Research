"""Interval root enclosure for one explicit ramp, using full finite mixtures."""
import json,math,sys
from pathlib import Path
from mpmath import mp,iv
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'r692'))
from joint import evaluator,MU
from certify import setup,distribution
def gammaQ(x):
    if x.b<=0:return iv.mpf(1)
    assert x.a>0
    poly=iv.mpf(1)
    for j in range(63,0,-1):poly=1+x*poly/j
    return iv.exp(-x)*poly
def main():
    c,wt,t0=setup();mp.dps=80;iv.dps=55
    row=json.loads(Path('../r700/refined.json').read_text())[0]
    a=iv.mpf(str(row['a']));b=iv.mpf(str(row['b']));L=iv.log(64)
    beta=1+a*iv.sqrt(L/64);w=(1+a*a)*L/2+iv.log(L)+b
    z=64/beta-1
    I=64/beta*iv.exp(z)*math.factorial(63)*gammaQ(z)/z**64
    K=1+iv.exp(-w)*(I-1)
    ev,_,_,_=evaluator(64,1e-16)
    approx=64*MU+ev(row['a'],row['b'])[0]+float(iv.log(K).mid)
    lo=mp.mpf(str(approx))-mp.mpf('1e-7');hi=mp.mpf(str(approx))+mp.mpf('1e-7')
    def residual(q):
        q=iv.mpf(mp.nstr(q,75));t1=q+w;rate=1-beta/64
        A=sum((ww*iv.exp(-beta*(t1-c*k)/64)*rate**(-64)*gammaQ(rate*(t1-c*k)) for k,ww in enumerate(wt)),iv.mpf(0))
        _,sq=distribution(q,c,wt);_,st=distribution(t1,c,wt)
        return A+sq-st-1
    assert residual(lo)>0 and residual(hi)<0
    for _ in range(18):
        mid=(lo+hi)/2;r=residual(mid)
        if r.a>0:lo=mid
        elif r.b<0:hi=mid
        else:break
    score=iv.mpf([mp.nstr(lo,75),mp.nstr(hi,75)])-iv.log(K)
    out=dict(n=64,a=str(row['a']),b=str(row['b']),root_bracket=[mp.nstr(lo,70),mp.nstr(hi,70)],score=str(score))
    Path('ramp_certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
if __name__=='__main__':main()
