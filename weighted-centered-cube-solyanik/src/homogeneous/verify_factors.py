"""Independent adaptive-quadrature cross-checks, not interval certificates."""
import json
import math
from pathlib import Path
from scipy.integrate import quad, dblquad

ROOT=Path(__file__).resolve().parent
PS=(1,1.25,1.5,2,3,4,8,16,'inf')
ORDERS=[32,64,128,256,512,1024]


def independent(n,p,tau):
    volume=2**n*math.gamma(1+1/p)**n/math.gamma(1+n/p)
    exponent=(tau-n)/p
    if n==2:
        integral,error=quad(lambda y:(1+y**p)**exponent,0,1,
                            epsabs=1e-11,epsrel=1e-11)
    else:
        integral,error=dblquad(lambda z,y:(1+y**p+z**p)**exponent,
                               0,1,lambda y:0,lambda y:1,
                               epsabs=1e-11,epsrel=1e-11)
    kappa=n*2**n*integral/tau
    value=(n/tau)*volume**(1-tau/n)*(2**n)**(tau/n)/kappa
    return value,error


def main():
    rows=json.loads((ROOT/'grid.json').read_text())
    expected={(n,str(p),j) for n in (1,2,3) for p in PS for j in range(1,128)}
    assert len(rows)==3429 and {(r['n'],str(r['p']),r['j']) for r in rows}==expected
    assert all(r['status']=='completed' and r['orders']==ORDERS and
               len(r['coefficients'])==6 and
               all(math.isfinite(v) and v>0 for v in r['coefficients']) for r in rows)
    assert all(r['coefficients']==[1.0]*6 for r in rows if r['n']==1 or r['p']=='inf')
    selected=[]
    for row in rows:
        if row['n']==1 or row['p']=='inf' or row['j'] not in (1,32,64,96,127):
            continue
        n,p,tau=row['n'],row['p'],row['tau']
        value,estimated_error=independent(n,p,tau)
        difference=abs(value-row['coefficients'][-1])
        assert estimated_error<1e-9 and difference<2e-9
        selected.append(dict(n=n,p=p,j=row['j'],tau=tau,factor=value,
                             batch_difference=difference,quad_error_estimate=estimated_error))
    assert len(selected)==80
    # Three independent closed forms prevent a common normalization error.
    exact_cases=[(2,2,1,math.sqrt(math.pi)/(2*math.asinh(1))),
                 (2,1,1,1/(math.sqrt(2)*math.log(2))),
                 (3,1,1,6**(-2/3)/math.log(4/3))]
    for n,p,tau,closed_form in exact_cases:
        value,_=independent(n,p,tau)
        assert abs(value-closed_form)<1e-11
        assert value>1
    result=dict(status='passed',complete_rows=len(rows),adaptive_cases=len(selected),
                closed_form_cases=len(exact_cases),
                max_difference=max(r['batch_difference'] for r in selected),
                max_error_estimate=max(r['quad_error_estimate'] for r in selected))
    (ROOT/'selected.json').write_text(json.dumps(selected,indent=2))
    (ROOT/'verification.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result))


if __name__=='__main__':main()
