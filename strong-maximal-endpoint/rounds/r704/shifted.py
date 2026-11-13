"""Each worker runs one model at a time and sets all imported model constants."""
import json,sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'r700'))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'r692'))
import joint
import refine
import profile
ORIGINAL_MU=joint.MU
def run(job):
    n,t,row=job
    c=ORIGINAL_MU+t*np.sqrt(np.log(n)/n);p=1/(1+np.exp(c))
    for module in [joint,refine,profile]:
        module.MU=c;module.P=p
    ramp=refine.run(row)
    ev,_,_,_=joint.evaluator(n,1e-16)
    ramp_score=ev(ramp['a'],ramp['b'])[0]+c*n
    values=[profile.run((n,order,18)) for order in [256,512]]
    correction=values[-1]['score']-ramp_score
    return dict(n=n,t=t,c=c,ramp_a=ramp['a'],ramp_b=ramp['b'],
                ramp_success=ramp['success'],ramp_residual=ramp['residual'],
                ramp_score=float(ramp_score),nonlinear_score=float(values[-1]['score']),
                gap=float(correction),scaled_gap=float(np.log(n)*correction),
                nonlinear_order_stability=float(abs(values[-1]['score']-values[0]['score'])),
                nonlinear_residual=values[-1]['average_residual'])
if __name__=='__main__':
    old={r['n']:r for r in json.loads(Path('../r700/refined.json').read_text())}
    jobs=[(n,t,old[n]) for n in [64,1024,16384,262144] for t in [-1,-np.sqrt(ORIGINAL_MU)/2,0,.5]]
    with ProcessPoolExecutor(12) as pool:out=list(pool.map(run,jobs))
    Path('shifted.json').write_text(json.dumps(out,indent=2)+'\n')
    for row in out:print(row,flush=True)
