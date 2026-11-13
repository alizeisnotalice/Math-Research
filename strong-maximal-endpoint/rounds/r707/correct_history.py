"""Replace affected r706 high-dimensional diagnostics using identical search bounds."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import json,sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from scipy.optimize import minimize_scalar
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'r706'))
import absolute
from stable_gamma import Special

def run(old):
    absolute.sp=Special(64);n=old['n'];z,H,eta=absolute.normal(n)
    identical=absolute.soft(n,old['c_soft'])
    sol=minimize_scalar(lambda c:-absolute.soft(n,c)['centered'],bounds=old['search_bounds'],method='bounded',options={'xatol':2e-10,'maxiter':120})
    best=absolute.soft(n,float(sol.x));dn=absolute.SIGMA*n**.5*z
    return dict(n=n,previous_c_soft=old['c_soft'],old_centered=old['centered'],same_c_stable_centered=identical['centered'],
        same_c_root_correction=identical['centered']-old['centered'],new_c_soft=float(sol.x),optimizer_success=bool(sol.success),
        new_centered=best['centered'],root_residual=best['root_residual'],
        corrected_constant_diagnostic=best['centered']+dn-absolute.CQ*z*z)
if __name__=='__main__':
    rows=[r for r in json.loads((HERE.parent/'r706/absolute.json').read_text()) if r['n']>=1048576]
    with ProcessPoolExecutor(3) as pool:out=list(pool.map(run,rows))
    (HERE/'corrected_r706.json').write_text(json.dumps(out,indent=2)+'\n')
    for r in out:print(r,flush=True)
