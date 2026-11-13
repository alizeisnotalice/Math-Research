"""Rerun distribution diagnostics through the positive-integral tail implementation."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import json,sys,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scipy import special as sp
import refined
from stable_gamma import Special
HERE=Path(__file__).resolve().parent

def run(row):
    start=time.time();n=row['n'];c=row['c_eta'];base=refined.full_soft(n,c,'stable64');check=refined.full_soft(n,c,'stable32')
    out=dict(n=n,c_eta=c,stable_centered=base['centered'],stable_root_residual=base['residual'],
             refined_soft_error=base['refined_soft_error'],scaled_refined_soft_error=np.log(n)*base['refined_soft_error'],
             stable_order_stability=abs(check['centered']-base['centered']),
             original_minus_stable=row['centered']-base['centered'],
             edgeworth4_minus_stable=row.get('edgeworth4_centered',base['centered'])-base['centered'])
    if 'nonlinear_centered' in row:
        for p in ['r700','r692','r704']:sys.path.insert(0,str(HERE.parent/p))
        import joint,refine,profile
        special=Special(64)
        for mod in [joint,refine,profile]:mod.MU=c;mod.P=sp.expit(-c);mod.sp=special
        old={r['n']:r for r in json.loads((HERE.parent/'r700/refined.json').read_text())}
        ramp=refine.run(old[n]);ev,_,_,_=joint.evaluator(n,1e-16)
        rc=ev(ramp['a'],ramp['b'])[0]+n*(c-refined.MU)
        vals=[profile.run((n,o,18)) for o in [256,512]];sc=vals[-1]['score']-n*refined.MU
        out.update(ramp_centered=rc,ramp_success=ramp['success'],ramp_residual=ramp['residual'],
              nonlinear_centered=sc,nonlinear_quadrature_stability=abs(vals[1]['score']-vals[0]['score']),
              nonlinear_mean_residual=vals[-1]['average_residual'],
              original_nonlinear_minus_corrected=row['nonlinear_centered']-sc,
              nonlinear_minus_resummed=sc-row['resummed_score_centered'],
              nonlinear_minus_expanded=sc-row['expansion_score_centered'],
              scaled_expansion_remainder=row['z']**2*(sc-row['expansion_score_centered']),
              scaled_slab_remainder=np.log(n)*(rc-base['centered']-row['resummed_slab_gain']))
    out['seconds']=time.time()-start
    return out
if __name__=='__main__':
    rows=json.loads((HERE/'refined.json').read_text());out=[]
    with ProcessPoolExecutor(12) as pool:
        for future in as_completed([pool.submit(run,r) for r in rows]):
            item=future.result();out.append(item);print(json.dumps(item),flush=True)
            (HERE/'corrected.json').write_text(json.dumps(sorted(out,key=lambda r:r['n']),indent=2)+'\n')
