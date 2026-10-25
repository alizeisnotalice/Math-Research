"""Two registered true-union scale-separation calibrations, no optimization."""
from pathlib import Path
import json,time,numpy as np
from round005_exposure import A,geometry
from round005_calibrate import face_keys
from round006_actual import actual
P=Path(__file__).resolve().parent

def main():
    st=time.monotonic();d=json.loads((P/'round006_mixed_breakpoint_results.json').read_text());best=max(d['rows'],key=lambda r:r['prediction']);patches,_,_=geometry();phi=np.array(best['phi']);direction=np.array(best['direction']);pool,meta=face_keys(phi,patches);rows=[]
    for eps in [1e-7,1e-8]:
        row,keys=actual(A*(1+eps*direction),phi,eps,pool,power=10,seed=20261025608);row['prediction']=best['prediction'];rows.append(row);print(row,flush=True)
    (P/'round006_small_epsilon_results.json').write_text(json.dumps(dict(status='two numerical scale-separation checks, not exact derivative proof',best_direction=best,face_key_metadata=meta,rows=rows,elapsed=time.monotonic()-st),indent=2)+'\n')
if __name__=='__main__':main()
