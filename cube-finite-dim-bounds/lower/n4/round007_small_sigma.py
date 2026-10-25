"""Two finite actual joint checks of the best dual-satellite directional model."""
from pathlib import Path
import json,time,numpy as np
from round005_exposure import A,geometry
from round007_three_orbit import evaluate,face_keys
from round007_actual_search import save_candidate,OLD
P=Path(__file__).resolve().parent

def main():
    st=time.monotonic();m=json.loads((P/'round007_model_results.json').read_text());best=max(m['rows'],key=lambda r:r['prediction']);patches,_,_=geometry();p1=np.array(best['phi1']);p2=np.array(best['phi2']);d=np.array(best['direction']);pool,meta=face_keys([p1,p2],patches);rows=[]
    for sigma in [1e-7,1e-8]:
        row,keys,selected=evaluate(A*(1+sigma*d),p1,p2,sigma,best['share1'],pool=pool,power=10,seed=20261025711);row['joint_model_prediction']=best['prediction'];rows.append(row);print(row,flush=True)
        if row['ratio']>OLD+1e-11 and row['increment_quotient']>1e-5:save_candidate(row,selected,'small-sigma');raise RuntimeError('positive finite union: certificate priority')
    (P/'round007_small_sigma_results.json').write_text(json.dumps(dict(status='two actual finite joint checks, not strict derivative proof',case=best,face_metadata=meta,rows=rows,elapsed=time.monotonic()-st),indent=2)+'\n')
if __name__=='__main__':main()
