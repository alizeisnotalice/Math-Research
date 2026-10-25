"""Registered final union refresh; retain256-center and all original witness keys."""
from pathlib import Path
import json,time,numpy as np
from scipy.stats import qmc
from round005_exposure import geometry
from round005_calibrate import face_keys
from round004_satellite_search import capture
from round006_actual import actual,save_candidate,OLD
P=Path(__file__).resolve().parent

def main():
    st=time.monotonic();data=json.loads((P/'round006_actual_results.json').read_text());screen=json.loads((P/'round006_screen_results.json').read_text());patches,_,_=geometry();out=[]
    for result,initial in zip(data['runs'],screen['runs']):
        x=np.array(result['x']);a=x[:4];phi=x[4:8];eps=x[8];pool,meta=face_keys(np.array(initial['phi']),patches)
        pool|=capture(a,phi,eps,qmc.Sobol(4,scramble=True,seed=20261025606).random_base2(8))
        row,keys=actual(a,phi,eps,pool,power=12,seed=20261025609);row['old_optimizer_quotient']=result['quotient'];out.append(row);print('FINAL REFRESH',row,flush=True)
        if row['ratio']>OLD+1e-7:save_candidate(row,keys);raise RuntimeError('positive actual union: certificate priority')
    (P/'round006_final_results.json').write_text(json.dumps(dict(status='full4D final finite unions with explicit singleton and old/new captures',rows=out,elapsed=time.monotonic()-st),indent=2)+'\n')
if __name__=='__main__':main()
