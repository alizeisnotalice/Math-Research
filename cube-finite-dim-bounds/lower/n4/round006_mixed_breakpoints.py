"""No new phases: finite breakpoint period maximization at stored positive-theta phases."""
from pathlib import Path
from fractions import Fraction as F
import json,time
import numpy as np
from round005_exposure import A,geometry,phase_model,prediction
from round006_direction_breakpoints import candidates,axis_value
P=Path(__file__).resolve().parent

def main():
    st=time.monotonic();data=json.loads((P/'round006_screen_results.json').read_text());patches,_,geom=geometry();J=float(np.prod(A));old=geom['baseline_joint_ratio'];rows=[];positive=None
    selected=sorted((r for r in data['scored'] if F(r['theta'])>0),key=lambda r:-r['best_prediction'])
    for i,row in enumerate(selected):
        model=phase_model(np.array(row['phi']),patches);direction=[];axes=[]
        for j in range(4):
            cuts,pats=candidates(model[j]);vals=[axis_value(model[j],d,J,old) for d in cuts];k=int(np.argmax(vals));direction.append(cuts[k]);axes.append(dict(candidates=len(cuts),maximizer=cuts[k]))
        val=prediction(model,np.array(direction),old)+float(F(row['theta']));out=dict(source_masses=row['source_masses'],phi=row['phi'],theta=row['theta'],direction=direction,prediction=val,axes=axes);rows.append(out)
        if i%16==0:print(i,out,flush=True)
        if val>1e-6:positive=out;(P/'round006_direction_packet.json').write_text(json.dumps(out,indent=2)+'\n');break
    result=dict(status='finite numerical breakpoint direction maximization, exact theta; no new phases',registered_phase_limit=171,processed=len(rows),positive=positive,rows=rows,elapsed=time.monotonic()-st)
    (P/'round006_mixed_breakpoint_results.json').write_text(json.dumps(result,indent=2)+'\n');print('FINAL',len(rows),'best',max(r['prediction'] for r in rows),'positive',positive,flush=True)
if __name__=='__main__':main()
