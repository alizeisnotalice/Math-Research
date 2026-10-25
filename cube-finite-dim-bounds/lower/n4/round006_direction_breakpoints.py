"""Finite breakpoint maximization of the separable period-direction model."""
from pathlib import Path
from fractions import Fraction as F
import time,json
import numpy as np
from round005_exposure import A,geometry,phase_model,prediction
P=Path(__file__).resolve().parent

def axis_value(part,d,J,old):
    w,cg,gp,gm,R,mp,mm,fp,fm,seam=part
    plus=np.max(np.where(mp,gp-R*d,-np.inf),axis=1);minus=np.max(np.where(mm,gm-R*d,-np.inf),axis=1)
    born=np.maximum(0,np.max(np.where(fp|fm,cg-R*d,-np.inf),axis=1));inside=np.maximum(plus,minus)
    gain=inside+np.maximum(0,born-np.maximum(inside,-born))
    L=np.where(seam,plus,0);right=np.where(seam,-minus,0)
    gap=np.maximum(0,right-L)-np.maximum(0,np.minimum(right,born)-np.maximum(L,-born))
    return float(J*(w@np.where(seam,-gap,gain))+old*d)

def candidates(part):
    w,cg,gp,gm,R,mp,mm,fp,fm,seam=part;cuts={-.5,0.,.5}
    patterns=set()
    for i in range(len(w)):
        lines=set()
        for g,mask in [(gp,mp),(gm,mm),(cg,fp|fm)]:
            for k in np.flatnonzero(mask[i]):lines.add((round(float(g[i,k]),14),round(float(R[i,k]),14)))
        patterns.add(tuple(sorted(lines)))
    for lines in patterns:
        for g,r in lines:
            if r and -.5<=g/r<=.5:cuts.add(round(g/r,13))
        for g,r in lines:
            for h,s in lines:
                for num,den in [(g-h,r-s),(g+h,r+s)]:
                    if abs(den)>1e-13:
                        x=num/den
                        if -.5<=x<=.5:cuts.add(round(x,13))
    return sorted(cuts),len(patterns)

def main():
    start=time.monotonic();screen=json.loads((P/'round006_screen_results.json').read_text());patches,_,geom=geometry();old=geom['baseline_joint_ratio'];J=float(np.prod(A));rows=[]
    for row in screen['runs']:
        model=phase_model(np.array(row['phi']),patches);direction=[];axes=[]
        for j in range(4):
            cuts,npats=candidates(model[j]);vals=[axis_value(model[j],d,J,old) for d in cuts];k=int(np.argmax(vals));direction.append(cuts[k]);axes.append(dict(candidates=len(cuts),line_patterns=npats,maximizer=cuts[k],value=vals[k]))
        value=prediction(model,np.array(direction),old)+float(F(row['theta']))
        out=dict(source_masses=row['source_masses'],phi=row['phi'],theta=row['theta'],old_optimizer_prediction=row['prediction'],direction=direction,prediction=value,axes=axes,numerical_breakpoint_rounding=13)
        rows.append(out);print('BREAKPOINT',out,flush=True)
        if value>1e-6:(P/'round006_direction_packet.json').write_text(json.dumps(out,indent=2)+'\n')
    result=dict(status='finite numerical affine-breakpoint maximization; not exact derivative certificate',rows=rows,elapsed=time.monotonic()-start)
    (P/'round006_breakpoint_results.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
