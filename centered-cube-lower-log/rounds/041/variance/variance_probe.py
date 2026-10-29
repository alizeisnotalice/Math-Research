"""Small exact continuous-bin box checks and true common-D 1D diagnostics.

Does not read/transform tail payloads and creates no new statistical CI.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import math
import sys
import time
import numpy as np
from source_tables import Source
from validate_source import reference_events

HERE=Path(__file__).resolve().parent


def rational_case(p,b,plus,center):
    k,r,w,v=4,4,4,8
    delta=F(1,w*k*v);h=delta/6;ell=F(1,6);length=1+ell
    atoms=[]
    for j,label in enumerate([0,1,2,3]):
        u=(F(label,4)+(F(b)+F(1,2))/r)%1 if j<2 else (F(label,4)+(F(1,2)-b)/r)%1
        if (j<2)==plus:atoms.append((F(j)+u)/k)
    extended=[a+s for s in [-1,0,1,2]for a in atoms]
    ca,cz=(F(0),F(1,2))if plus else(F(1,2),F(1))
    left,right=center-h/2,center+h/2
    cuts={left,right}
    for a in extended:
        for side in [-1,1]:
            for ds in [-1,1]:
                x=a-F(1,2)+side*ell/2+ds*delta/2
                if left<x<right:cuts.add(x)
    for edge in [ca,cz]:
        for shift in [-1,0,1,2]:
            for side in [-1,1]:
                x=edge+shift-F(1,2)+side*ell/2
                if left<x<right:cuts.add(x)
    clip=lambda x:max(F(0),min(F(1),x))
    def kernel(x):
        y=x+F(1,2)
        overlap=sum((max(F(0),min(y+ell/2,cz+s)-max(y-ell/2,ca+s))
                     for s in [-1,0,1,2]),F(0))
        peak=sum((clip((y-a+ell/2+delta/2)/delta)-clip((y-a-ell/2+delta/2)/delta)
                  for a in extended),F(0))
        return (F(1,2)+(1-p)*overlap+p/k*peak)/length
    cuts=sorted(cuts);values=[kernel(x)for x in cuts]
    mean=sum(((z-a)*(u+t)/2 for a,z,u,t in zip(cuts,cuts[1:],values,values[1:])),F(0))/h
    second=sum(((z-a)*(u*u+u*t+t*t)/3 for a,z,u,t in zip(cuts,cuts[1:],values,values[1:])),F(0))/h
    variance=second-mean*mean
    y=center+F(1,2);support=(delta+h)/2
    perend=[sum(abs(a-(y+side*ell/2))<=support for a in extended)for side in [-1,1]]
    assert max(perend)<=1
    partial=sum(perend)
    lip=((1-p)+p/k*partial/delta)/length
    old=lip*lip*h*h/4;new=lip*lip*h*h/12
    assert variance<=new<=old and old==3*new
    return {"p":str(p),"table":b,"class":"plus"if plus else"minus","center":str(center),
        "delta":str(delta),"h":str(h),"N_partial":partial,
        "mean":str(mean),"variance":str(variance),"local_Lipschitz_cap":str(lip),
        "old_range_quarter_cap":str(old),"new_Lipschitz_twelfth_cap":str(new),
        "variance_le_new_cap_exact":True,"affine_cap_equality":variance==new}


def rational_controls():
    rows=[]
    for p in [F(1,2),F(1)]:
        for b in [0,3]:
            for plus in [True,False]:
                j=0 if plus else 2;label=j
                phase=(F(b)+F(1,2))/4 if plus else(F(1,2)-b)/4
                atom=(F(j)+(F(label,4)+phase)%1)/4
                centers=[(atom+side*F(1,12)-F(1,2))%1 for side in [-1,1]]
                centers.extend([F(1,7),F(1,12)])
                rows.extend(rational_case(p,b,plus,c)for c in centers)
    assert any(x["affine_cap_equality"]and F(x["variance"])>0 for x in rows)
    return rows


def true_box_controls():
    rows=[];nodes,weights=np.polynomial.legendre.leggauss(24);length=7/6
    for r in [4,32]:
        src=Source(r)
        for b in [0,r-1]:
            evaluated=src.evaluate(length,b);events=reference_events(src,length,b)
            for plus in [True,False]:
                name="plus"if plus else"minus"
                partial=evaluated["transition_count_"+name]
                for n in [0,1,2]:
                    indices=np.flatnonzero(partial==n)
                    if not len(indices):continue
                    q=int(indices[0]);left,right=q/src.Q,(q+1)/src.Q
                    cuts=np.unique(np.r_[left,events[(left<events)&(events<right)],right])
                    x=np.concatenate([(a+z)/2+(z-a)/2*nodes for a,z in zip(cuts,cuts[1:])])
                    ww=np.concatenate([(z-a)*weights/2 for a,z in zip(cuts,cuts[1:])]);ww/=ww.sum()
                    z=src.point_class_kernel(x,length,b,plus)
                    mean=float(ww@z);variance=float(ww@((z-mean)**2))
                    osc=evaluated["range_"+name][q]
                    old=osc**2/4;new=osc**2/12
                    assert variance<=new+1e-20
                    rows.append({"R":r,"table":b,"class":name,"bin":q,"N_partial":n,
                        "reference_mean":mean,"reference_variance":variance,
                        "old_quarter_cap":old,"candidate_twelfth_cap":new,
                        "variance_over_candidate_cap":variance/new,
                        "normalizer":evaluated["normalizer"]})
    assert any(x["N_partial"]==0 and x["reference_variance"]>0 for x in rows)
    return rows


def main():
    target=HERE/"CANDIDATE_VARIANCE_RESULTS.json"
    if target.exists():raise RuntimeError("Refusing replacement.")
    started=time.time();rational=rational_controls();actual=true_box_controls()
    out={"scope":"small continuous-uniform-bin variance checks only; no pilot transformation or CI",
        "script_SHA256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_SHA256":hashlib.sha256((HERE/"source_tables.py").read_bytes()).hexdigest(),
        "python":sys.executable,"numpy":np.__version__,"rational_true_box_cases":rational,
        "common_D_true_box_numerical_cases":actual,"ordinary_FP_not_interval_certified":True,
        "seconds":time.time()-started}
    target.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({"rational_cases":len(rational),"nonzero_affine_equalities":sum(
        x["affine_cap_equality"]and F(x["variance"])>0 for x in rational),
        "common_D_cases":len(actual),"max_actual_over_candidate":max(
            x["variance_over_candidate_cap"]for x in actual),"seconds":out["seconds"]},indent=2))


if __name__=="__main__":main()
