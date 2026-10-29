"""Frozen deterministic 1D controls; no high-dimensional simulation or randomness."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import math
import sys
import time
import numpy as np
from source_tables import Source, field_labels

HERE=Path(__file__).resolve().parent


def selected_bins(src,L,b):
    ell=L-1;points=[]
    for j in [0,63,126,127,191,254]:
        for side in [-1,1]:
            x=(src.atoms[b,j]+side*ell/2-.5)%1
            q=int(math.floor(x*src.Q))
            points.extend([(q-1)%src.Q,q,(q+1)%src.Q])
    for edge in [0,src.a,1]:
        for side in [-1,1]:
            q=int(math.floor(((edge+side*ell/2-.5)%1)*src.Q))
            points.extend([(q-1)%src.Q,q,(q+1)%src.Q])
    return sorted(set(points+[0,src.Q//7,src.Q//2,src.Q-1]))


def cuts_for_bin(src,L,b,q):
    left,right=q/src.Q,(q+1)/src.Q
    cuts=[left,right]
    ell=L-1
    for plus in [True,False]:
        cuts.extend(x-.5 for x in src.background_breakpoints(ell,plus)
                    if left<x-.5<right)
    for atom in src.atoms[b]:
        for shift in [-1,0,1,2]:
            for side in [-1,1]:
                for ds in [-1,1]:
                    x=atom+shift-.5+side*ell/2+ds*src.delta/2
                    if left<x<right:cuts.append(x)
    return np.unique(cuts)


def rational_peak_mean(src,L,b,q,plus):
    ell=F(L)-1;delta=F(1,src.w*src.K*src.v)
    left,right=F(q,src.Q),F(q+1,src.Q)
    base=[]
    for j in range(src.K):
        if (j<src.cut)!=plus:continue
        label=(int(src.labels[j])+(b if j<src.cut else -b))%4
        base.append((F(j)+F(2*label+1,8))/src.K)
    atoms=[a+s for s in [-1,0,1,2]for a in base]
    cuts={left,right}
    for atom in atoms:
        for side in [-1,1]:
            for ds in [-1,1]:
                x=atom-F(1,2)+side*ell/2+ds*delta/2
                if left<x<right:cuts.add(x)
    cuts=sorted(cuts)
    clip=lambda x:max(F(0),min(F(1),x))
    def value(x):
        y=x+F(1,2)
        return sum((clip((y-a+ell/2+delta/2)/delta)
                    -clip((y-a-ell/2+delta/2)/delta)for a in atoms),F(0))
    values=[value(x)for x in cuts]
    mean=sum(((z-a)*(u+v)/2 for a,z,u,v in zip(cuts,cuts[1:],values,values[1:])),F(0))/(right-left)
    return mean


def main():
    target=HERE/"SOURCE_VALIDATION.json"
    if target.exists():raise RuntimeError("Refusing replacement.")
    started=time.time();src=Source();L=7/6
    nodes,weights=np.polynomial.legendre.leggauss(24)
    table_rows=[];rational_rows=[]
    for b in range(src.R):
        evaluated=src.evaluate(L,b);selected=selected_bins(src,L,b)
        ref_rows=[]
        for q in selected:
            cuts=cuts_for_bin(src,L,b,q)
            x=np.concatenate([(a+z)/2+(z-a)/2*nodes for a,z in zip(cuts,cuts[1:])])
            w=np.concatenate([(z-a)/2/src.binwidth*weights for a,z in zip(cuts,cuts[1:])])
            for plus in [True,False]:
                key="plus"if plus else"minus"
                z=src.point_class_kernel(x,L,b,plus)
                mean=float(w@z);variance=float(w@((z-mean)**2))
                approx=float(evaluated["mean_"+key][q])
                error=abs(approx-mean)
                point_values=src.point_class_kernel(np.unique(np.concatenate([x,cuts])),L,b,plus)
                observed_range=float(np.ptp(point_values))
                range_cap=float(evaluated["range_"+key][q])
                assert error<evaluated["mean_error_abs"]+2e-13
                assert observed_range<=range_cap+2e-13
                assert variance<=range_cap**2/4+1e-20
                ref_rows.append({"bin":q,"class":key,"reference_mean":mean,
                    "computed_mean":approx,"absolute_mean_error":error,
                    "reference_variance":variance,"observed_range":observed_range,
                    "range_cap":range_cap})
        for q in selected[:4]+selected[-4:]:
            for plus in [True,False]:
                exact=rational_peak_mean(src,F(7,6),b,q,plus)
                approx=src._peak_bin_mean_and_transition_count(
                    np.array([src.bins_center[q]+.5]),L-1,b,plus)[0][0]
                error=abs(float(exact)-approx)
                # Count error: atom positions/lengths first represented in FP.
                assert error<2e-10
                rational_rows.append({"b":b,"bin":q,"class":"plus"if plus else"minus",
                    "exact_extra_peak_count_mean":str(exact),"computed":float(approx),
                    "count_absolute_error":error,"kernel_error":error*src.p/(src.K*L)})
        table_rows.append({"b":b,"selected_bin_count":len(selected),
            "normalizer":evaluated["normalizer"],
            "max_relative_variance_class_bound":float(evaluated["relative_variance_class_max"].max()),
            "analytic_mean_error_abs":evaluated["analytic_mean_error_abs"],
            "floating_point_guard_abs":evaluated["floating_point_guard_abs"],
            "max_observed_mean_error":max(x["absolute_mean_error"]for x in ref_rows),
            "max_observed_range_over_cap":max(x["observed_range"]/x["range_cap"]for x in ref_rows),
            "reference_rows":ref_rows})
    # Check copied field_labels text and actual label values, without importing
    # a module in another agent's directory or creating its __pycache__.
    import ast,inspect
    old=HERE.parent/"round034_block_mixture"/"block_probe.py"
    text=old.read_text();node=next(n for n in ast.parse(text).body
        if isinstance(n,ast.FunctionDef)and n.name=="field_labels")
    copied=ast.get_source_segment(text,node)
    assert copied==inspect.getsource(field_labels).rstrip()
    data={"scope":"1D deterministic true-box mean/range controls; no high-dimensional tail samples",
        "python":sys.executable,"numpy":np.__version__,"L":"7/6",
        "source_SHA256":hashlib.sha256((HERE/"source_tables.py").read_bytes()).hexdigest(),
        "validation_SHA256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "field_labels_verbatim_round034":True,"n":src.n,"K":src.K,"Q":src.Q,
        "pchi":src.pchi,"background_lower":src.background_lower,
        "Dsup_bound":src.Dsup_bound,"Dprime_bound":src.Dprime_bound,
        "tables":table_rows,"rational_peak_rows":rational_rows,
        "roundoff_interval_certified":False,"seconds":time.time()-started}
    target.write_text(json.dumps(data,separators=(",",":")))
    print(json.dumps({"seconds":data["seconds"],"Q":src.Q,"pchi":src.pchi,
        "source_SHA256":data["source_SHA256"],"tables":[{k:v for k,v in x.items()
        if k!="reference_rows"} for x in table_rows],
        "max_rational_peak_kernel_error":max(x["kernel_error"]for x in rational_rows)},indent=2))


if __name__=="__main__":main()
