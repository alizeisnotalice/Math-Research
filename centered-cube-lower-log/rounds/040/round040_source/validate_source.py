"""Deterministic R4 compatibility and R4/8/16/32 true-box controls; zero MC."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
import sys
import time
import numpy as np
from source_tables import Source

HERE=Path(__file__).resolve().parent


def exact_geometry(src):
    half=F(1,2*src.w*src.v)
    minimum=F(1)
    mass_tests=0;left=0;right=0
    for b in range(src.R):
        atoms=[]
        for j,num in enumerate(src.atom_offset_numerators[b]):
            u=F(int(num),2*src.R)
            lo,hi=u-half,u+half
            assert 0<=lo<=hi<=1
            # Actual fraction of the uniform source box staying in its cell.
            assert (min(F(1),hi)-max(F(0),lo))/(2*half)==1
            left+=lo==0;right+=hi==1;mass_tests+=1
            atoms.append((F(j)+u)/src.K)
        for segment in [atoms[:src.cut],atoms[src.cut:]]:
            gaps=[z-a for a,z in zip(segment,segment[1:])]+[segment[0]+1-segment[-1]]
            minimum=min(minimum,*gaps)
    assert minimum>=F(1,4*src.K)
    assert src.delta+src.binwidth<float(minimum)
    return {"box_mass_in_cell_exactly_one_tests":mass_tests,
            "box_touch_left_cell_edge_count":left,"box_touch_right_cell_edge_count":right,
            "minimum_class_internal_circular_gap":str(minimum),
            "partial_candidate_at_most_one_per_endpoint":True,
            "source_mass_exact":"1","peak_class_mass_exact":[str(F(src.cut,src.K)),str(F(src.K-src.cut,src.K))]}


def reference_events(src,L,b):
    ell=L-1
    atoms=src.atoms[b]
    events=[]
    for shift in [-1,0,1,2]:
        for side in [-1,1]:
            for ds in [-1,1]:
                events.extend((atoms+shift-.5+side*ell/2+ds*src.delta/2).tolist())
    for plus in [True,False]:
        events.extend(x-.5 for x in src.background_breakpoints(ell,plus))
    return np.asarray(events)


def reference_bins(src,L,b):
    ell=L-1;queries={0,src.Q//7,src.Q-1}
    for j in [0,src.cut-1,src.cut,src.K-1]:
        for side in [-1,1]:
            x=(src.atoms[b,j]+side*ell/2-.5)%1
            queries.add(int(math.floor(x*src.Q)))
    for side in [-1,1]:
        queries.add(int(math.floor(((src.a+side*ell/2-.5)%1)*src.Q)))
    return sorted(queries)


def exact_peak_bin_mean(src,b,q,plus):
    ell=F(1,6);delta=F(1,src.w*src.K*src.v)
    left,right=F(q,src.Q),F(q+1,src.Q)
    base=[(F(j)+F(int(src.atom_offset_numerators[b,j]),2*src.R))/src.K
          for j in range(src.K)if (j<src.cut)==plus]
    atoms=[a+s for s in [-1,0,1,2]for a in base]
    cuts={left,right}
    for atom in atoms:
        for side in [-1,1]:
            for ds in [-1,1]:
                x=atom-F(1,2)+side*ell/2+ds*delta/2
                if left<x<right:cuts.add(x)
    cuts=sorted(cuts);clip=lambda x:max(F(0),min(F(1),x))
    def height(x):
        y=x+F(1,2)
        return sum((clip((y-a+ell/2+delta/2)/delta)
                    -clip((y-a-ell/2+delta/2)/delta)for a in atoms),F(0))
    heights=[height(x)for x in cuts]
    return sum(((z-a)*(u+v)/2 for a,z,u,v in zip(cuts,cuts[1:],heights,heights[1:])),F(0))/(right-left)


def compatibility():
    sys.dont_write_bytecode=True
    oldpath=HERE.parent/"round038_highdim/source_tables.py"
    spec=importlib.util.spec_from_file_location("source38_readonly",oldpath)
    oldmod=importlib.util.module_from_spec(spec);spec.loader.exec_module(oldmod)
    old,new=oldmod.Source(),Source(4)
    assert np.array_equal(old.atoms,new.atoms)
    assert np.array_equal(old.labels,new.labels)
    tests=0
    for length in [1.,10/9,7/6,11/9]:
        for b in range(4):
            a,z=old.evaluate(length,b),new.evaluate(length,b)
            assert a.keys()==z.keys()
            for key in a:
                if isinstance(a[key],np.ndarray):assert np.array_equal(a[key],z[key]),key
                else:assert a[key]==z[key],key
            tests+=1
    return {"old_source_SHA256":hashlib.sha256(oldpath.read_bytes()).hexdigest(),
            "R4_atoms_labels_bit_identical":True,"evaluate_full_Q_bit_identical_tests":tests,
            "lengths":["1","10/9","7/6","11/9"]}


def main():
    target=HERE/"SOURCE_VALIDATION.json"
    if target.exists():raise RuntimeError("Refusing replacement.")
    started=time.time();compat=compatibility();rows=[];L=7/6
    nodes,weights=np.polynomial.legendre.leggauss(24)
    for R in [4,8,16,32]:
        src=Source(R);geometry=exact_geometry(src)
        max_mean=0.;max_range_ratio=0.;max_mass=0.;max_rel=0.;observations=0
        fraction_rows=[]
        for b in range(R):
            d=src.evaluate(L,b);events=reference_events(src,L,b)
            max_mass=max(max_mass,abs(d["normalizer"]-1),
                         abs(d["mean_plus"].mean()-src.pchi),
                         abs(d["mean_minus"].mean()-(1-src.pchi)))
            max_rel=max(max_rel,float(d["relative_variance_class_max"].max()))
            for q in reference_bins(src,L,b):
                left,right=q/src.Q,(q+1)/src.Q
                cuts=np.unique(np.r_[left,events[(left<events)&(events<right)],right])
                x=np.concatenate([(a+z)/2+(z-a)/2*nodes for a,z in zip(cuts,cuts[1:])])
                w=np.concatenate([(z-a)*weights/2 for a,z in zip(cuts,cuts[1:])]);w/=w.sum()
                for plus in [True,False]:
                    name="plus"if plus else"minus"
                    actual=src.point_class_kernel(x,L,b,plus)
                    mean=float(w@actual);variance=float(w@((actual-mean)**2))
                    error=abs(mean-d["mean_"+name][q]);cap=d["range_"+name][q]
                    observed=np.ptp(src.point_class_kernel(np.unique(np.r_[x,cuts]),L,b,plus))
                    assert error<d["mean_error_abs"]+2e-13
                    assert observed<=cap+2e-13 and variance<=cap**2/4+1e-20
                    max_mean=max(max_mean,error);max_range_ratio=max(max_range_ratio,float(observed/cap))
                    observations+=1
            if b in {0,R//2,R-1}:
                for plus in [True,False]:
                    j=0 if plus else src.cut
                    q=int(math.floor(((src.atoms[b,j]+(L-1)/2-.5)%1)*src.Q))
                    exact=exact_peak_bin_mean(src,b,q,plus)
                    nominal=src._peak_bin_mean_and_transition_count(
                        np.array([src.bins_center[q]+.5]),L-1,b,plus)[0][0]
                    error=abs(float(exact)-nominal)*src.p/(src.K*L)
                    assert error<1e-12
                    fraction_rows.append({"b":b,"bin":q,"class":"plus"if plus else"minus",
                        "exact_extra_peak_mean":str(exact),"kernel_absolute_error":error})
        record={"R":R,"geometry":geometry,"full_Q_evaluate_tables":R,
            "reference_class_bin_checks":observations,"normalizer_and_class_mass_max_residual":max_mass,
            "max_reference_mean_error":max_mean,"max_observed_range_over_cap":max_range_ratio,
            "max_relative_variance_class_bound":max_rel,
            "analytic_mean_error_abs":d["analytic_mean_error_abs"],
            "floating_point_guard_abs":d["floating_point_guard_abs"],
            "pchi":src.pchi,"background_lower":src.background_lower,
            "fraction_peak_checks":fraction_rows}
        rows.append(record);print(json.dumps({k:v for k,v in record.items()if k not in ["geometry","fraction_peak_checks"]}),flush=True)
    output={"scope":"deterministic 1D source controls; no MC or full feature cache",
        "L":"7/6","R4_compatibility":compat,"R_checks":rows,
        "source_SHA256":hashlib.sha256((HERE/"source_tables.py").read_bytes()).hexdigest(),
        "validation_SHA256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "python":sys.executable,"numpy":np.__version__,"ordinary_FP_not_interval_certified":True,
        "seconds":time.time()-started}
    target.write_text(json.dumps(output,indent=2)+"\n")
    print("COMPLETE",output["seconds"],output["source_SHA256"],flush=True)


if __name__=="__main__":main()
