"""Continuous-uniform 1-Lipschitz MGF/convex-order search; bounded deterministic work.

Piecewise affine laws are integrated directly, not sampled. Randomness only
chooses a fixed seeded finite list of rational functions, not tail observations.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import math
import sys
import time
import numpy as np

HERE=Path(__file__).resolve().parent
SEED=202610290411


def log_sinhc(x):
    x=abs(float(x))
    if x==0:return 0.
    if x<1e-4:return math.log1p(x*x/6+x**4/120+x**6/5040)
    return math.log(math.sinh(x)/x)


def stoploss(knots,values,threshold):
    total=F(0)
    for a,z,u,v in zip(knots,knots[1:],values,values[1:]):
        lo,hi=min(u,v),max(u,v)
        if lo>=threshold:total+=(z-a)*((u+v)/2-threshold)
        elif hi>threshold:total+=(z-a)*(hi-threshold)**2/(2*(hi-lo))
    return total


def uniform_stoploss(t):
    if t<=-F(1,2):return -t
    if t>=F(1,2):return F(0)
    return (F(1,2)-t)**2/2


def inspect(name,knots,values):
    assert knots[0]==0 and knots[-1]==1
    assert all(abs(v-u)<=z-a for a,z,u,v in zip(knots,knots[1:],values,values[1:]))
    mean=sum(((z-a)*(u+v)/2 for a,z,u,v in zip(knots,knots[1:],values,values[1:])),F(0))
    centered=[v-mean for v in values]
    variance=sum(((z-a)*(u*u+u*v+v*v)/3 for a,z,u,v in zip(knots,knots[1:],centered,centered[1:])),F(0))
    assert variance<=F(1,12)
    max_stoploss_excess=F(0)
    for t in [F(j,40)for j in range(-24,25)]:
        excess=stoploss(knots,centered,t)-uniform_stoploss(t)
        assert excess<=0
        max_stoploss_excess=max(max_stoploss_excess,excess)
    rows=[]
    for t in [-64.,-16.,-4.,-1.,-.1,0.,.1,1.,4.,16.,64.]:
        terms=[math.log(float(z-a))+t*float((u+v)/2)+log_sinhc(t*float(v-u)/2)
               for a,z,u,v in zip(knots,knots[1:],centered,centered[1:])]
        logmgf=float(np.logaddexp.reduce(terms))
        loguniform=log_sinhc(t/2)
        assert logmgf<=loguniform+2e-13
        assert loguniform<=t*t/24+2e-13
        rows.append({"t":t,"log_MGF":logmgf,"log_uniform_MGF":loguniform,
                     "proposed_quadratic_bound":t*t/24})
    return {"name":name,"knots":list(map(str,knots)),"values":list(map(str,values)),
            "mean":str(mean),"variance":str(variance),"stoploss_49_rational_thresholds_exact":True,
            "MGF_rows":rows}


def main():
    target=HERE/"uniform_results.json"
    if target.exists():raise RuntimeError("Refusing replacement.")
    started=time.time()
    definitions=[("affine",[F(0),F(1)],[F(0),F(1)]),
        ("reverse_affine",[F(0),F(1)],[F(1),F(0)]),
        ("tent",[F(0),F(1,2),F(1)],[F(0),F(1,2),F(0)]),
        ("absolute_value",[F(0),F(1,2),F(1)],[F(1,2),F(0),F(1,2)]),
        ("clipped",[F(0),F(1,3),F(1)],[F(0),F(1,3),F(1,3)]),
        ("asymmetric_positive_tail",[F(0),F(3,4),F(1)],[F(0),F(0),F(1,4)]),
        ("asymmetric_zigzag",list(map(F,[0,F(1,10),F(3,10),F(3,5),1])),
         list(map(F,[0,F(1,10),-F(1,10),F(1,5),-F(1,10)]))),
        ("eight_oscillations",[F(i,16)for i in range(17)],[F(i%2,16)for i in range(17)])]
    rng=np.random.default_rng(SEED)
    for j in range(64):
        knots=[F(i,32)for i in range(33)];values=[F(0)]
        for slope in rng.integers(-8,9,size=32):values.append(values[-1]+F(int(slope),8*32))
        definitions.append((f"seeded_piecewise_{j}",knots,values))
    cases=[inspect(*x)for x in definitions]
    # This is a counterexample to using only oscillation, or to using an
    # affine function under a two-endpoint law instead of continuous uniform.
    counter={"law":"centered Bernoulli at -1/2,+1/2", "stoploss_at_zero":"1/4",
        "uniform_stoploss_at_zero":"1/8","positive_excess_exact":"1/8",
        "log_MGF_t1":math.log(math.cosh(.5)),"log_uniform_MGF_t1":log_sinhc(.5)}
    assert counter["log_MGF_t1"]>counter["log_uniform_MGF_t1"]
    out={"scope":"finite piecewise-affine uniform-law search; no source/tail MC or CI",
        "seed_for_function_list":SEED,"cases":cases,"out_of_scope_counterexample":counter,
        "script_SHA256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "python":sys.executable,"numpy":np.__version__,"seconds":time.time()-started}
    target.write_text(json.dumps(out,separators=(",",":")))
    print(json.dumps({"functions":len(cases),"exact_stoploss_checks":len(cases)*49,
        "MGF_checks":len(cases)*11,"counterexample":counter,"seconds":out["seconds"]},indent=2))


if __name__=="__main__":main()
