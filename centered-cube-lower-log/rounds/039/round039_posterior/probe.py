"""Deterministic PB/Binomial posterior-prior checks, before general proof.

Log-domain positive DP; no Gaussian or saddlepoint approximation. The finite
rational counterexample is separately checked with Python integers/Fraction.
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


def pb_logpmf(p):
    p=np.asarray(p,dtype=float)
    out=np.full(len(p)+1,-np.inf);out[0]=0.
    for i,x in enumerate(p):
        before=out[:i+1].copy()
        with np.errstate(divide="ignore"):
            a,b=np.log1p(-x),np.log(x)
        out[:i+1]=before+a
        out[1:i+2]=np.logaddexp(out[1:i+2],before+b)
    return out


def logsumexp(v):
    return float(np.logaddexp.reduce(v))


def inspect_case(name,p,eta=.5):
    p=np.asarray(p,dtype=np.longdouble);n=len(p)
    posterior=pb_logpmf(p);prior=pb_logpmf(np.full(n,eta))
    delta=float(p.sum()-n*eta);d=delta/math.sqrt(n)
    dispersion=float(np.sum((p-eta)**2));tau=dispersion/n
    variance=float(np.sum(p*(1-p)))
    logratio=posterior-prior
    winner=int(np.argmax(logratio));maximum=float(logratio[winner])
    # Candidate to audit: C=20/3, exponent (5/7)*Delta^2/D, only in range.
    in_range=bool(np.all((p>=.3)&(p<=.7)) and .49<=eta<=.51)
    log_candidate=(math.log(20/3)+(5/7)*delta*delta/dispersion
                   if dispersion else 0.)
    if in_range:assert maximum<=log_candidate+1e-9
    segment_rows=[]
    cumulative=np.exp(np.logaddexp.accumulate(prior))
    for r in [2,4,8,16,32,64]:
        cuts=[int(np.searchsorted(cumulative,np.longdouble(b)/r,side="left"))
              for b in range(1,r)]+[n]
        lower=[0]+[v+1 for v in cuts[:-1]]
        logmax=-math.inf;empty=0
        for lo,hi in zip(lower,cuts):
            if lo>hi:empty+=1;continue
            value=logsumexp(posterior[lo:hi+1])-logsumexp(prior[lo:hi+1])
            logmax=max(logmax,value)
        assert logmax<=maximum+1e-9
        segment_rows.append({"R":r,"empty_segments":empty,"max_log_segment_ratio":logmax})
    # Explicitly inspect all upper-tail integer intervals for the false simple cap.
    post_tail=np.exp(np.logaddexp.accumulate(posterior[::-1])[::-1])
    prior_tail=np.exp(np.logaddexp.accumulate(prior[::-1])[::-1])
    simple_multiplier=2*math.exp(2*d*d) if abs(d)<10 else math.inf
    excess=post_tail-simple_multiplier*prior_tail-2/math.sqrt(n)
    tailwinner=int(np.argmax(excess))
    return {"name":name,"n":n,"eta":eta,"p_min":float(p.min()),"p_max":float(p.max()),
        "in_frozen_probability_range":in_range,"delta":delta,"d":d,"D":dispersion,
        "dispersion_per_coordinate":tau,"variance_per_coordinate":variance/n,
        "log_sup_pmf_ratio":maximum,"sup_ratio_argmax":winner,
        "log_conditional_candidate_bound":log_candidate,
        "candidate_applies":in_range,"candidate_verified_numerically":maximum<=log_candidate+1e-9,
        "pmf_mass_residual":max(abs(math.expm1(logsumexp(posterior))),abs(math.expm1(logsumexp(prior)))),
        "segment_checks":segment_rows,
        "simple_2exp2d2_plus_2_over_sqrtn_tail_countercheck":{
            "tail_lower_cut":tailwinner,"excess_over_candidate":float(excess[tailwinner]),
            "posterior_probability":float(post_tail[tailwinner]),
            "prior_probability":float(prior_tail[tailwinner]),"multiplier":simple_multiplier}}


def exact_counterexample():
    n=256;k=147
    a=F(sum(math.comb(n,j)*9**j*7**(n-j)for j in range(k,n+1)),16**n)
    b=F(sum(math.comb(n,j)for j in range(k,n+1)),2**n)
    residual=a-16*b-F(1,8)
    assert residual>0
    # d=1, and 2e^2<16 follows e<11/4, so this also refutes the
    # smaller multiplier 2 exp(2d^2) with additive 2/sqrt(n)=1/8.
    assert F(2)*F(11,4)**2<16
    return {"n":n,"p":"9/16","eta":"1/2","d":"1","interval":"[147,256]",
        "posterior_tail":str(a),"prior_tail":str(b),"tested_multiplier":"16",
        "additive":"1/8","exact_positive_residual":str(residual),
        "residual_float":float(residual),"also_refutes":"2 exp(2d^2)*prior+2/sqrt(n)"}


def main():
    target=HERE/"RESULTS.json"
    if target.exists():raise RuntimeError("Refusing repeated deterministic screen.")
    started=time.time();cases=[]
    definitions=[("n1_center",[.5]),("n1_shift",[.7]),
        ("homogeneous_shift_n64",np.full(64,.5+1/math.sqrt(64))),
        ("homogeneous_shift_n256",np.full(256,.5+1/math.sqrt(256))),
        ("homogeneous_shift_n1024",np.full(1024,.5+1/math.sqrt(1024))),
        ("homogeneous_shift_n4096",np.full(4096,.5+1/math.sqrt(4096))),
        ("balanced_heterogeneous_n256",np.tile([.3,.7],128)),
        ("balanced_heterogeneous_n1024",np.tile([.3,.7],512)),
        ("shifted_heterogeneous_n1024",np.tile([.35,.65],512)+1/math.sqrt(1024)),
        ("asymmetric_heterogeneous_eta49",np.tile([.30,.68],512)),
        ("outside_range_concentrated_n128",np.r_[np.zeros(64),np.ones(64)]),
        ("outside_range_nearly_concentrated_n128",np.tile([.001,.999],64))]
    for name,p in definitions:
        eta=.49 if name=="asymmetric_heterogeneous_eta49"else .5
        row=inspect_case(name,p,eta);cases.append(row)
        print(json.dumps({k:row[k]for k in ["name","n","d","dispersion_per_coordinate",
            "log_sup_pmf_ratio","log_conditional_candidate_bound","candidate_applies"]}),flush=True)
    exact=exact_counterexample()
    out={"scope":"deterministic exact-law finite screen; no asymptotic/C_n claim",
        "cases":cases,"exact_counterexample":exact,
        "source_SHA256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "python":sys.executable,"numpy":np.__version__,"longdouble_bits":int(np.finfo(np.longdouble).nmant+1),
        "seconds":time.time()-started}
    target.write_text(json.dumps(out,separators=(",",":")))
    print("COMPLETE",out["seconds"],"exact_counter_residual",exact["residual_float"],flush=True)


if __name__=="__main__":main()
