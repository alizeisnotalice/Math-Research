"""Exact small true-box positive-polynomial controls, then numerical log MGF/tail.

No high-dimensional MC, no normal approximation, no frozen-source modification.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
import sys
import time
import numpy as np

HERE=Path(__file__).resolve().parent


def quantiles(n,r):
    counts=[math.comb(n,s)for s in range(n+1)]
    cdf=list(itertools.accumulate(counts))
    upper=[next(s for s,c in enumerate(cdf)if r*c>=b*2**n)for b in range(1,r)]+[n]
    lower=[0]+[x+1 for x in upper[:-1]]
    prior=[sum((F(counts[s],2**n)for s in range(a,z+1)),F(0))if a<=z else F(0)
           for a,z in zip(lower,upper)]
    return lower,upper,prior


def class_kernel(x,length,b,plus):
    k=4;delta=F(1,128);ell=length-1
    if length==1:return F(1,2)
    y=x+F(1,2);ca,cz=(F(0),F(1,2))if plus else(F(1,2),F(1))
    overlap=sum((max(F(0),min(y+ell/2,cz+s)-max(y-ell/2,ca+s))for s in [-1,0,1,2]),F(0))
    clip=lambda a:max(F(0),min(F(1),a));peak=F(0)
    for j,label in enumerate([0,1,2,3]):
        if(j<2)!=plus:continue
        phase=(F(b)+F(1,2))/4 if plus else(F(1,2)-b)/4
        atom=(F(j)+(F(label,4)+phase)%1)/k
        for shift in [-1,0,1,2]:
            a=atom+shift
            peak+=clip((y-a+ell/2+delta/2)/delta)-clip((y-a-ell/2+delta/2)/delta)
    return (F(1,2)+F(1,2)*overlap+F(1,8)*peak)/length


def segment_value(kernels,lo,hi):
    if lo>hi:return F(0)
    coeff=[F(1)]
    for minus,plus in kernels:
        nxt=[F(0)]*(len(coeff)+1)
        for s,v in enumerate(coeff):nxt[s]+=v*minus;nxt[s+1]+=v*plus
        coeff=nxt
    return sum(coeff[lo:hi+1],F(0))


def logmeanexp(x):
    return float(np.logaddexp.reduce(x)-math.log(len(x)))


def component_checks(values,ratios,rlist,q0,upper_mean):
    if not any(values):
        return {"zero_component":True,"S_bound":0.,"mean":"0"}
    assert all(v>0 for v in values)
    h=np.log(np.array(values,dtype=float));mean=sum(values)/len(values)
    assert upper_mean>=mean
    c=np.log(np.array(ratios,dtype=float));s=float(c@c)
    r=np.array(rlist,dtype=float);q0=float(q0)
    constant=4*(1+float(max(rlist,default=F(0))))/(1-q0)**2
    surrogate=constant*float(np.log1p(r).sum())
    assert s<=surrogate+1e-13
    mgf=[]
    for t in [-16.,-4.,-1.,.5,1.,4.,16.]:
        lhs=logmeanexp(t*(h-h.mean()));rhs=t*t*s/8
        assert lhs<=rhs+1e-12
        mgf.append({"t":t,"log_MGF":lhs,"Hoeffding_bound":rhs})
    tails=[]
    for multiplier in [F(1),F(101,100),F(11,10),F(3,2),F(2)]:
        threshold=upper_mean*multiplier
        actual=F(sum(v>threshold for v in values),len(values))
        if multiplier<=1:bound=1.
        elif s==0:bound=0.
        else:bound=math.exp(-2*math.log(float(multiplier))**2/s)
        assert float(actual)<=bound+1e-13
        tails.append({"threshold_over_mean_upper":str(multiplier),"exact_probability":str(actual),"McD_bound":bound})
    return {"zero_component":False,"mean":str(mean),"mean_upper":str(upper_mean),
        "S_bound":s,"C_log_product_surrogate":surrogate,"log_MGF_checks":mgf,"tail_checks":tails}


def small_case(n,length):
    r=4;delta=F(1,128);width=delta/6
    lower,upper,pi=quantiles(n,r);centers=[]
    for i in range(n):
        b,j=i%r,i%4;label=j
        phase=(F(b)+F(1,2))/r if j<2 else(F(1,2)-b)/r
        atom=(F(j)+(F(label,4)+phase)%1)/4
        centers.append((atom+(1 if i%2==0 else -1)*(length-1)/2-F(1,2))%1)
    # Each coordinate independently takes both bin endpoints. This tests the
    # arbitrary-independent-law McD claim, not the uniform-bin 1/12 lemma.
    z=[[[[(class_kernel(c+side*width/2,length,b,False),
           class_kernel(c+side*width/2,length,b,True))][0]
          for b in range(r)]for side in [-1,1]]for c in centers]
    hidden=list(itertools.product([0,1],repeat=n))
    components=[[segment_value([z[i][o[i]][b]for i in range(n)],lower[b],upper[b])
                 for o in hidden]for b in range(r)]
    fixed=[[math.prod(sum(z[i][o[i]][b])for i in range(n))for o in hidden]for b in range(r)]
    ratio_by_table=[];r_by_table=[];q_by_table=[]
    for b in range(r):
        ratios=[];rs=[];qs=[]
        for i in range(n):
            candidates=[]
            for s in range(2):
                a=z[i][0][b][s];v=z[i][1][b][s]
                mean=(a+v)/2;osc=abs(a-v)
                q=osc/mean;assert q<1
                candidates.append(q)
            q=max(candidates);qs.append(q);rs.append(q*q/4);ratios.append(1/(1-q))
            # Before logarithms, check positive multilinear ratio exactly.
            if lower[b]<=upper[b]:
                for idx,o in enumerate(hidden):
                    changed=list(o);changed[i]=1-changed[i];jdx=hidden.index(tuple(changed))
                    a,v=components[b][idx],components[b][jdx]
                    exact_kernel_ratio=max(max(z[i][0][b][s],z[i][1][b][s])/
                                           min(z[i][0][b][s],z[i][1][b][s])for s in range(2))
                    assert a/v<=exact_kernel_ratio<=ratios[-1]
        ratio_by_table.append(ratios);r_by_table.append(rs);q_by_table.append(max(qs,default=F(0)))
    models={}
    for name,terms in [("joint",components),("actual_prior",[[pi[b]*v for v in fixed[b]]for b in range(r)])]:
        means=[sum(x)/len(x)for x in terms]
        # Explicit nonnormalized upper weights, as in prior-probability enclosures.
        uppers=[m*F(101,100)if m else F(0)for m in means]
        records=[component_checks(vals,ratio_by_table[b],r_by_table[b],q_by_table[b],uppers[b])
                 for b,vals in enumerate(terms)]
        total_upper=sum(uppers);sumchecks=[]
        for multiplier in [F(1),F(101,100),F(11,10),F(3,2),F(2)]:
            threshold=total_upper*multiplier
            target=[sum(terms[b][i]for b in range(r))>threshold for i in range(len(hidden))]
            union=[]
            for i in range(len(hidden)):
                events=[terms[b][i]>threshold*uppers[b]/total_upper for b in range(r)if uppers[b]]
                union.append(any(events));assert not target[i]or union[i]
            if multiplier<=1:bound=1.
            else:
                gap=math.log(float(multiplier))
                bound=min(1.,sum(math.exp(-2*gap*gap/x["S_bound"])for x in records if x["S_bound"]>0))
            actual=F(sum(target),len(hidden));assert float(actual)<=bound+1e-13
            sumchecks.append({"threshold_over_sum_mean_upper":str(multiplier),
                "threshold_allocation_covers_pointwise_exactly":True,"exact_tail":str(actual),"sum_McD_bound":bound})
        models[name]={"component_checks":records,"sum_tail_checks":sumchecks}
    return {"n":n,"K":4,"R":4,"L":str(length),"delta":str(delta),"query_bin_width":str(width),
        "query_law":"independent equiprobable bin endpoints", "hidden_states":len(hidden),
        "prior":[str(x)for x in pi],"empty_segments":[b for b in range(r)if lower[b]>upper[b]],
        "max_r":str(max(x for row in r_by_table for x in row)),
        "max_q":str(max(q_by_table)),"exact_positive_polynomial_coordinate_ratios_pass":True,
        "models":models}


def source_budget():
    B=F(29879,44945280);eta_floor=F(49,100);error=F(1,10**8)
    # tilde_mu-e >= mu_true-2e; L<=2 => L*(tilde_mu-e)>=eta_floor-4e.
    mean_floor_scaled=eta_floor-4*error
    rbound=B*B/(4*mean_floor_scaled**2)
    assert rbound<F(1,10**6)
    h=F(1,48960);analytic=F(1,510)*F(880,7)*h*h/12
    assert analytic+F(1,10**13)<error
    return {"B":str(B),"assumed_total_mean_error_budget":str(error),
        "scaled_mean_lower_bound":str(mean_floor_scaled),"quarter_cap_r_upper":str(rbound),
        "quarter_cap_r_upper_float":float(rbound),"r0":"1/1000000","q0":"1/500",
        "C_at_r0":"1000001/249001", # 4(1+1e-6)/(.998)^2
        "analytic_midpoint_error_upper":str(analytic),"certificate_exact":True,
        "FP_qualification":"budget covers analytic error plus nominal guard; complete FP error not certified"}


def main():
    target=HERE/"RESULTS.json"
    if target.exists():raise RuntimeError("Refusing replacement.")
    started=time.time()
    cases=[small_case(n,L)for n in [1,3,5]for L in [F(1),F(7,6)]]
    budget=source_budget()
    out={"scope":"small exact true-box polynomial controls plus numerical log-MGF; zero highdim MC",
        "cases":cases,"source_r0_budget":budget,
        "source_SHA256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "python":sys.executable,"numpy":np.__version__,"seconds":time.time()-started}
    target.write_text(json.dumps(out,separators=(",",":")))
    print(json.dumps({"case_count":len(cases),"hidden_states_total":sum(c["hidden_states"]for c in cases),
        "r0_certificate":budget,"seconds":out["seconds"]},indent=2))


if __name__=="__main__":main()
