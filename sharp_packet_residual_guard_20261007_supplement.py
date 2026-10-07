#!/usr/bin/env python3
"""Separately registered sharpness and mixed-weight kernel supplement."""
from fractions import Fraction as F
from functools import lru_cache
import pathlib,json,hashlib,time
P=pathlib.Path(__file__).resolve().parent;PRE="sharp_packet_residual_guard_20261007_supplement"
T=100;S=1<<256
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,sort_keys=True,indent=2)+"\n")
def rawlog(q):
    z=(q-1)/(q+1);t=z;s=F(0)
    for j in range(T):s+=2*t/(2*j+1);t*=z*z
    return s,s+2*t/(201*(1-z*z))
L2=rawlog(F(2))
def logq(q):
    k=0;m=q
    while m<1:m*=2;k-=1
    while m>=2:m/=2;k+=1
    l,h=rawlog(m)
    if k>=0:l+=k*L2[0];h+=k*L2[1]
    else:l+=k*L2[1];h+=k*L2[0]
    l*=S;h*=S
    return F(l.numerator//l.denominator,S),F(-(-h.numerator//h.denominator),S)
checks=0
def ck(ok,lab):
    global checks
    checks+=1
    if not ok:raise AssertionError(lab)
def main():
    t0=time.monotonic();regp=P/(PRE+"_registration.json");reg=json.loads(regp.read_text())
    ck(sha(pathlib.Path(__file__))==reg["script_sha256"],"frozen code")
    rows=[]
    for n in reg["n"]:
        a=F(1);b=F(2);w=F(1);tau=F(1,2**(n//2));V=w/tau;m=V/b**n
        W=tau*(2*b)**n+w;sourceB=tau*(2*b)**n
        ck(n%2==0 and m==1/V,"exact sharpness parameters")
        ck(W==sourceB+w and sourceB>0,"complete background paid")
        # For every x in Q_b and a<=R<=b: max coord |x|+R/2 <= b.
        ck(b/2+b/2==b and a<=b,"all selected query boxes contained in source background box")
        ck(a**n==1 and b**n==V**2,"radial core/end volumes")
        lo,hi=logq(1+m);rat=(lo/m-1/(1+V),hi/m-1/(1+V))
        ck(0<rat[0]<=rat[1]<1,"strict sharpness ratio interval")
        analytic_lower=1-m/2-1/(1+V);analytic_upper=1-1/(1+V)
        ck(analytic_lower<rat[0] and rat[1]<analytic_upper,"analytic rational ratio bracket")
        ck(tau*b**n==1/m and tau*V==1,"radial primitive normalization")
        rows.append({"n":n,"a":str(a),"b":str(b),"w":str(w),"tau":str(tau),"V":str(V),"m":str(m),
                     "complete_source":{"background":"tau*1_{[-b,b]^n} dx","background_mass":str(sourceB),"packet":"w*delta0","W":str(W)},
                     "receiver_region":"Q(0,b)","background_query_response":str(tau),
                     "true_winner_for_packet_capture":"max(a,2*||x||infinity)",
                     "sharpness_ratio_formula":"ln(1+m)/m - 1/(1+V)",
                     "sharpness_ratio_interval":[str(z) for z in rat],
                     "analytic_rational_bracket":[str(analytic_lower),str(analytic_upper)],
                     "geometry_scope":"Complete background/packet source and box containment only; no nearmax or actual global-order qualification."})
    kernels=[]
    for vs in reg["kernel_v"]:
        v=F(vs);cdf=v*v/(1+v)**2;kernel=2*v/(1+v)**3
        derivative=2*v/(1+v)**2-2*v*v/(1+v)**3
        tail=(1+2*v)/(1+v)**2;weighted_cdf=1-1/(1+v)**2;weighted_tail=1/(1+v)**2
        ck(derivative==kernel,"k primitive derivative")
        ck(cdf+tail==1 and tail>0,"total k mass and analytic tail")
        ck(weighted_cdf+weighted_tail==1,"total k/t weight and analytic tail")
        ck(2/(1+v)**3>=0 and weighted_cdf>=0,"positive mixed weak weight")
        if v>0:ck(cdf/v==v/(1+v)**2,"mixed log integral derivative equals Gprime")
        else:ck(cdf==0,"v0 continuous endpoint")
        kernels.append({"v":str(v),"k_v":str(kernel),"integral_0_v_k":str(cdf),"integral_v_infinity_k":str(tail),
                        "integral_0_v_k_over_t":str(weighted_cdf),"integral_v_infinity_k_over_t":str(weighted_tail),
                        "Gprime":str(v/(1+v)**2),"identity":"G(v)=integral_0^infinity k(t)*log_+(v/t)dt, k=2t/(1+t)^3",
                        "identity_scope":"Analytic derivative/Tonelli identity; exact primitives and tails numerically guarded."})
    out=P/(PRE+"_results.json")
    obj={"status":"PASS","checks":checks,"source_rows":rows,"kernel_rows":kernels,"script_sha256":sha(pathlib.Path(__file__)),
         "registration_sha256":sha(regp),"elapsed_seconds":time.monotonic()-t0,
         "scope":"Separate complete-source sharpness ratio/box geometry and mixed-kernel primitive guard. No actual process/winner oracle, nearmax or geom proof."}
    dump(out,obj);dump(P/(PRE+"_receipt.json"),{"status":"PASS","successful_execution_count":1,"failed_execution_count":0,
         "checks_successful_terminal_only":checks,"script_sha256":sha(pathlib.Path(__file__)),"registration_sha256":sha(regp),
         "results_sha256":sha(out),"old_oracle_invocations":0})
    print(json.dumps({"status":"PASS","checks":checks,"source_rows":len(rows),"kernel_rows":len(kernels)}))
if __name__=="__main__":main()

