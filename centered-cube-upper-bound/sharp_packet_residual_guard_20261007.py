#!/usr/bin/env python3
"""New scalar sharp-packet arithmetic guard, no source or process simulation."""
from fractions import Fraction as F
from functools import lru_cache
import pathlib,json,hashlib,time
P=pathlib.Path(__file__).resolve().parent;PRE="sharp_packet_residual_guard_20261007"
T=100;BITS=256;SCALE=1<<BITS
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,obj):p.write_text(json.dumps(obj,sort_keys=True,indent=2)+"\n")
def basic(q):
    z=(q-1)/(q+1);z2=z*z;term=z;lo=F(0)
    for j in range(T):lo+=2*term/(2*j+1);term*=z2
    return lo,lo+2*term/((2*T+1)*(1-z2))
LN2=basic(F(2))
@lru_cache(None)
def ln(q):
    if q<=0:raise ValueError("log0 prohibited")
    if q==1:return F(0),F(0)
    k=0;m=q
    while m<1:m*=2;k-=1
    while m>=2:m/=2;k+=1
    l,h=basic(m)
    if k>=0:l+=k*LN2[0];h+=k*LN2[1]
    else:l+=k*LN2[1];h+=k*LN2[0]
    l*=SCALE;h*=SCALE
    return F(l.numerator//l.denominator,SCALE),F(-(-h.numerator//h.denominator),SCALE)
def clean(c):return {x:y for x,y in c.items() if x!=1 and y!=0}
def interval(c,r):
    l=h=r
    for arg,coef in c.items():
        a,b=ln(arg)
        if coef<0:a,b=b,a
        l+=coef*a;h+=coef*b
    return l,h
def serialize(c):return {str(k):str(v) for k,v in sorted(c.items())}
def fmt(ab):return [str(z) for z in ab]
checks=0
def ck(v,s):
    global checks
    checks+=1
    if not v:raise AssertionError(s)
def kernel(t):return t/(1+t)**2
def main():
    start=time.monotonic();regp=P/(PRE+"_registration.json");reg=json.loads(regp.read_text())
    ck(sha(pathlib.Path(__file__))==reg["script_sha256"],"frozen code")
    for e in reg["skills_read"]:ck(sha(pathlib.Path(e["path"]))==e["sha256"],"read skill")
    rows=[];rounds=[];integrals=[]
    for rd in reg["rounds"]:
        before=checks;rrows=[]
        for vs in rd["v"]:
            v=F(vs);Gc=clean({1+v:F(1)});Gr=-v/(1+v);Gi=interval(Gc,Gr)
            ck(v>0 and Gi[0]>0,"positive G")
            # Exact differentiation of log(1+v)-v/(1+v).
            derivative=F(1)/(1+v)-F(1)/(1+v)**2
            ck(derivative==v/(1+v)**2,"G derivative rational identity")
            ck(kernel(v)==derivative,"primitive derivative")
            # Compressed finite Darboux enclosure for integral_0^v kernel(t)dt.
            nodes={F(0),v}
            for j in range(1,5):nodes.add(v*j/4)
            power=F(1,4)
            while power<v:nodes.add(power);power*=2
            if v>1:nodes.add(F(1))
            nodes=sorted(nodes);lo=hi=F(0);pieces=[]
            for a,b in zip(nodes,nodes[1:]):
                vals=[kernel(a),kernel(b)]
                if a<=1<=b:vals.append(F(1,4))
                mn=min(vals);mx=max(vals);l=(b-a)*mn;h=(b-a)*mx
                ck(0<=l<=h,"Darboux piece")
                lo+=l;hi+=h
                pieces.append({"left":str(a),"right":str(b),"kernel_min":str(mn),"kernel_max":str(mx),"integral_lower":str(l),"integral_upper":str(h)})
            ck(lo<Gi[0] and Gi[1]<hi,"certified G inside independent integral enclosure")
            ck(F(1)/(1+v)-F(1)==Gr,"primitive constant equals G rational term")
            integrals.append({"v":str(v),"G_interval":fmt(Gi),"derivative":str(derivative),
                              "integral_interval":fmt((lo,hi)),"Darboux_pieces":pieces,
                              "mixed_kernel_identity":"G(v)=integral_0^v t/(1+t)^2 dt = v^2 integral_0^1 s/(1+vs)^2 ds",
                              "identity_scope":"Analytic primitive and exact substitution, finite enclosure only as consistency guard."})
            for uname,u in [("min",max(F(1),v)),("saturation",1+v),("double",2*(1+v))]:
                pmax=v/u;pc=1-F(1)/u
                for pname,p in [("zero",F(0)),("halfmax",pmax/2),("max",pmax),("switch",min(pmax,pc))]:
                    ck(u>=max(1,v) and 0<=p<=pmax<=1,"feasible")
                    ck(pc>=0,"switch domain")
                    if p==1:
                        branch="p1; min is log u";hc=clean({u:F(1)});hr=-p
                        ck(u>=v and pmax==1,"p1 qualification")
                    elif p<=pc:
                        branch="minus log(1-p)";hc=clean({1/(1-p):F(1)});hr=-p
                        ck(u*(1-p)>=1,"exact branch order")
                    else:
                        branch="log u";hc=clean({u:F(1)});hr=-p
                        ck(u*(1-p)<1,"exact branch order")
                    gap=dict(Gc)
                    for q,c in hc.items():gap[q]=gap.get(q,F(0))-c
                    gap=clean(gap);rr=Gr-hr;gapint=interval(gap,rr);hi_=interval(hc,hr)
                    exactgapzero=not gap and rr==0
                    if exactgapzero:
                        ck(gapint==(0,0),"symbolic saturation equality")
                    else:ck(gapint[0]>0,"strict sharp upper bound")
                    if not hc and hr==0:
                        ck(hi_==(0,0),"exact h zero");sign="zero"
                    elif hi_[1]<0:sign="negative";ck(p>0,"negative h retained")
                    else:ck(hi_[0]>0,"positive h interval");sign="positive"
                    if uname=="saturation" and p==pmax:
                        ck(exactgapzero,"registered saturating optimizer")
                    row={"round":rd["name"],"v":str(v),"u_case":uname,"u":str(u),"p_case":pname,"p":str(p),
                         "pmax":str(pmax),"switch_p":str(pc),"branch":branch,"h_log_coefficients":serialize(hc),
                         "h_rational_term":str(hr),"h_interval":fmt(hi_),"h_sign":sign,"G_interval":fmt(Gi),
                         "G_minus_h_log_coefficients":serialize(gap),"G_minus_h_rational_term":str(rr),
                         "G_minus_h_interval":fmt(gapint),"symbolic_saturation":exactgapzero,
                         "p1_log0_evaluated":False}
                    rows.append(row);rrows.append(row)
        rounds.append({"name":rd["name"],"v":rd["v"],"rows":len(rrows),"checks":checks-before,
                       "sign_counts":{s:sum(x["h_sign"]==s for x in rrows) for s in ("zero","negative","positive")},
                       "saturation_rows":sum(x["symbolic_saturation"] for x in rrows),"p1_rows":sum(x["p"]=="1" for x in rrows)})
    zeros=[]
    for u in (F(1),F(2)):
        v=p=F(0)
        ck(v/u==p and interval({},0)==(0,0),"v0 h=G=0")
        ck(v/(1+v)**2==0,"v0 derivative")
        zeros.append({"v":"0","u":str(u),"p":"0","G":"0","h":"0","derivative":"0","integral":"0"})
    ck(len(rows)==108 and len(integrals)==9,"registered counts")
    result={"status":"PASS","checks":checks,"rows":rows,"rounds":rounds,"zero_endpoints":zeros,"integrals":integrals,
            "script_sha256":sha(pathlib.Path(__file__)),"registration_sha256":sha(regp),"elapsed_seconds":time.monotonic()-start,
            "log_method":{"terms":T,"outward_dyadic_bits":BITS,"tail":"2*z^201/[201(1-z^2)]","range":"q=2^k*m, 1<=m<2",
                          "equality":"Exact identical log arguments merged and rational constants cancelled before enclosure."},
            "scope":"Finite scalar envelope, derivative algebra and specified positive integral enclosure only. No source, winner, history, actual law, global numerical optimization, nearmax or geom proof."}
    op=P/(PRE+"_results.json");dump(op,result)
    dump(P/(PRE+"_receipt.json"),{"status":"PASS","successful_execution_count":1,"failed_execution_count":0,
                                  "checks_successful_terminal_only":checks,"script_sha256":sha(pathlib.Path(__file__)),
                                  "registration_sha256":sha(regp),"results_sha256":sha(op),
                                  "elapsed_seconds":result["elapsed_seconds"],"old_oracle_invocations":0})
    print(json.dumps({"status":"PASS","checks":checks,"rows":len(rows),"rounds":rounds,"elapsed_seconds":result["elapsed_seconds"]}))
if __name__=="__main__":main()

