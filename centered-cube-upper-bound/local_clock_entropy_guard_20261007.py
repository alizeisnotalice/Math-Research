#!/usr/bin/env python3
"""New registered finite posterior entropy guard; no old max/oracle."""
from fractions import Fraction as F
import pathlib,json,hashlib,gzip,time
from functools import lru_cache
P=pathlib.Path(__file__).resolve().parent; PRE="local_clock_entropy_guard_20261007"
TERMS=100; BITS=256; SCALE=1<<BITS
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,obj):p.write_text(json.dumps(obj,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
def outward(lo,hi):
    return F((lo*SCALE).numerator//(lo*SCALE).denominator,SCALE),F(-(-((hi*SCALE).numerator)//(hi*SCALE).denominator),SCALE)
def basiclog(q):
    z=(q-1)/(q+1);z2=z*z;term=z;s=F(0)
    for j in range(TERMS):
        s+=2*term/(2*j+1);term*=z2
    tail=2*term/((2*TERMS+1)*(1-z2))
    return s,s+tail
LN2=basiclog(F(2))
@lru_cache(None)
def log_interval(q):
    if q<=0:raise ValueError("log domain")
    if q==1:return F(0),F(0)
    k=0;m=q
    while m<1:m*=2;k-=1
    while m>=2:m/=2;k+=1
    lo,hi=basiclog(m)
    if k>=0:lo+=k*LN2[0];hi+=k*LN2[1]
    else:lo+=k*LN2[1];hi+=k*LN2[0]
    return outward(lo,hi)
def addscaled(total,c,pair):
    lo,hi=pair
    if c<0:lo,hi=hi,lo
    return total[0]+c*lo,total[1]+c*hi
def val(pair):return [str(pair[0]),str(pair[1])]
checks=0
def check(ok,label):
    global checks
    checks+=1
    if not ok:raise AssertionError(label)
def main():
    start=time.monotonic();rp=P/(PRE+"_registration.json");reg=json.loads(rp.read_text())
    check(sha(pathlib.Path(__file__))==reg["script_sha256"],"frozen code")
    for e in reg["skills_read"]:
        check(sha(pathlib.Path(e["path"]))==e["sha256"],"skill file")
    allrows=[];rounds=[]
    for rnd in reg["posterior_rounds"]:
        before=checks;ws=rnd["weights"];N=len(ws);rows=[];zero_count=0;p1_count=0
        for mask in range(1,1<<N):
            ids=[i for i in range(N) if mask>>i&1];tot=sum(ws[i] for i in ids)
            ps={i:F(ws[i],tot) for i in ids}
            for name in reg["rates"]:
                if name=="all1":aa={i:F(1) for i in ids}
                elif name=="alternating":aa={i:F(i%2) for i in ids}
                else:aa={i:F(i+1,N+1) for i in ids}
                abar=sum((ps[i]*aa[i] for i in ids),F(0))
                check(sum(ps.values(),F(0))==1 and all(p>0 for p in ps.values()),"posterior")
                check(all(a>=0 for a in aa.values()),"rates")
                flow={f"logp:{j}":-ps[j]*(aa[j]-abar) for j in ids}
                jump={f"logp:{j}":F(0) for j in ids};rem={}
                jump_records=[];H=(F(0),F(0));D=(F(0),F(0))
                for j in ids:H=addscaled(H,-ps[j],log_interval(ps[j]))
                for i in ids:
                    p=ps[i];a=aa[i];rate=a*(1-p)
                    check(rate>=0,"mass-tilted rate")
                    if p==1:
                        p1_count+=1
                        check(rate==0 and all(c==0 for c in flow.values()),"p1 no deletion/domain")
                        rem[f"log1minus:{i}"]=F(0)
                        jump_records.append({"removed":i,"rate":"0","endpoint":"p=1; deletion not evaluated"})
                        continue
                    conditional={j:ps[j]/(1-p) for j in ids if j!=i}
                    check(sum(conditional.values(),F(0))==1,"post-deletion normalisation")
                    # Independent expansion of H(p_without_i)-H(p).
                    diff={j:ps[j] for j in ids}
                    for j in conditional:diff[j]-=conditional[j]
                    for j,c in diff.items():jump[f"logp:{j}"]+=rate*c
                    rem[f"log1minus:{i}"]=rate
                    jump_records.append({"removed":i,"rate":str(rate),"conditional_p":{str(j):str(pj) for j,pj in conditional.items()},
                                         "entropy_logp_difference_coefficients":{str(j):str(c) for j,c in diff.items()},
                                         "entropy_log1minus_coefficient":"1"})
                    if a==0:
                        check(rate==0,"a0 exact zero")
                        continue
                    D=addscaled(D,-rate,log_interval(1-p))
                    check(p+(1-p)*log_interval(1-p)[0]>0,"individual positive slack")
                for j in ids:check(flow[f"logp:{j}"]+jump[f"logp:{j}"]==0,"formal logp cancellation")
                check(all(rem[f"log1minus:{i}"]==aa[i]*(1-ps[i]) for i in ids),"generator remainder")
                iszero=all(aa[i]==0 or ps[i]==1 for i in ids)
                if iszero:
                    zero_count+=1;check(D==(0,0),"zero dissipation exact")
                else:check(D[0]>0,"positive dissipation")
                check(D[0]<=D[1] and D[1]<=abar,"certified D <= abar")
                row={"original_N":N,"mask":mask,"retained_labels":ids,"retained_weights":[ws[i] for i in ids],
                     "p":{str(i):str(ps[i]) for i in ids},"rate_family":name,"a":{str(i):str(aa[i]) for i in ids},
                     "abar":str(abar),"entropy_interval":val(H),"D_interval":val(D),
                     "abar_minus_D_interval":val((abar-D[1],abar-D[0])),
                     "flow_logp_coefficients":{k:str(v) for k,v in flow.items()},
                     "jump_logp_coefficients":{k:str(v) for k,v in jump.items()},
                     "remaining_log1minus_coefficients":{k:str(v) for k,v in rem.items()},
                     "jump_records":jump_records,"p1_endpoint":len(ids)==1,"exact_zero_D":iszero}
                rows.append(row)
        allrows.extend(rows);rounds.append({"N":N,"weights":ws,"nonempty_subsets":(1<<N)-1,"rows":len(rows),
                                           "checks":checks-before,"zero_dissipation_rows":zero_count,"p1_rows":p1_count})
    products=[]
    for n in reg["product_n"]:
        before=checks;a=F(1);b=F(2);d=a/(8*n);N=n**n;extent=d*(n-1);weight=F(1,N)
        check(0<=extent<=b/2,"entire product support in Qb")
        check(0<=extent<=a/2,"entire product support even in Qa")
        check(d>0 and n>=2,"distinct spacing")
        check(d*n>a*F(1,16),"spacing formula")
        check(weight*N==1,"complete equal source mass")
        # Coordinate k ranges 0..n-1: exact half-open cell index floor((d*k)/d)=k.
        cell_coordinates=[]
        for k in range(n):
            x=d*k;idx=(x/d).numerator//(x/d).denominator
            check(idx==k and d*k<=x<d*(k+1),"distinct coordinate cell")
            check(abs(x)<=b/2,"coordinate captured in same Qb")
            cell_coordinates.append({"k":k,"coordinate":str(x),"cell_index":idx})
        l=log_interval(F(n));H=(n*l[0],n*l[1])
        check(log_interval(F(N))[0]<=H[1] and H[0]<=log_interval(F(N))[1],"log N = n log n overlapping certified intervals")
        check(H[0]>0,"positive product entropy")
        products.append({"n":n,"a":str(a),"b":str(b),"d":str(d),"source_atoms_symbolic":"d*k, k in {0,...,n-1}^n",
                         "N":str(N),"each_atom_weight":str(weight),"W":"1","receiver":["0"]*n,
                         "max_coordinate":str(extent),"coordinate_cell_certificate":cell_coordinates,
                         "same_Qb_full_capture":True,"distinct_cells":True,"entropy_interval":val(H),
                         "checks":checks-before,"explicit_atom_enumeration":False,
                         "scope":"Compressed exact product geometry and entropy only; no maximum or nearmax qualification."})
    check(len(allrows)==879 and len(products)==3,"registered counts")
    ap=P/(PRE+"_posterior_rows.json.gz")
    raw=json.dumps(allrows,ensure_ascii=False,sort_keys=True,separators=(",",":")).encode()
    ap.write_bytes(gzip.compress(raw,mtime=0))
    result={"status":"PASS","checks":checks,"posterior_rows":len(allrows),"rounds":rounds,"product_sources":products,
            "posterior_rows_artifact":{"path":str(ap),"sha256":sha(ap),"uncompressed_sha256":hashlib.sha256(raw).hexdigest()},
            "script_sha256":sha(pathlib.Path(__file__)),"registration_sha256":sha(rp),
            "interval_method":{"terms":TERMS,"outward_dyadic_bits":BITS,"log_range":"[1,2) times integer power of 2",
                               "tail":"2*z^(2T+1)/((2T+1)*(1-z^2)); signed range-reduction coefficient"},
            "elapsed_seconds":time.monotonic()-start,"old_oracle_invocations":0,
            "scope":"New finite posterior generator arithmetic and compressed product source geometry. Not nearmax, receiver integral, clock path simulation, entropy time-uniform budget or geom closure."}
    outp=P/(PRE+"_results.json");dump(outp,result)
    dump(P/(PRE+"_receipt.json"),{"status":"PASS","successful_execution_count":1,"failed_execution_count":0,
                                  "checks_successful_terminal_only":checks,"script_sha256":sha(pathlib.Path(__file__)),
                                  "registration_sha256":sha(rp),"results_sha256":sha(outp),
                                  "rows_artifact_sha256":sha(ap),"elapsed_seconds":result["elapsed_seconds"],
                                  "old_oracle_invocations":0})
    print(json.dumps({"status":"PASS","checks":checks,"posterior_rows":len(allrows),"product_sources":3,"elapsed_seconds":result["elapsed_seconds"]}))
if __name__=="__main__":main()

