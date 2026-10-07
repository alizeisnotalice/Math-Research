#!/usr/bin/env python3
"""Registered exact mesh guard. Does not run the saved continuous oracle."""
from fractions import Fraction as F
import pathlib, json, hashlib, gzip, time
P=pathlib.Path(__file__).resolve().parent
PREFIX="high_m_packet_mesh_guard_20261007"
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,obj): p.write_text(json.dumps(obj,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
def artifact(name,obj):
    p=P/name
    raw=json.dumps(obj,ensure_ascii=False,sort_keys=True,separators=(",",":")).encode()
    p.write_bytes(gzip.compress(raw,mtime=0))
    return {"path":str(p),"sha256":sha(p),"uncompressed_sha256":hashlib.sha256(raw).hexdigest(),"bytes":len(p.read_bytes())}
checks=0
def check(ok,label):
    global checks
    checks+=1
    if not ok: raise AssertionError(label)
def mesh(n,a,b,xi):
    delta=a*xi/(2*n)
    w=(b-a)/delta
    J=-(-w.numerator//w.denominator)
    nodes=[min(b,a+j*delta) for j in range(J+1)]
    return delta,J,nodes
def mesh_checks(n,a,b,xi,delta,J,nodes):
    check(n>0 and 0<xi<=1 and 0<a<b,"parameters")
    check(delta==a*xi/(2*n),"delta")
    w=(b-a)/delta
    check(J>=w and J-1<w,"ceiling")
    check(len(nodes)==J+1 and nodes[0]==a and nodes[-1]==b,"endpoints")
    ratios=[]
    for j,s in enumerate(nodes):
        check(a<=s<=b,"window")
        check(s==min(b,a+j*delta),"node construction")
        if j:
            r=s/nodes[j-1]; ratios.append(r)
            check(s>nodes[j-1],"strict monotonicity")
            check(s-nodes[j-1]<=delta,"adjacent step")
            check(r<=1+xi/(2*n),"adjacent ratio")
            check(r**n<=1+xi,"exact power envelope")
    return ratios
def main():
    start=time.monotonic()
    regp=P/(PREFIX+"_registration.json")
    reg=json.loads(regp.read_text())
    for entry in reg["frozen_inputs"]:
        check(sha(pathlib.Path(entry["path"]))==entry["sha256"],"frozen input "+entry["path"])
    summaries=[]; synthetic_artifacts=[]; total_nodes=0; total_adjacent=0
    for n in reg["synthetic_dimensions"]:
        windows=[("unit",F(1),F(2)),("half",F(1,2),F(3,2)),("narrow",F(1),F(1)+F(1,16*n))]
        for label,a,b in windows:
            for xs in reg["xi"]:
                before=checks; xi=F(xs); delta,J,nodes=mesh(n,a,b,xi)
                ratios=mesh_checks(n,a,b,xi,delta,J,nodes)
                obj={"n":n,"window":label,"a":str(a),"b":str(b),"xi":xs,"delta":str(delta),"J":J,
                     "nodes":[str(s) for s in nodes],"adjacent_ratios":[str(s) for s in ratios],
                     "all_exact_power_checks_pass":True}
                name=f"{PREFIX}_mesh_n{n}_{label}_xi{xi.numerator}_{xi.denominator}.json.gz"
                saved=artifact(name,obj); synthetic_artifacts.append(saved)
                maxr=max(ratios)
                summaries.append({k:obj[k] for k in ("n","window","a","b","xi","delta","J")}|
                                 {"nodes":len(nodes),"adjacent":J,"checks":checks-before,
                                  "max_adjacent_ratio":str(maxr),"max_ratio_power":str(maxr**n),"artifact":saved})
                total_nodes+=len(nodes); total_adjacent+=J
    oldreg=json.loads((P/reg["original_registration"]).read_text())
    oldresults=json.loads((P/reg["original_results"]).read_text())
    profilemap={pathlib.Path(e["path"]).name:e for e in oldresults["profiles"]}
    selection_artifacts=[]; selected_summary=[]; total_queries=0
    for sel in reg["saved_selection"]:
        before=checks
        # Independently verify first nonempty selection in frozen manifest/record order.
        found=None
        for pi,plan in enumerate(oldreg["plans"]):
            if plan["round"]!=sel["round"]: continue
            filename=next(fn for fn in profilemap if plan["case_id"] in fn)
            obj=json.loads(gzip.decompress((P/filename).read_bytes()))
            for ri,record in enumerate(obj["records"]):
                if not record["empty_response"] and F(record["M"])>0:
                    found=(pi,ri,filename,obj,record); break
            if found is not None: break
        pi,ri,filename,obj,record=found
        check(pi==sel["plan_index"] and ri==sel["record_index"] and filename==sel["profile"],"first nonempty selection")
        case=obj["case"]; plan=oldreg["plans"][pi]
        check(case==plan,"complete source manifest equality")
        check(case["case_id"]==sel["case_id"] and record["M"]==sel["saved_M"] and record["R"]==sel["saved_R"],"frozen receiver")
        source=case["source"]; n=case["n"]; a=F(case["a"]); b=F(case["b"])
        x=list(map(F,record["x"])); atoms=[list(map(F,y)) for y in source["atoms"]]
        weights=list(map(F,source["weights"])); W=sum(weights,F(0))
        check(len(atoms)==case["N"]==len(weights) and W==1 and all(w>0 for w in weights),"complete labelled source mass")
        D=[2*max(abs(xi-yi) for xi,yi in zip(x,y)) for y in atoms]
        check([str(d) for d in D]==record["D"],"distance equality with saved profile")
        M=F(record["M"]); runs=[]
        for xs in reg["xi"]:
            xi=F(xs); delta,J,nodes=mesh(n,a,b,xi)
            mesh_checks(n,a,b,xi,delta,J,nodes)
            queries=[]; last_mass=F(0)
            for s in nodes:
                labels=[i for i,d in enumerate(D) if d<=s]
                mass=sum((weights[i] for i in labels),F(0)); response=mass/s**n
                check(0<=mass<=W and mass>=last_mass,"full cumulative mass")
                check(response<=M,"new finite response below frozen continuous M")
                check(all((D[i]<=s)==(i in labels) for i in range(len(D))),"closed capture all labels")
                check(response*s**n==mass,"response identity")
                queries.append({"s":str(s),"captured_labels":labels,"boundary_labels":[i for i,d in enumerate(D) if d==s],
                                "cumulative_mass":str(mass),"response":str(response)})
                last_mass=mass
            meshM=max(F(q["response"]) for q in queries)
            meshR=next(q["s"] for q in queries if F(q["response"])==meshM)
            check(meshM<=M<=(1+xi)*meshM,"saved continuous mesh envelope")
            runs.append({"xi":xs,"delta":str(delta),"J":J,"meshM":str(meshM),"smallest_mesh_winner":meshR,
                         "saved_continuous_M":str(M),"saved_continuous_R":record["R"],"queries":queries})
            total_queries+=len(queries)
        saved=artifact(PREFIX+f"_saved_round{sel['round']}.json.gz",
                       {"selection":sel,"complete_case":case,"receiver":record["x"],"kind":record["kind"],
                        "exact_D":[str(d) for d in D],"saved_continuous_M":record["M"],"saved_continuous_R":record["R"],
                        "source_mass":str(W),"runs":runs})
        selection_artifacts.append(saved)
        selected_summary.append({"round":sel["round"],"case_id":sel["case_id"],"n":n,"N":len(atoms),"checks":checks-before,
                                 "saved_M":record["M"],"saved_R":record["R"],
                                 "comparisons":[{"xi":r["xi"],"J":r["J"],"meshM":r["meshM"],"meshR":r["smallest_mesh_winner"]} for r in runs],
                                 "artifact":saved})
    check(total_nodes==27252 and total_adjacent==27225,"registered synthetic counts")
    check(total_queries==597,"registered new fixed query count")
    result={"status":"PASS","checks":checks,"synthetic_meshes":len(summaries),"synthetic_nodes":total_nodes,
            "synthetic_adjacent_checks":total_adjacent,"new_saved_source_queries":total_queries,"saved_receiver_mesh_comparisons":9,
            "summaries":summaries,"selected_summaries":selected_summary,"elapsed_seconds":time.monotonic()-start,
            "script_sha256":sha(pathlib.Path(__file__)),"registration_sha256":sha(regp),
            "continuous_oracle_invocations":0,"old_MC_invocations":0,
            "scope":"Exact finite mesh arithmetic and new fixed-query cumulative masses compared with frozen saved continuous M. No numerical continuous maximization, near-identity limit, actual FIRST or Lebesgue integration claim."}
    resultp=P/(PREFIX+"_results.json"); dump(resultp,result)
    receipt={"status":"PASS","successful_execution_count":1,"failed_execution_count":0,"checks_successful_terminal_only":checks,
             "results_sha256":sha(resultp),"script_sha256":sha(pathlib.Path(__file__)),"registration_sha256":sha(regp),
             "elapsed_seconds":result["elapsed_seconds"],"artifacts":synthetic_artifacts+selection_artifacts,
             "continuous_oracle_invocations":0,"old_MC_invocations":0}
    dump(P/(PREFIX+"_receipt.json"),receipt)
    print(json.dumps({k:result[k] for k in ("status","checks","synthetic_meshes","synthetic_nodes","synthetic_adjacent_checks","new_saved_source_queries","elapsed_seconds")}))
if __name__=="__main__": main()

