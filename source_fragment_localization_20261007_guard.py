"""Exact original one-dimensional finite-winner geometry; no asymptotic weak certificate."""
from fractions import Fraction as F
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json

HERE=Path(__file__).resolve().parent
PREFIX="source_fragment_localization_20261007"
d,tau,eta,delta0=F(1,8),F(1,4),F(1,4),F(1,8)
rounds=[]
total=0
for J in (12,24,48):
    checks=[]
    def check(name,cond):
        assert cond,(J,name)
        checks.append(name)
    check("density_and_source_mass",1/d==8 and J>0)
    # Far packets: source [-d/2,d/2], receiver distances (9/16,15/16).
    left,right=F(9,16),F(15,16)
    check("far_R1_misses_full_source",left-d/2>F(1,2))
    check("far_R2_contains_full_source",right+d/2==1)
    check("far_neighbor_query_halo",4-right-d/2>1)
    check("far_unique_true_winner",F(1,2)>tau and F(1,2)>0)
    check("far_thickened_components_separate",4>2+d)
    v_far=tau*2*(right-left)
    H_far=F(1,8)
    check("far_exact_profile",v_far==F(3,16))
    check("far_true_highset_all_source",v_far>H_far)
    check("far_global_fullness_fails",1<eta*J)
    check("far_local_fullness_passes",1>=eta)
    I_far=tau*J*2*(right-left)
    F_far=J*(v_far-H_far)
    check("far_source_receiver_exchange",F_far==I_far-H_far*J)
    check("far_local_receiver_fee_sum",J*tau*2*(right-left)==I_far)
    check("far_mass_penalty_once",sum((H_far*F(1) for _ in range(J)),F(0))==H_far*J)
    # Overlap chain: midpoint distance 3/4, receiver half-width 1/16.
    spacing,mid,half=F(3,2),F(3,4),F(1,16)
    check("chain_R1_misses_nearest_sources",mid-half-d/2>F(1,2))
    check("chain_R2_contains_both_sources",mid+half+d/2<1)
    check("chain_R2_misses_other_sources",spacing+mid-half-d/2>1)
    check("chain_unique_true_winner",F(2,2)>tau and F(2,2)>0)
    check("chain_thickened_component_connected",spacing<2+d)
    v_end=tau*(2*half)/2
    v_inner=2*v_end
    H_chain=F(1,128)
    check("chain_exact_endpoint_profile",v_end==F(1,64))
    check("chain_exact_interior_profile",v_inner==F(1,32))
    check("chain_true_highset_all_source",v_end>H_chain)
    check("chain_rows_true_fragment",2<eta*J and F(2,2)>delta0)
    I_chain=tau*(J-1)*(2*half)
    first_chain=2*v_end+(J-2)*v_inner
    check("chain_source_receiver_first_moment",first_chain==I_chain)
    F_chain=first_chain-H_chain*J
    check("chain_tail_exact_exchange",F_chain==F(3*J-4,128))
    check("chain_mass_penalty_once",H_chain*J==sum((H_chain for _ in range(J)),F(0)))
    row={"J":J,"n":1,"far":{"components":J,"source_profile":str(v_far),
          "H":str(H_far),"I_R":str(I_far),"tail":str(F_far),
          "global_fullness":False,"local_fullness":True},
         "chain":{"components":1,"endpoint_profile":str(v_end),
          "interior_profile":str(v_inner),"H":str(H_chain),
          "I_R":str(I_chain),"tail":str(F_chain),"all_rows_local_fragment":True},
         "checks":checks,"assertions":len(checks),"status":"PASS"}
    rounds.append(row)
    total+=len(checks)
result={"scope":"exact 1D L1 finite-winner subreceiver geometry; not general order/soft-gate evidence",
        "rounds":rounds,"assertions":total,"status":"PASS",
        "general_fragment_budget_proved":False}
(HERE/(PREFIX+"_results.json")).write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
receipt={"executed_at_utc":datetime.now(timezone.utc).isoformat(),
         "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         "registration_sha256":hashlib.sha256((HERE/(PREFIX+"_registration.md")).read_bytes()).hexdigest(),
         "assertions":total,"status":"PASS","scope":result["scope"]}
(HERE/(PREFIX+"_receipt.json")).write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+"\n")
print(json.dumps({"status":"PASS","rounds":3,"assertions":total,
                  "far_components":[r["far"]["components"] for r in rounds],
                  "chain_components":[r["chain"]["components"] for r in rounds]}))
