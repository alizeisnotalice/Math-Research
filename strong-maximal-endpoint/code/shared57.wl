(* Exact three-fiber counterexample to pairwise common-action certification. *)
Clear[r];
q1[r_] := Piecewise[{{1, 0 < r <= 2}, {2/r, r > 2}}];
q2[r_] := Piecewise[{{0, 0 < r <= 1}, {2 - 2/r, 1 < r <= 2}, {1, 2 < r <= 3}, {3/r, r > 3}}];
q3[r_] := Piecewise[{{1, 0 < r <= 1}, {1/r, 1 < r <= 2}, {2 - 3/r, 2 < r <= 3}, {3/r, r > 3}}];
qs = {q1[r], q2[r], q3[r]};
ms = MaxValue[{#, r > 0}, r] & /@ qs;
common = MaxValue[{Mean[qs], r > 0}, r];
pairLoss = Table[MinValue[{(ms[[i]] - qs[[i]]) + (ms[[j]] - qs[[j]]), r > 0}, r], {i, 3}, {j, 3}];
tripleLoss = MinValue[{Total[ms - qs], r > 0}, r];
Print["Wolfram version: ", $Version];
Print["fiber maxima = ", ms];
Print["common value and regret = ", {common, Mean[ms] - common}];
Print["pair loss matrix = ", pairLoss];
Print["triple loss = ", tripleLoss];
Print["Phi2 cost = ", 8 + 4 Log[2]];
Clear[m, del, i, j];
prefix[m_, del_, i_, j_] := If[j < i, 1, If[j == i, 1 - del/j, 1]];
Print["M=5 prefix table (delta=1/2): ", Table[prefix[5, 1/2, i, j], {i, 1, 5}, {j, 1, 5}]];
Print["all proper-subset common maxima check M=5: ",
 And @@ Table[And @@ Table[prefix[5, 1/2, i, j] == 1, {i, Delete[Range[5], j]}], {j, 1, 5}]];
Print["aggregate shared value M=5, delta=1/2: ", 1 - (1/2)/5^2];
Print["Phi2 cost M=5, delta=1/2: ", FullSimplify[5^2 - 1/2 + (5 - 1) (1 + 1/2) Log[1 + 1/2]]];

phi2[x_]:=If[x<=1,x,x(1+Log[x])];
gm3={{1,1,0},{0,2,1},{1,0,2}};
checks=<|"three_fiber_common_regret"->TrueQ[common==8/9&&Mean[ms]-common==1/9],"pair_losses_zero"->TrueQ[pairLoss==ConstantArray[0,{3,3}]],"triple_loss"->TrueQ[tripleLoss==1/3],"actual_three_slab_Phi2_cost"->TrueQ[FullSimplify[Total[Map[phi2,gm3,{2}],2]-(8+4Log[2])]===0]|>;
Do[gm=Table[1-If[i==j,dd,0]+If[i<mm&&j==i+1,dd,0],{i,mm},{j,mm}];
 pr=(Accumulate[#]/Range[mm])&/@gm;expect=Table[If[i==j,1-dd/j,1],{i,mm},{j,mm}];
 checks["prefix_integrals_M"<>ToString[mm]<>"_delta"<>ToString[dd]]=TrueQ[pr==expect];
 checks["common_regret_M"<>ToString[mm]<>"_delta"<>ToString[dd]]=TrueQ[Max[Mean[pr]]==1-dd/mm^2];
 checks["actual_cost_M"<>ToString[mm]<>"_delta"<>ToString[dd]]=TrueQ[FullSimplify[Total[Map[phi2,gm,{2}],2]-(mm^2-dd+(mm-1)(1+dd)Log[1+dd])]===0],{mm,{2,3,5,8,12}},{dd,{1/2,1/3}}];
If[!And@@Values[checks],Print[checks];Exit[2]];
out=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/shared57.json"];
Export[out,<|"wolfram"->$Version,"checks"->checks,"all_pass"->True,"scope"->"Exact finite input integrals/costs; general fixed-arity proof is analytic"|>,"RawJSON"];
