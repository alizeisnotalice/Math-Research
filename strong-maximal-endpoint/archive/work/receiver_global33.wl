(* r033: global receiver-c calibration and one asymmetric two-block check.
   Uses the r032 split integrator as a library.  Floating evidence only. *)
root33 = DirectoryName[$InputFileName]; args33 = $ScriptCommandLine;
seed33 = If[Length[args33] >= 2, args33[[2]], FileNameJoin[{root33, "variable29_archived_theta0.json"}]];
out33 = If[Length[args33] >= 3, args33[[3]], FileNameJoin[{root33, "receiver_global33.json"}]];
$r032LibraryMode = True;
Get[FileNameJoin[{root33, "receiver_continuous32.wl"}]];
$r032LibraryMode = False;

(* No FOSD is assumed: all choices below are evaluated directly. *)
mode33 = If[Length[args33] >= 4, args33[[4]], "full"];
runGlobal[cj_?NumericQ, no_Integer] := evaluateReceiverC[ConstantArray[cj, ng], baseScale, no, 32];
If[mode33 === "asymonly",
 bestC = c + .10; trialScale = 1.0276924765700282; scaleDrop = baseScale-trialScale,
 coarseC = c + {.05, .10, .20};
 coarse = (<|"receiver_c" -> #, "offset" -> (# - c), "result" -> runGlobal[#, 128]|> &) /@ coarseC;
 pickMax[rr_] := First@Ordering[rr[[All, "result", "min_average"]], -1];
 bestCoarse = coarse[[pickMax[coarse]]];
 refC = DeleteDuplicates[Select[bestCoarse["receiver_c"] + {-.025, .025}, # > 0. &]];
 refine = (<|"receiver_c" -> #, "offset" -> (# - c), "result" -> runGlobal[#, 128]|> &) /@ refC;
 allRuns = Join[coarse, refine]; best = allRuns[[pickMax[allRuns]]]; bestC = best["receiver_c"];
 excess = best["result", "min_average"] - 1.00001;
 scaleDrop = Min[.02, Max[0., .8 excess/.04]];
 trialScale = baseScale (1. - scaleDrop);
 globalTrial256 = evaluateReceiverC[ConstantArray[bestC, ng], trialScale, 256, 32];
];

(* Coordinate parameters: u[[j]] uses (c1,c2)=(cA,cB), while its mirror
   u[[j+ng-1]] uses (c1,c2)=(cB,cA).  Pair symmetry is thereby preserved. *)
probVecBlocks[t1_?NumericQ, t2_?NumericQ, ca_List, cb_List] := Module[{ans = ConstantArray[0., ng], z, q1, q2},
 ans[[1]] = polyC[t1, ca[[1]]] polyC[t2, ca[[1]]] Exp[-t1 - t2];
 Do[z = u[[j]]; q1 = polyC[t1 - z[[1]], ca[[j]]]; q2 = polyC[t2 - z[[2]], cb[[j]]];
    ans[[j]] += q1 q2 Exp[z[[1]] + z[[2]] - t1 - t2];
    z = u[[j + ng - 1]]; q1 = polyC[t1 - z[[1]], cb[[j]]]; q2 = polyC[t2 - z[[2]], ca[[j]]];
    ans[[j]] += q1 q2 Exp[z[[1]] + z[[2]] - t1 - t2], {j, 2, ng}]; ans];
constantBlocks[t1_?NumericQ, left_?NumericQ, right_?NumericQ, h_?NumericQ, ca_List, cb_List] := Module[{cost = 0., av = ConstantArray[0., ng], z, fac},
 If[right <= left, Return[{cost, av}]]; cost = Exp[tau - t1] (Exp[-left] - Exp[-right]) phi[h];
 z = {0., 0.}; fac = polyC[t1, ca[[1]]] Exp[-t1]; av[[1]] = fac (gbcdfC[right, ca[[1]]] - gbcdfC[left, ca[[1]]]);
 Do[z = u[[j]]; fac = polyC[t1-z[[1]], ca[[j]]] Exp[z[[1]]-t1]; av[[j]] += fac (gbcdfC[right-z[[2]], cb[[j]]] - gbcdfC[left-z[[2]], cb[[j]]]);
    z = u[[j+ng-1]]; fac = polyC[t1-z[[1]], cb[[j]]] Exp[z[[1]]-t1]; av[[j]] += fac (gbcdfC[right-z[[2]], ca[[j]]] - gbcdfC[left-z[[2]], ca[[j]]]), {j,2,ng}]; {cost,av}];
softBlocks[t1_?NumericQ, left_?NumericQ, right_?NumericQ, amp_, ca_List, cb_List] := Module[{x,wt,cost=0.,av=ConstantArray[0.,ng],h},
 If[right<=left,Return[{cost,av}]];{x,wt}=quad[32,left,right];
 Do[h=profile[likeAmp[x[[i]],amp]];cost+=wt[[i]] Exp[tau-t1-x[[i]]] phi[h];av+=wt[[i]] h probVecBlocks[t1,x[[i]],ca,cb],{i,Length[x]}];{cost,av}];
runBlocks[ca_List, cb_List, srcScale_?NumericQ, no_Integer] := Module[{ww=srcScale w0,tx,tw,cost=0.,av=ConstantArray[0.,ng],amp,roots,z,a},
 {tx,tw}=outerQuad[no]; Do[amp=sourceAmp[tx[[i]],ww];roots=boundary[#,amp]&/@{1.,N[n],phiPrime};
  {z,a}=constantBlocks[tx[[i]],roots[[1]],roots[[2]],1.,ca,cb];cost+=tw[[i]]z;av+=tw[[i]]a;
  {z,a}=constantBlocks[tx[[i]],roots[[3]],hi,cap,ca,cb];cost+=tw[[i]]z;av+=tw[[i]]cap a;
  {z,a}=softBlocks[tx[[i]],roots[[2]],roots[[3]],amp,ca,cb];cost+=tw[[i]]z;av+=tw[[i]]a,{i,Length[tx]}];
 <|"source_scale"->srcScale,"scaled_cost"->cost,"actual_logJ"->Log[cost]-tau,"averages_by_group"->av/coeff,"min_average"->Min[av/coeff]|>];

If[mode33 === "asymonly",
 bindGroup = 24; ca = ConstantArray[bestC, ng]; cb = ConstantArray[bestC, ng];
 globalBlock128 = runBlocks[ca, cb, trialScale, 128];
 ca[[bindGroup]] += .05;
 asymmetric128 = runBlocks[ca, cb, trialScale, 128];
 correction = <|"status" -> "r033 corrected mirror-swapped coordinate-parameter check; floating numerical evidence only",
  "selected_global_receiver_c" -> bestC, "source_scale" -> trialScale, "outer_nodes" -> 128, "soft_inner_order" -> 32,
  "group" -> bindGroup, "mapping" -> "positive u[j]: (c1,c2)=(cA,cB); mirror u[j+32]: (c1,c2)=(cB,cA)",
  "point_positive" -> u[[bindGroup]], "point_mirror" -> u[[bindGroup+ng-1]],
  "global_ca_equals_cb_baseline" -> globalBlock128,
  "asymmetric_c1_raised" -> <|"c1" -> ca[[bindGroup]], "c2" -> cb[[bindGroup]], "result" -> asymmetric128|>,
  "interpretation" -> "This replaces the earlier invalid two-block mapping. Compare same-source, same-128-node direct quadratures only."|>;
 Export[out33, correction, "RawJSON"]; Print[<|"baseline_min"->globalBlock128["min_average"],"asym_min"->asymmetric128["min_average"]|>]; Exit[]];

bestAv = globalTrial256["averages_by_group"]; bindGroup = First@Ordering[bestAv, 1];
ca = ConstantArray[bestC, ng]; cb = ConstantArray[bestC, ng]; ca[[bindGroup]] += .05;
asymmetric128 = runBlocks[ca, cb, trialScale, 128];

result33 = <|"status" -> "r033 global receiver-c calibration with one asymmetric block check; floating numerical evidence only",
 "source_c" -> c, "base_source_scale" -> baseScale, "coarse_global_c_N128_soft32" -> coarse,
 "refinement_global_c_N128_soft32" -> refine, "selected_global_c" -> bestC,
 "selected_global_offset" -> bestC-c, "selected_at_base_scale_N128_soft32" -> best["result"],
 "one_source_scale_trial" -> <|"scale" -> trialScale,"scale_drop"->scaleDrop,"global_N256_soft32"->globalTrial256|>,
 "asymmetric_block_check" -> <|"group"->bindGroup,"block_A_c"->ca[[bindGroup]],"block_B_c"->cb[[bindGroup]],
   "point_A"->u[[bindGroup]],"point_B"->u[[bindGroup+ng-1]],"N128_soft32"->asymmetric128|>,
 "interpretation" -> "No first-order stochastic-order claim in c is used. Compare direct quadratures only; a global receiver c may dominate a local or asymmetric perturbation."|>;
Export[out33, result33, "RawJSON"];
Print[<|"best_c"->bestC,"trial_min"->globalTrial256["min_average"],"asym_min"->asymmetric128["min_average"]|>];
