(* Exact analytic finite-Pareto cone volume for r036 accepted receiver anchors. *)
root=DirectoryName[$InputFileName];repo=FileNameJoin[{DirectoryName[root],"Math-Research-strong-maximal","strong-maximal-endpoint"}];seed=Import[FileNameJoin[{repo,"rounds","r030","receiver_seed.json"}],"RawJSON"];run=Import[FileNameJoin[{root,"r036_geometry_curve.json"}],"RawJSON"];out=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],FileNameJoin[{root,"r036_geometry_volume.json"}]];
v=N[seed["v_nonnegative"]];d=N[seed["receiver_d"]];n=seed["n"];c0=seed["c"];cr=c0+.1075;
uof[dd_]:=Join[{{0.,0.}},Transpose[{Rest[dd]+Rest[v],Rest[dd]-Rest[v]}],Transpose[{Rest[dd]-Rest[v],Rest[dd]+Rest[v]}]];
coneV[uu_]:=Module[{s,xx,yy,grp,best,tot=0.,i},s=SortBy[uu,First];xx=DeleteDuplicates[s[[All,1]]];
 Do[grp=Select[s,#[[1]]==xx[[i]]&];best=Min[Select[s,#[[1]]<=xx[[i]]&][[All,2]]];
  tot+=(Exp[-xx[[i]]]-If[i<Length[xx],Exp[-xx[[i+1]]],0.])Exp[-best],{i,Length[xx]}];tot];
dOut=d+.15(v/16.)^2;uBase=uof[d];uOut=uof[dOut];vb=coneV[uBase];vo=coneV[uOut];jb=run["baseline_N128_soft32","result","actual_logJ"];jo=run["best_trial_N256_soft64","result","actual_logJ"];
(* The receiver c is common to both rows; any c-dependent physical prefactor
   cancels in their relative score.  Record the candidate B-ratio separately,
   not as an asserted absolute normalization. *)
logBShift=n(Log[1+Exp[-cr]]-Log[1+Exp[-c0]]);
res=<|"status"->"finite Pareto cone volume for r036 anchors; analytic formula, floating anchor locations","all_65_groups_strict_in_stored_runs"->True,"baseline_volume"->vb,"outer_up_volume"->vo,
 "baseline_raw_log_score_logV_minus_logJ"->Log[vb]-jb,"outer_up_raw_log_score_logV_minus_logJ"->Log[vo]-jo,"raw_score_difference_outer_minus_base"->(Log[vo]-jo)-(Log[vb]-jb),
 "receiver_c"->cr,"source_c"->c0,"candidate_common_c_log_B_ratio_times_n"->logBShift,"normalization_status"->"A receiver-c physical prefactor is not established here; it cancels in the reported same-c relative difference and raw scores are not compared to scalar Q.",
 "r028_scope_note"->"fixed global receiver c with full receive is covered by the r028 obstruction, so no scalar-Q improvement is claimed."|>;
Export[out,res,"RawJSON"];Print[<|"Vbase"->vb,"Vout"->vo,"rawDelta"->res["raw_score_difference_outer_minus_base"]|>];
