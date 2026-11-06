(* r041 refinement: independent transpose and quadrature checks.  This loads
   only the library definitions, including all oriented anchor averages;
   the base scan is disabled and its JSON is never written. *)
rootR41=DirectoryName[$InputFileName];argsR41=$ScriptCommandLine;repoR41=DirectoryName[rootR41];outR41=If[Length[argsR41]>=2,argsR41[[2]],FileNameJoin[{$TemporaryDirectory,"unbalanced_refine41.json"}]];seedR41=If[Length[argsR41]>=3,argsR41[[3]],FileNameJoin[{repoR41,"rounds","r030","receiver_seed.json"}]];
$unbalanced41LibraryMode=True;Block[{$ScriptCommandLine={"unbalanced_blocks41.wl",FileNameJoin[{$TemporaryDirectory,"unused41.json"}],seedR41}},Get[FileNameJoin[{rootR41,"unbalanced_blocks41.wl"}]]];$unbalanced41LibraryMode=False;
(* First coordinate block 16/48, then its t1<->t2 transpose. *)
a164832=run41["n16_48_N256_soft32",16,48,0.,256,32];
a481632=run41["n48_16_N256_soft32",48,16,0.,256,32];
a164864=run41["n16_48_N256_soft64",16,48,0.,256,64];
a481664=run41["n48_16_N256_soft64",48,16,0.,256,64];
mirrorR41=Join[{1},Range[ng+1,2ng-1],Range[2,ng]];
transposeCheckR41[x_,y_]:=<|
 "cost_abs_difference"->Abs[x["scaled_cost"]-y["scaled_cost"]],
 "actual_logJ_abs_difference"->Abs[x["actual_logJ"]-y["actual_logJ"]],
 "cone_volume_abs_difference"->Abs[x["accepted_cone_volume"]-y["accepted_cone_volume"]],
 "score_abs_difference"->Abs[x["accepted_cone_log_score"]-y["accepted_cone_log_score"]],
 "max_mirror_per_anchor_average_difference"->Max[Abs[x["anchor_averages"]-y["anchor_averages"][[mirrorR41]]]],
 "accepted_indices_mirror_match"->(x["accepted_anchor_indices"]===Sort[mirrorR41[[y["accepted_anchor_indices"]]]])|>;
quadCheckR41[x_,y_]:=<|
 "cost_abs_difference"->Abs[x["scaled_cost"]-y["scaled_cost"]],
 "actual_logJ_abs_difference"->Abs[x["actual_logJ"]-y["actual_logJ"]],
 "score_abs_difference"->Abs[x["accepted_cone_log_score"]-y["accepted_cone_log_score"]],
 "max_per_anchor_average_difference"->Max[Abs[x["anchor_averages"]-y["anchor_averages"]]],
 "accepted_indices_match"->(x["accepted_anchor_indices"]===y["accepted_anchor_indices"] )|>;
resR41=<|"status"->"r041 256-node unequal-block quadrature and transpose validation; floating only",
 "wolfram"->$Version,"parameters"-><|"blocks"->{{16,48},{48,16}},"q"->0.,"outer_nodes"->256,"soft_orders"->{32,64}|>,
 "n16_48"-><|"soft32"->a164832,"soft64"->a164864|>,"n48_16"-><|"soft32"->a481632,"soft64"->a481664|>,
 "transpose_checks"-><|"soft32"->transposeCheckR41[a164832,a481632],"soft64"->transposeCheckR41[a164864,a481664]|>,
 "soft_order_checks"-><|"n16_48"->quadCheckR41[a164832,a164864],"n48_16"->quadCheckR41[a481632,a481664]|>,
 "note"->"The 65 anchor averages are explicit; comparing 16/48 to 48/16 after mirror permutation checks t1,t2 exchange, not a group-average symmetry assumption."|>;
Export[outR41,resR41,"RawJSON"];Print[<|"output"->outR41,"transpose"->resR41["transpose_checks"],"quadrature"->resR41["soft_order_checks"]|>];
