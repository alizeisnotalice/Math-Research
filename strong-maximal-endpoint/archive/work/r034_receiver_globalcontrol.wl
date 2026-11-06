(* r034 necessary control: uniform receiver-c refinement at fixed source H. *)
root=DirectoryName[$InputFileName];args=$ScriptCommandLine;
repo=FileNameJoin[{DirectoryName[root],"Math-Research-strong-maximal","strong-maximal-endpoint"}];
seed=If[Length[args]>=2,args[[2]],FileNameJoin[{repo,"rounds","r030","receiver_seed.json"}]];
out=If[Length[args]>=3,args[[3]],FileNameJoin[{root,"r034_receiver_globalcontrol.json"}]];
$r032LibraryMode=True;Get[FileNameJoin[{repo,"code","receiver_continuous32.wl"}]];$r032LibraryMode=False;
srcScale=1.0276924765700282; global0=c+.10;
run[offset_?NumericQ]:=Module[{z=evaluateReceiverC[ConstantArray[global0+offset,ng],srcScale,128,32]},<|"offset_from_r033_global"->offset,"receiver_c"->global0+offset,"result"->z|>];
coarse=run /@ {.005,.01,.0125,.025};best=First@Ordering[coarse[[All,"result","min_average"]],-1];
(* One local 1D refinement beyond the four required controls. *)
refOffset=coarse[[best,"offset_from_r033_global"]]+If[coarse[[best,"offset_from_r033_global"]]>=.025,.005,.0025];refine=run[refOffset];
all=Join[coarse,{refine}];winner=all[[First@Ordering[all[[All,"result","min_average"]],-1]]];
answer=<|"status"->"r034 fixed-source uniform receiver-c control; floating numerical evidence only","source_scale"->srcScale,"r033_global_receiver_c"->global0,
 "coarse_N128_soft32"->coarse,"one_refinement_N128_soft32"->refine,"winner"->winner,
 "local_coordinate_best_min_reference"->1.002940553706471,"interpretation"->"Uniform c control is required before attributing any gain to a coordinate-specific perturbation. No FOSD claim."|>;
Export[out,answer,"RawJSON"];Print[<|"winner_offset"->winner["offset_from_r033_global"],"winner_min"->winner["result","min_average"]|>];
