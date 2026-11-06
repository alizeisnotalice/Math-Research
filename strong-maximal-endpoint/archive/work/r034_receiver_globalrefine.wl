(* Single r034 global-c refinement at the bracket midpoint +0.0075. *)
root=DirectoryName[$InputFileName];args=$ScriptCommandLine;repo=FileNameJoin[{DirectoryName[root],"Math-Research-strong-maximal","strong-maximal-endpoint"}];
seed=If[Length[args]>=2,args[[2]],FileNameJoin[{repo,"rounds","r030","receiver_seed.json"}]];out=If[Length[args]>=3,args[[3]],FileNameJoin[{root,"r034_receiver_globalrefine.json"}]];
$r032LibraryMode=True;Get[FileNameJoin[{repo,"code","receiver_continuous32.wl"}]];$r032LibraryMode=False;
srcScale=1.0276924765700282;cReceiver=c+.1075;res=evaluateReceiverC[ConstantArray[cReceiver,ng],srcScale,128,32];
Export[out,<|"status"->"single r034 global-c refinement; floating evidence only","source_scale"->srcScale,"offset_from_r033_global"->.0075,"receiver_c"->cReceiver,"N128_soft32"->res|>,"RawJSON"];Print[<|"min"->res["min_average"],"group"->First@Ordering[res["averages_by_group"],1]|>];
