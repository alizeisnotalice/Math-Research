ClearAll[e,w,f];
data=Get[FileNameJoin[{DirectoryName[$InputFileName],"full_charts_results.wl"}]];
params={{1/16,1,1},{1/32,3/4,2}};
rows=Flatten[Table[With[{r=data["rows"][[i]],p=params[[j]]},<|"chart"->r["chart"],"term"->r["term"],"sample"->j,"value"->ToString[r["value"]/.{e->p[[1]],f->p[[2]],w->p[[3]]},InputForm]|>],{j,1,2},{i,1,Length[data["rows"]]}],1];
Print["Wolfram ",$Version," PID ",$ProcessID," samples ",Length[rows]];
Export[FileNameJoin[{DirectoryName[$InputFileName],"full_charts_samples.json"}],rows,"RawJSON"];
