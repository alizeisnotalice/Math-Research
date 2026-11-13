mu=1+ProductLog[1/E];
out=<|"cross_entropy_series"->ToString[Series[z/e-(e^(-2)-1)Log[1+e z],{e,0,2}],InputForm],
"source_entropy_limit"->ToString[N[Log[2 Pi E mu]/2,30],InputForm],
"cross_entropy_limit"->ToString[N[Log[2 Pi]/2+mu/2,30],InputForm],
"relative_entropy_limit"->ToString[N[(mu-1-Log[mu])/2,30],InputForm]|>;
Print[out];Export[FileNameJoin[{DirectoryName[$InputFileName],"wolfram.json"}],out,"RawJSON"];
