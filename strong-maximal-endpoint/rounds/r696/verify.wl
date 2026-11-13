mu=1+ProductLog[1/E];
ratio=-(Log[mu])/2+x^2(1-1/mu)/2;
candidate=FullSimplify[ratio/.x^2->mu];
out=<|"gaussian_relative_entropy"->ToString[candidate,InputForm],
"value"->ToString[N[candidate,30],InputForm]|>;
Print[out];Export[FileNameJoin[{DirectoryName[$InputFileName],"wolfram.json"}],out,"RawJSON"];
