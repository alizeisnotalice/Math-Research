root=DirectoryName[DirectoryName[$InputFileName]];
Clear[s,t];
weights={{1,1/2,1/4},{1,1,1/2}};
tails={(8/3)Exp[-t]-2Exp[-2t]+(1/3)Exp[-4t],2t Exp[-t]+Exp[-2t]};
checks=MapThread[Function[{a,tail},Module[{densityLaplace},
 densityLaplace=LaplaceTransform[-D[tail,t],t,s];
 <|"weights"->ToString[a,InputForm],"tail"->ToString[tail,InputForm],
 "density_laplace_residual"->ToString[FullSimplify[densityLaplace-Times@@(1/(1+# s)&/@a),Assumptions->s>0],InputForm],
 "tail_at_zero"->ToString[tail/.t->0,InputForm]|>
]],{weights,tails}];
result=<|"wolfram_version"->$Version,"checks"->checks|>;
Export[FileNameJoin[{root,"rounds","r006","symbolic_checks.json"}],result,"RawJSON"];
Print[ExportString[result,"RawJSON"]];
