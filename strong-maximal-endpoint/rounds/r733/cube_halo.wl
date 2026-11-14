base=DirectoryName[$InputFileName];
run[task_]:=Module[{nn=task[[1]],wp=task[[2]],fun,val},
 fun[s_?NumericQ]:=Exp[s]LaguerreL[nn-1,1,-2s];
 val=1+2NIntegrate[fun[s],{s,0,N[Log[2],wp+20]},WorkingPrecision->wp+20,AccuracyGoal->wp,PrecisionGoal->wp,Method->{"GlobalAdaptive","SymbolicProcessing"->0},MaxRecursion->20];
 <|"n"->nn,"target_digits"->wp,"guard_digits"->20,"log_halo"->ToString[N[Log[val],wp],InputForm]|>
];
LaunchKernels[12];DistributeDefinitions[run];res=ParallelMap[run,Flatten[Table[{nn,wp},{nn,{1,2,4,8,16,32,64,128}},{wp,{40,70}}],1]];CloseKernels[];
Export[FileNameJoin[{base,"cube_halo_wolfram.json"}],<|"status"->"complete","version"->$Version,"kernels"->12,"method"->"Independent positive Laguerre integral, not finite hypergeometric sum","results"->res|>,"RawJSON"];
Print[<|"status"->"complete","cases"->Length[res]|>];Quit[];
