folder=DirectoryName[$InputFileName];
moment[k_,prec_]:=Module[{mu,r,val},mu=N[1+ProductLog[Exp[-1]],prec];r=mu-1;
val=(NIntegrate[Log[(1-v)/v]^k,{v,0,1/2,1},WorkingPrecision->prec,AccuracyGoal->prec-15,PrecisionGoal->prec-15]+NIntegrate[Log[(1+w)/w]^k,{w,0,r},WorkingPrecision->prec,AccuracyGoal->prec-15,PrecisionGoal->prec-15])/mu;
<|"k"->k,"precision"->prec,"raw_moment"->ToString[val,InputForm]|>];
LaunchKernels[8];DistributeDefinitions[moment];cases=Flatten[Table[{k,prec},{k,1,4},{prec,{50,80}}],1];res=ParallelMap[Function[z,TimeConstrained[MemoryConstrained[moment@@z,268435456,$Failed],120,$Aborted]],cases];CloseKernels[];
mu=N[1+ProductLog[Exp[-1]],100];
checks=<|"positive_density_transform"->FullSimplify[Exp[x]/(1+Exp[x])^2,Assumptions->Element[x,Reals]],"negative_branch_peak"->FullSimplify[(Exp[c]/(rr(Exp[c]-1)^2))/.c->Log[(1+rr)/rr],Assumptions->rr>0],"mean_formula"->ToString[mu+Log[mu]-1,InputForm],"window_length_ratios"->Table[ToString[N[Exp[u]Log[(1+Exp[-u])/(1-Exp[-u])],80],InputForm],{u,{1,2,4}}]|>;
checks["positive_density_transform"]=ToString[checks["positive_density_transform"],InputForm];checks["negative_branch_peak"]=ToString[checks["negative_branch_peak"],InputForm];Export[FileNameJoin[{folder,"density_wolfram.json"}],<|"version"->$Version,"moments"->res,"checks"->checks|>,"RawJSON"];Print[<|"moment_cases"->Length[res],"successful"->Count[res,_Association]|>];Quit[];
