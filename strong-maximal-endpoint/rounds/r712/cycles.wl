folder=DirectoryName[$InputFileName];
cycle[t_,prec_]:=Module[{mu,r,fun,pts,out},
 mu=N[1+ProductLog[Exp[-1]],prec];r=mu-1;
 fun[v_?NumericQ]:=Module[{L},If[v==0,Return[1-Log[Abs[t](1+r)]]];L=(1+r)(1+t(1/v-1));-Log[Abs[v]]+Which[L==0,0,L==1,1,True,1+(1-L)Log[Abs[1-L]]/L]];
 pts=Sort[DeleteDuplicates[Join[{-r,0,1},Select[If[t==1,{}, {-t/(1-t)}],-r<#<1&],Select[{(1+r)t/((1+r)t-r)},-r<#<1&]]]];
 out=NIntegrate[fun[v],Evaluate[Prepend[pts,v]],WorkingPrecision->prec,AccuracyGoal->prec-15,PrecisionGoal->prec-15,MaxRecursion->30]/(1+r);
 <|"theta"->ToString[t,InputForm],"precision"->prec,"mean"->ToString[out,InputForm],"deficit"->ToString[2mu-out,InputForm]|>];
LaunchKernels[12];DistributeDefinitions[cycle];cases=Flatten[Table[{sg 10^j,prec},{sg,{-1,1}},{j,{-6,-3,-1,0,1,3,6}},{prec,{50,80}}],2];
res=ParallelMap[Function[z,TimeConstrained[MemoryConstrained[cycle@@z,268435456,$Failed],180,$Aborted]],cases,Method->"CoarsestGrained"];CloseKernels[];
mat=DiagonalMatrix[{v1,v2}]+(IdentityMatrix[2]-DiagonalMatrix[{v1,v2}]).{{0,a},{b,0}};
out=<|"version"->$Version,"determinant"->ToString[Expand[Det[mat]],InputForm],"determinant_identity"->FullSimplify[Det[mat]-(v1 v2-a b(1-v1)(1-v2))],"mu"->ToString[N[1+ProductLog[Exp[-1]],100],InputForm],"results"->res|>;
Export[FileNameJoin[{folder,"cycles_wolfram.json"}],out,"RawJSON"];Print[<|"cases"->Length[res],"success"->Count[res,_Association]|>];Quit[];
