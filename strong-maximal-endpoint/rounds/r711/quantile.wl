folder=DirectoryName[$InputFileName];
quant[n_,prec_]:=Module[{mu,p,w,cdf,lo,hi,mid,target},
 mu=N[1+ProductLog[Exp[-1]],prec];p=(mu-1)/mu;
 w=Table[Binomial[n,j]p^j(1-p)^(n-j),{j,0,n}];
 cdf[x_?NumericQ]:=Sum[If[x>mu j,w[[j+1]]GammaRegularized[n+1,0,x-mu j],0],{j,0,n}];
 target=1/2+1/n;lo=mu n-10;hi=mu n+20;
 If[!(cdf[lo]<target<cdf[hi]),Return[<|"n"->n,"precision"->prec,"status"->"bracket_failed"|>]];
 Do[mid=(lo+hi)/2;If[cdf[mid]<target,lo=mid,hi=mid],{Ceiling[prec Log[2,10]]-12}];mid=(lo+hi)/2;
 <|"n"->n,"precision"->prec,"status"->"converged","q"->ToString[mid,InputForm],"offset"->ToString[mid-mu n,InputForm],"cdf_residual"->ToString[cdf[mid]-target,InputForm],"width"->ToString[hi-lo,InputForm]|>];
LaunchKernels[8];DistributeDefinitions[quant];
cases=Flatten[Table[{n,prec},{n,{4,8,16,32,64,128,256,512}},{prec,{50,80}}],1];
r=ParallelMap[Function[z,TimeConstrained[MemoryConstrained[quant@@z,536870912,$Failed],180,$Aborted]],cases,Method->"CoarsestGrained"];CloseKernels[];
out=<|"version"->$Version,"results"->r,"median_offset_limit"->ToString[N[(3+(1+ProductLog[Exp[-1]]))/6,100],InputForm],"integer_star_efficiency"->ToString[N[Log[3]/3,100],InputForm],"star_penalty"->ToString[N[1+ProductLog[Exp[-1]]-Log[3]/3,100],InputForm],"exact_skewness_identity"->FullSimplify[2+(m-1)(2-m)-m(3-m)]|>;
Export[FileNameJoin[{folder,"quantile_wolfram.json"}],out,"RawJSON"];Print[<|"cases"->Length[r],"converged"->Count[r,a_Association/;a["status"]=="converged"]|>];Quit[];
