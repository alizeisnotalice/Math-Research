base=DirectoryName[$InputFileName];mu=1+ProductLog[Exp[-1]];
calc[n_]:=Module[{p=100,m,th,a,eps,sc,cdf,cdf0,med,med0,c0,cross,chi,rev},m=mu n-1;th=1+n^(-1/4);a=m/th;eps=Log[n]/Sqrt[n];sc=th {1,1-2/Sqrt[a],1+1/Sqrt[a]};
 cdf[x_?NumericQ]:=(1-eps) GammaRegularized[N[a,p],0,x/N[sc[[1]],p]]+eps(GammaRegularized[N[a,p],0,x/N[sc[[2]],p]]/3+2GammaRegularized[N[a,p],0,x/N[sc[[3]],p]]/3);
 cdf0[x_?NumericQ]:=GammaRegularized[N[a,p],0,x/N[th,p]];
 med=x/.FindRoot[cdf[x]==1/2,{x,N[m,p]},WorkingPrecision->p,AccuracyGoal->65,PrecisionGoal->65];med0=x/.FindRoot[cdf0[x]==1/2,{x,N[m,p]},WorkingPrecision->p,AccuracyGoal->65,PrecisionGoal->65];
 c0=Sqrt[2Pi mu](1/2-(CDF[NormalDistribution[],2]/3+2CDF[NormalDistribution[],-1]/3));
 cross=-(a-1)eps(Log[1-2/Sqrt[a]]/3+2Log[1+1/Sqrt[a]]/3);
 chi=((1-4/a)^(-a)-1)/3+2((1-1/a)^(-a)-1)/3;rev=eps^2 chi/(2(1-eps)^2);
 <|"n"->n,"precision"->p,"epsilon"->ToString[N[eps,60],InputForm],"all_scales_gt_1"->TrueQ[And@@Thread[sc>1]],"median_shift_over_log_n"->ToString[N[(med-med0)/Log[n],60],InputForm],"median_minus_mu_n"->ToString[N[med-mu n,60],InputForm],"predicted_limit"->ToString[N[c0,60],InputForm],"reverse_KL_upper"->ToString[N[rev,60],InputForm],"entropy_cross_term"->ToString[N[cross,60],InputForm],"cdf_median_residual"->ToString[N[Abs[cdf[med]-1/2],60],InputForm]|>
];
LaunchKernels[12];DistributeDefinitions[mu,calc];res=ParallelMap[calc,2^Range[12,34,2]];CloseKernels[];
Export[FileNameJoin[{base,"gamma_obstruction_wolfram.json"}],<|"status"->"complete","version"->$Version,"kernels"->12,"scope"->"Abstract legal rank densities only, not actual rectangular catalogues", "results"->res|>,"RawJSON"];Print[<|"cases"->Length[res],"all_scales_gt_1"->And@@Lookup[res,"all_scales_gt_1"]|>];Quit[];
