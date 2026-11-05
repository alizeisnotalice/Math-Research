Clear[n,c,q]; mu=1+ProductLog[Exp[-1]];sigma=Sqrt[mu];
logI[n_,wp_]:=N[Log[n]+n-1+LogGamma[n]+Log[GammaRegularized[n,n-1]]-n Log[n-1],wp];
(* Original weights: independent of the positive-tilted binomial implementation. *)
logavg[n_Integer,cc_?NumericQ,qq_?NumericQ,wp_Integer]:=Module[{pp=N[1/(1+Exp[cc]),wp],rate=(n-1)/n},
 -qq/n-n Log[rate]+Log[Sum[Binomial[n,k]pp^k(1-pp)^(n-k)Exp[cc k/n]GammaRegularized[n,rate Max[qq-cc k,0]],{k,0,n}]]];
root=DirectoryName[DirectoryName[$InputFileName]];
inputs=Import[FileNameJoin[{root,"rounds","r017","checks.json"}],"RawJSON"];
rows=Table[Module[{cc,qq,approx,residuals,roots,ii},
 cc=N[mu-sigma Sqrt[Log[nn]/nn]/2,110];
 (* Import is an Association, so access the observed boundary directly. *)
 qq=SelectFirst[inputs["rows"],#["n"]==nn&]["shifted_q"];
 residuals=Table[logavg[nn,N[cc,wp],SetPrecision[qq,wp],wp],{wp,{60,100}}];
 roots=Table[x/.FindRoot[logavg[nn,N[cc,wp],x,wp]==0,{x,SetPrecision[qq,wp]},WorkingPrecision->wp,AccuracyGoal->Floor[wp/2],PrecisionGoal->Floor[wp/2]],{wp,{60,100}}];
 ii=Table[logI[nn,wp],{wp,{60,100}}];
 approx=N[Log[2Pi nn]/2-1-nn Log[1-1/nn]+Log[GammaRegularized[nn,nn-1]]+1/(12nn)-1/(360nn^3),100];
 <|"n"->nn,"input_q_is_MATLAB_double"->True,"independent_root_80"->ToString[N[roots[[2]],80],InputForm],"root_precision_doubling_stable_25_digits"->TrueQ[Abs[roots[[1]]-roots[[2]]]<10^-25],"log_average_60"->ToString[N[residuals[[1]],30],InputForm],"log_average_100"->ToString[N[residuals[[2]],30],InputForm],"precision_doubling_stable_40_digits"->TrueQ[Abs[residuals[[1]]-residuals[[2]]]<10^-40],"exact_log_I_80"->ToString[N[ii[[2]],80],InputForm],"Stirling_log_I_difference"->ToString[N[ii[[2]]-approx,30],InputForm]|>],{nn,{256,1024}}];
result=<|"wolfram_version"->$Version,"working_precisions"->{60,100},"method"->"full original-binomial upper incomplete Gamma sum","rows"->rows|>;
Export[FileNameJoin[{root,"rounds","r017","symbolic_checks.json"}],result,"RawJSON"];
Print[ExportString[result,"RawJSON"]];
