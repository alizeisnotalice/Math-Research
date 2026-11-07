(* Independent 60-digit counterparts of the MATLAB kernel grid. *)
ClearAll[case62];
case62[{nn_Integer,off_Integer}]:=Module[{wp=60,mu,ll,rho,cc,zz,ks,ww,vals,x,yg,den,der,cdf},
 mu=N[1+ProductLog[1/E],wp];ll=Log[N[nn,wp]];rho=Sqrt[ll/nn];cc=mu+(off/4)rho;
 zz=z/.FindRoot[Erfc[z/Sqrt[2]]/2==Sqrt[mu/nn]z,{z,Sqrt[ll]},WorkingPrecision->wp,AccuracyGoal->40,PrecisionGoal->40];
 ks=Range[0,nn];ww=N[Table[Binomial[nn,k]Exp[-cc k]/(1+Exp[-cc])^nn,{k,ks}],wp];
 vals=Table[x=nn mu-Sqrt[mu nn]zz+sh ll;yg=x-cc ks;
 den=Total[MapThread[If[#2>0,#1 Exp[(nn-1)Log[#2]-#2-LogGamma[nn]],0]&,{ww,yg}]];
 der=Total[MapThread[If[#2>0,#1 Exp[(nn-1)Log[#2]-#2-LogGamma[nn]]((nn-1)/#2-1),0]&,{ww,yg}]];
 cdf=Total[MapThread[If[#2>0,#1 GammaRegularized[nn,0,#2],0]&,{ww,yg}]];
 <|"shift"->sh,"density"->ToString[N[den,35],InputForm],"N_density_over_logN"->ToString[N[nn den/ll,35],InputForm],"CDF"->ToString[N[cdf,35],InputForm],"score"->ToString[N[der/den,35],InputForm]|>,{sh,-4,4}];
 <|"N"->nn,"offset_quarters"->off,"precision"->wp,"weight_mass_residual"->ToString[N[Total[ww]-1,30],InputForm],"critical"->vals|>];
outdir=If[Length[$ScriptCommandLine]>=2,Last[$ScriptCommandLine],FileNameJoin[{$TemporaryDirectory,"strong-maximal-r062-wolfram"}]];
If[!DirectoryQ[outdir],CreateDirectory[outdir,CreateIntermediateDirectories->True]];
LaunchKernels[4];ParallelEvaluate[$HistoryLength=0];DistributeDefinitions[case62,outdir];
tasks=Flatten[Table[{nn,off},{off,-8,8},{nn,{4096,2048,1024,512,256,128,64,32}}],1];
start=AbsoluteTime[];rows=ParallelMap[Function[t,Module[{r=case62[t]},Export[FileNameJoin[{outdir,"wl_"<>ToString[t[[1]]]<>"_"<>ToString[t[[2]]]<>".json"}],r,"RawJSON"];r]],tasks,Method->"FinestGrained"];
result=<|"status"->"complete","wolfram"->$Version,"subkernels"->$KernelCount,"seconds"->N[AbsoluteTime[]-start],"rows"->rows,"scope"->"High-precision Gamma-Binomial kernel grid, independent of MATLAB; not a proof of an endpoint upper bound"|>;
Export[FileNameJoin[{outdir,"wolfram_upper62.json"}],result,"RawJSON"];CloseKernels[];Print[KeyTake[result,{"status","subkernels","seconds"}]];
