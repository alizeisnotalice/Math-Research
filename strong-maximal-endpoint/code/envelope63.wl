ClearAll[case63];
case63[{nn_Integer,off_Integer,nodes_Integer}]:=Module[{wp=precision63,mu,ll,cc,pc,mc,vc,z,dd,ss,ks,ww,lg,lgf,ts,den,der,yy,pdf,delta,v,U,dt,m1,m2,score,start,trap},
start=AbsoluteTime[];mu=N[1+ProductLog[1/E],wp];ll=Log[N[nn,wp]];cc=mu+off Sqrt[ll/nn];pc=1/(1+Exp[cc]);mc=1+cc pc;vc=1+cc^2 pc(1-pc);
z=zz/.FindRoot[Erfc[zz/Sqrt[2]]/2==Sqrt[mu/nn]zz,{zz,Sqrt[ll]},WorkingPrecision->wp,AccuracyGoal->60,PrecisionGoal->60];dd=Sqrt[mu nn]z;ss=nn mu-dd;
ks=Range[Max[0,Floor[nn pc-Sqrt[110 nn]]],Min[nn,Ceiling[nn pc+Sqrt[110 nn]]]];lg=LogGamma[N[nn,wp]];lgf=LogGamma[N[nn+1,wp]];ww=Exp[lgf-LogGamma[N[ks+1,wp]]-LogGamma[N[nn-ks+1,wp]]-cc ks-nn Log[1+Exp[-cc]]];
ts=nn mc+Subdivide[-3Sqrt[ll],3Sqrt[ll],nodes-1] Sqrt[nn vc];dt=ts[[2]]-ts[[1]];
den={};der={};Do[yy=tt-cc ks;pdf=Map[If[#>0,Exp[(nn-1)Log[#]-#-lg],0]&,yy];AppendTo[den,ww.pdf];AppendTo[der,ww.(pdf((nn-1)/yy-1))],{tt,ts}];
delta=Map[Max[#,0]&,ts-ss+20ll];v=delta/nn;Do[v=v-(v+(nn-1)Log[1+v]-delta)/(1+(nn-1)/(1+v)),{12}];U=Exp[v]-1;
trap[x_]:=dt(Total[x]-(First[x]+Last[x])/2);m1=trap[den U];m2=trap[den U^2];score=trap[der U];
<|"N"->nn,"offset"->off,"nodes"->nodes,"working_precision"->wp,"EU_over_dN"->ToString[N[m1/(dd/nn),35],InputForm],"EU2_over_dN2"->ToString[N[m2/(dd/nn)^2,35],InputForm],"E_U_location_score"->ToString[N[score,35],InputForm],"cap_inactive"->TrueQ[Max[U]<1],"seconds"->N[AbsoluteTime[]-start]|>];
outdir=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],FileNameJoin[{$TemporaryDirectory,"strong-maximal-r063"}]];
precision63=If[Length[$ScriptCommandLine]>=3,ToExpression[$ScriptCommandLine[[3]]],80];
tier63=If[Length[$ScriptCommandLine]>=4,$ScriptCommandLine[[4]],"small"];
If[!DirectoryQ[outdir],CreateDirectory[outdir,CreateIntermediateDirectories->True]];
ns63=If[tier63=="large",{131072,524288},{8192,32768}];LaunchKernels[12];ParallelEvaluate[$HistoryLength=0];DistributeDefinitions[case63,outdir,precision63];
tasks=Flatten[Table[{nn,off,nodes},{nn,ns63},{off,{-1,0,1}},{nodes,{1025,2049}}],2];
start=AbsoluteTime[];rows=ParallelMap[Function[t,Module[{r=case63[t]},Export[FileNameJoin[{outdir,"wl_"<>ToString[precision63]<>"_"<>StringRiffle[ToString/@t,"_"]<>".json"}],r,"RawJSON"];r]],tasks,Method->"FinestGrained"];
Export[FileNameJoin[{outdir,"envelope63_"<>tier63<>"_"<>ToString[precision63]<>".json"}],<|"status"->"complete","version"->$Version,"workers"->$KernelCount,"seconds"->N[AbsoluteTime[]-start],"rows"->rows,"scope"->"High precision numerical trapezoidal envelope moments, not an interval certificate"|>,"RawJSON"];CloseKernels[];Print["Completed envelope63 Wolfram checks"];
