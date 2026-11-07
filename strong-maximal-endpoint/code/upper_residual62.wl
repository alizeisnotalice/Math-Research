ClearAll[residual62];
residual62[{nn_,off_}]:=Module[{wp=precision62,mu,ll,cc,pc,mc,vc,ks,ww,ts,rows,yy,pdf,den,e1,e2,ell,sc,sec,ms,start,lg,lgf},
start=AbsoluteTime[];mu=N[1+ProductLog[1/E],wp];ll=Log[N[nn,wp]];cc=mu+off Sqrt[ll/nn];pc=1/(1+Exp[cc]);mc=1+cc pc;vc=1+cc^2 pc(1-pc);
ks=Range[Max[0,Floor[nn pc-Sqrt[110 nn]]],Min[nn,Ceiling[nn pc+Sqrt[110 nn]]]];
lg=N[LogGamma[N[nn,wp]],wp];lgf=N[LogGamma[N[nn+1,wp]],wp];
ww=Exp[lgf-LogGamma[N[ks+1,wp]]-LogGamma[N[nn-ks+1,wp]]-cc ks-nn Log[1+Exp[-cc]]];
ms=-2pc(1-pc)+cc pc(1-pc)(1-2pc);
ts=nn mc+Subdivide[-2Sqrt[ll],2Sqrt[ll],256] Sqrt[nn vc];
rows=Table[yy=tt-cc ks;pdf=Map[If[#>0,Exp[(nn-1)Log[#]-#-lg],0]&,yy];den=ww.pdf;e1=(ww ks).(pdf/yy);e2=(ww ks^2).(pdf/yy^2);ell=(ww.(pdf((nn-1)/yy-1)))/den;sc=nn pc-(nn-1)e1/den;sec=nn pc((nn+1)pc-1)-2nn pc(nn-1)e1/den+(nn-1)(nn-2)e2/den;{Abs[sc]/ll,Abs[sec+nn ms ell]/ll^2},{tt,ts}];
<|"N"->nn,"offset"->off,"precision"->wp,"support_size"->Length[ks],"omitted_mass_Hoeffding_bound"->ToString[N[2Exp[-220],30],InputForm],"central_parameter_score_over_logN_max"->ToString[N[Max[rows[[All,1]]],35],InputForm],"central_residual_over_logN2_max"->ToString[N[Max[rows[[All,2]]],35],InputForm],"seconds"->N[AbsoluteTime[]-start]|>];
precision62=If[Length[$ScriptCommandLine]>=2,ToExpression[$ScriptCommandLine[[2]]],80];
outpath62=If[Length[$ScriptCommandLine]>=3,$ScriptCommandLine[[3]],FileNameJoin[{$TemporaryDirectory,"strong-maximal-residual62.json"}]];
LaunchKernels[4];ParallelEvaluate[$HistoryLength=0];DistributeDefinitions[residual62,precision62];
rows=ParallelMap[residual62,{{8192,-2},{8192,0},{8192,2},{32768,0}},Method->"FinestGrained"];
Export[outpath62,<|"version"->$Version,"status"->"complete","rows"->rows,"scope"->"80-digit numerical residual verification with Hoeffding-truncated binomial support; not interval certification"|>,"RawJSON"];CloseKernels[];Print["Residual checks completed"];
