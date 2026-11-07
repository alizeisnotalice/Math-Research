ClearAll[rampCase63];
rampCase63[r_]:=Module[{wp=precision63,n,off,L,mu,c,pc,q,ks,w,lg,raw,g,g2,Q,Q2,pdf,pdf2,dp,pos,pos2,h,aa,ac,acc,res,start},
start=AbsoluteTime[];n=Round[r["N"]];off=Round[r["offset"]];mu=N[1+ProductLog[1/E],wp];L=Log[N[n,wp]];c=mu+off Sqrt[L/n];pc=1/(1+Exp[c]);q=ToExpression[r["q_decimal"]<>"`"<>ToString[wp]];
ks=Range[Max[0,Floor[n pc-Sqrt[110n]]],Min[n,Ceiling[n pc+Sqrt[110n]]]];lg=LogGamma[N[n,wp]];
w=Exp[LogGamma[N[n+1,wp]]-LogGamma[N[ks+1,wp]]-LogGamma[N[n-ks+1,wp]]-c ks-n Log[1+Exp[-c]]];raw=q-c ks;g=raw;g2=raw+n;
Q=GammaRegularized[n,#]&/@g;Q2=GammaRegularized[n,#]&/@g2;pdf=Exp[(n-1)Log[g]-g-lg];pdf2=Exp[(n-1)Log[g2]-g2-lg];dp=(n-1)/g-1;
pos=n (GammaRegularized[n+1,#]&/@g)-raw Q;pos2=n (GammaRegularized[n+1,#]&/@g2)-(raw+n)Q2;
h=Q+(pos-pos2)/n;aa=n pc-ks;ac=ks pdf+ks/n (Q-Q2);acc=-ks^2 pdf dp+ks^2/n (pdf-pdf2);res=w.((aa^2-n pc(1-pc))h+2aa ac+acc);
<|"N"->n,"offset"->off,"q_decimal"->r["q_decimal"],"precision"->wp,"curvature"->ToString[N[res,35],InputForm],"mean_residual"->ToString[N[w.h-1,35],InputForm],"seconds"->N[AbsoluteTime[]-start]|>];
inputPath63=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"rounds/r063/validation63.json"];
outputPath63=If[Length[$ScriptCommandLine]>=3,$ScriptCommandLine[[3]],FileNameJoin[{$TemporaryDirectory,"ramp63-wolfram.json"}]];
precision63=If[Length[$ScriptCommandLine]>=4,ToExpression[$ScriptCommandLine[[4]]],80];
data63=Import[inputPath63,"RawJSON"];
inputs=If[ListQ[data63],data63,data63["high_precision_ramp"]["rows"]];
LaunchKernels[Min[12,Length[inputs]]];ParallelEvaluate[$HistoryLength=0];DistributeDefinitions[rampCase63,precision63];start=AbsoluteTime[];
rows=ParallelMap[rampCase63,inputs,Method->"FinestGrained"];
Export[outputPath63,<|"status"->"complete","version"->$Version,"workers"->$KernelCount,"seconds"->N[AbsoluteTime[]-start],"rows"->rows,"scope"->"80-digit curvature at exact listed decimal q, not a certified exact mean-one optimizer"|>,"RawJSON"];CloseKernels[];Print["Completed ramp curvature verification"];
