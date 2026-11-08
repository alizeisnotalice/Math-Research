(* Explicit inputs.json, output directory, optional worker count. *)
args=Rest[$ScriptCommandLine];If[Length[args]<2,Print["Need inputs.json and output directory"];Exit[2]];
in66=Import[args[[1]],"RawJSON"];out66=args[[2]];
If[!DirectoryQ[out66],CreateDirectory[out66,CreateIntermediateDirectories->True]];
workers66=If[Length[args]>=3,ToExpression[args[[3]]],12];
ClearAll[second66];
second66[spec_]:=Module[{M=spec[[1]],factor=spec[[2]],wp,mu,c,p,r,psi,aa,bb,small,P,A,B,u,target,start},
start=AbsoluteTime[];wp=factor(2Ceiling[Log[10,M^2+1]]+80);mu=c0/.FindRoot[c0==1+Exp[-c0],{c0,N[5/4,wp]},WorkingPrecision->wp];c=mu;p=1/(1+Exp[c]);
psi[t_]:=(1+Exp[-c+I c t])/((1+Exp[-c])(1-I t));aa[t_]:=p+Exp[-c+I c t](-1+I t)/(1+Exp[-c+I c t]);bb[t_]:=-p(1-p)+Exp[-c+I c t](-1+I t)^2/(1+Exp[-c+I c t])^2;
u=-N[1/M,wp];small=Exp[M^2(Log[psi[u]]+I mu/M)];P=small(Exp[-I mu M]+Exp[-mu])/((1+Exp[-mu])(1-I M));A=aa[M]+M^2 aa[u];B=bb[M]+M^2 bb[u];target=p Exp[-mu/2];
<|"M"->ToString[M],"b"->ToString[M^2+1],"precision"->wp,"factor"->factor,"first_abs"->ToString[N[Abs[P A],35],InputForm],"second_abs_over_M"->ToString[N[Abs[P(A^2+B)]/M,35],InputForm],"limit"->ToString[N[target,35],InputForm],"first_limit_complex_residual"->ToString[N[Abs[P A+target],25],InputForm],"second_limit_complex_residual"->ToString[N[Abs[P(A^2+B)/M+I target],25],InputForm],"seconds"->N[AbsoluteTime[]-start]|>];
specs=({ToExpression[#["M"]],#["factor"]}&/@in66["second_specs"]);
LaunchKernels[workers66];DistributeDefinitions[second66];start=AbsoluteTime[];rows=ParallelMap[second66,specs,Method->"FinestGrained"];
Export[FileNameJoin[{out66,"second66.json"}],<|"version"->$Version,"workers"->$KernelCount,"seconds"->N[AbsoluteTime[]-start],"rows"->rows,"scope"->"Numerical corroboration of analytic limits; no directed rounding"|>,"RawJSON"];CloseKernels[];Print["Completed second-response checks"];
