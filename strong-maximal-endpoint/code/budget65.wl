(* Independent actual capped-ramp modular, with analytic concave-log tail bound. *)
args=Rest[$ScriptCommandLine];If[Length[args]<2,Print["Need inputs.json and output.json"];Exit[2]];
in=Import[args[[1]],"RawJSON"];outPath=args[[2]];
ClearAll[one65];
one65[row_]:=Module[{calc,n=Round[row["N"]],al,vals},
calc[wp_]:=Module[{a=ToExpression[row["slope_decimal"]<>"`"<>ToString[wp]],u0,em,der,int,tail,cap,v},
u0=Sqrt[N[n,wp]] (20+4Sqrt[Log[N[n,wp]]]);
em[u_]:=Log[1+a u/n]+(n-1)Log[1+Log[1+a u/n]]-u;
der=a/(n+a u0)(1+(n-1)/(1+Log[1+a u0/n]))-1;
int=NIntegrate[Exp[em[v]],{v,0,u0},WorkingPrecision->wp,AccuracyGoal->Floor[wp/2],PrecisionGoal->Floor[wp/2],MaxRecursion->25];
tail=Exp[em[u0]]/(-der);cap=Exp[Log[2]+(n-1)Log[1+Log[2]]-n/a];
<|"precision"->wp,"log_budget_factor"->ToString[N[Log[int+cap],35],InputForm],"soft_tail_relative_bound"->ToString[N[tail/(int+cap),20],InputForm],"tail_derivative"->ToString[N[der,20],InputForm]|>];
vals=calc/@{80,160};<|"N"->n,"profile"->row["profile"],"checks"->vals|>];
LaunchKernels[8];DistributeDefinitions[one65];t=AbsoluteTime[];rows=ParallelMap[one65,in,Method->"FinestGrained"];Export[outPath,<|"version"->$Version,"workers"->$KernelCount,"seconds"->N[AbsoluteTime[]-t],"rows"->rows,"scope"->"High precision quadrature and analytic tail formula evaluated without directed rounding"|>,"RawJSON"];CloseKernels[];Print["Completed modular checks"];
