(* Constants for the general C_n second scale, derived from the existing sandwich. *)
mu=1+ProductLog[1/E];
values=Table[<|"precision"->prec,"mu"->ToString[N[mu,prec],InputForm],"sigma"->ToString[N[Sqrt[mu],prec],InputForm],
 "next_scale_constant"->ToString[N[Sqrt[mu]Log[2 Pi mu]/2,prec],InputForm],
 "defining_identity_residual"->ToString[With[{u=N[mu,prec+40]},N[u-1-Exp[-u],prec]],InputForm]|>,{prec,{80,160}}];
out=<|"version"->$Version,"values"->values,
 "theorem_scope"->"(mu*n-log C_n)/sqrt(n log n) tends to sqrt(mu); existing general sandwich r021/r044/r688. Finer log n coefficient remains undetermined."|>;
Export[FileNameJoin[{DirectoryName[$InputFileName],"second_constant.json"}],out,"RawJSON"];Print[values[[1]]];Quit[];
