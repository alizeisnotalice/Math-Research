Clear[u,H,a2,a3,z,w];
z=u^(-1/2)-H u^(1/2)+a2 u^(3/2)+a3 u^(5/2); w=z^2;
(* Tail Mills polynomial 1-1/w+3/w^2-15/w^3; H=log L+a. *)
res=Normal@Series[w+2 Log[u w]+2H-1/u-2 Log[1-1/w+3/w^2-15/w^3],{u,0,2}];
coef2=Solve[Coefficient[res,u,1]==0,a2][[1]];
coef3=Solve[(Coefficient[res,u,2]/.coef2)==0,a3][[1]];
mu=1+ProductLog[Exp[-1]]; sigma=Sqrt[mu]; aa=Log[sigma Sqrt[2Pi]];
rows=Table[Module[{roots,hh,approx},roots=Table[zz/.FindRoot[Log[Erfc[zz/Sqrt[2]]/2]-Log[sigma]-Log[zz]+ll/2==0,{zz,N[Sqrt[ll],wp]},WorkingPrecision->wp,AccuracyGoal->Floor[wp/2],PrecisionGoal->Floor[wp/2]],{wp,{40,80}}];hh=Log[ll]+aa;approx=Sqrt[ll]-hh/Sqrt[ll]+(2hh-1-hh^2/2)/ll^(3/2);
<|"log_n"->ll,"z_80"->ToString[N[roots[[2]],65]],"precision_doubling_difference"->ToString[N[Abs[roots[[1]]-roots[[2]]],10]],"next_term_z_error"->ToString[N[roots[[2]]-approx,20]]|>],{ll,{100,1000,10000}}];
result=<|"wolfram_version"->$Version,"a2"->ToString[a2/.coef2,InputForm],"a3"->ToString[a3/.coef3,InputForm],"rows"->rows|>;
Export[FileNameJoin[{DirectoryName[$InputFileName],"root_symbolic_checks.json"}],result,"RawJSON"];
Print[ExportString[result,"RawJSON"]];
