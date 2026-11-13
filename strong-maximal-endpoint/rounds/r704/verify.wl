pc=(u-1)/u;
m2=FullSimplify[-2 pc(1-pc)+u pc(1-pc)(1-2pc)];
v1=FullSimplify[2u pc(1-pc)-u^2 pc(1-pc)(1-2pc)];
F=m2 t^2/2-v1 t/(2 Sqrt[u]);
peak=t/.First[Solve[D[F,t]==0,t]];
mu=1+ProductLog[1/E];
out=<|"mean_second_derivative"->ToString[m2,InputForm],
"variance_first_derivative"->ToString[v1,InputForm],
"peak_t"->ToString[FullSimplify[peak,Assumptions->u>1],InputForm],
"peak_gain"->ToString[FullSimplify[F/.t->peak,Assumptions->u>1],InputForm],
"log_coefficient"->ToString[FullSimplify[(3-u)/6+(u-1)/8],InputForm],
"values"->ToString[N[{peak,(u-1)/8,(9-u)/24,(u-1-Log[u])/2}/.u->mu,25],InputForm]|>;
Print[out];Export[FileNameJoin[{DirectoryName[$InputFileName],"wolfram.json"}],out,"RawJSON"];
