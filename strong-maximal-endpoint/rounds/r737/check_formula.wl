base=DirectoryName[$InputFileName];
qbound={FullSimplify[(5/6)^2*7/44>1/16],FullSimplify[2*(22/7)*(11/4)<25]};
LaunchKernels[12];
res=ParallelTable[Module[{n=h,V,cap,terms,peak,prec=80},terms=Table[LogGamma[n+1]-LogGamma[j+1]-LogGamma[n-j+1]+j Log[3/4]+(n-j) Log[1/4],{j,0,n/2}];terms=N[terms,prec];terms[[-1]]+=Log[1/2];peak=Max[terms];V=-(peak+Log[Total[Exp[terms-peak]]]);cap=N[n Log[2/Sqrt[3]],prec];<|"h"->h,"precision"->prec,"V"->ToString[V,InputForm],"fixed_capacity"->ToString[cap,InputForm],"centered_gap"->ToString[V-cap-Log[n]/2,InputForm]|>],{h,{100,200,400,800,1000,2000,4000,8000,10000,16000,24000,32000}}];
CloseKernels[];
Export[FileNameJoin[{base,"check_formula_wolfram.json"}],<|"version"->$Version,"kernels"->12,"quantile_lipschitz_rational_checks"->qbound,"spectral_diagnostics"->res|>,"RawJSON"];
Print[<|"quantile_checks"->qbound,"spectral_cases"->Length[res]|>];Quit[];
