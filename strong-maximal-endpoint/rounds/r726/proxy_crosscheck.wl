base=DirectoryName[$InputFileName];
main=Import[FileNameJoin[{base,"nonlinear_search.json"}],"RawJSON"];
neg=Import[FileNameJoin[{base,"negative_search.json"}],"RawJSON"];
mu=N[1+ProductLog[Exp[-1]],70];
num[x_String]:=Module[{s=StringSplit[ToLowerCase[x],"e"]},ToExpression[s[[1]]<>"`60"<>If[Length[s]==2,"*^"<>s[[2]],""]]];
num[x_]:=SetPrecision[x,60];
mean[c_]:=1+c/(1+Exp[c]);
phi[c_,t_]:=(1+Exp[-c+I t c])/((1+Exp[-c])(1-I t));
check[item_]:=Module[{r=item[[2]],kind=item[[1]],abz,t,eps,cp,cm,saved,vals,delta,cent,mx,score,tol},
 abz=r["counts_plus_minus_zero"];t=num[r["frequency_t"]];eps=num[r["epsilon"]];
 cp=num/@r[If[kind=="main","best_proxy_cplus","extended_best_cplus"]];
 cm=num/@r[If[kind=="main","best_proxy_cminus","extended_best_cminus"]];
 saved=num/@r[If[kind=="main","best_proxy_values","extended_proxy_values"]];
 vals=Table[abz[[1]] mean[cp[[j]]]+abz[[2]] mean[cm[[j]]]+abz[[3]] mu-eps Re[Exp[I (j-1)Pi/32] phi[cp[[j]],t]^abz[[1]] phi[cm[[j]],-t]^abz[[2]]],{j,64}];
 delta=Max[Abs[vals-saved]];tol=Max[10^-7,(51/100) Max[(10^(Floor[Log[10,Abs[#]]]-9)&/@saved)]];cent=vals-mu r["n"];mx=Max[cent];score=mx+Log[Mean[Exp[cent-mx]]]-Log[BesselI[0,eps]];
 <|"kind"->kind,"n"->r["n"],"pattern"->r["pattern"],"epsilon"->r["epsilon"],"max_absolute_proxy_difference"->ToString[delta,InputForm],"geometric_logscore_minus_mu_n_64nodes"->ToString[score,InputForm],"saved_decimal_comparison_tolerance"->ToString[tol,InputForm],"difference_within_saved_decimal_tolerance"->TrueQ[delta<tol]|>
];
LaunchKernels[12];DistributeDefinitions[mu,num,mean,phi,check];
res=ParallelMap[check,Join[({"main",#}&/@main["results"]),({"negative",#}&/@neg["negative_side_results"])] ];
Export[FileNameJoin[{base,"proxy_crosscheck.json"}],<|"status"->"complete","precision_digits"->60,"input_precision_note"->"Saved controls, frequencies, and scores have10 significant decimal digits;60-digit arithmetic recomputes those decimals only. Tolerance includes half a saved score rounding unit plus angle/control rounding allowance.","kernels"->Length[Kernels[]],"cases"->Length[res],"scope"->"Recompute geometric proxy at saved decimal controls and 64 exact angle nodes; not a median or interval certificate.","all_saved_decimal_tolerances_pass"->And@@Lookup[res,"difference_within_saved_decimal_tolerance"],"results"->res|>,"RawJSON"];
Print[<|"cases"->Length[res],"all_tolerances_pass"->And@@Lookup[res,"difference_within_saved_decimal_tolerance"]|>];CloseKernels[];Quit[];
