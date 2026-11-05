(* Independent high precision and exact checks; no numerical certificate claim. *)
base=DirectoryName[DirectoryName[$InputFileName]];
s=N[BesselI[0,2 Sqrt[2]],100];
s50=N[BesselI[0,2 Sqrt[2]],50];
(* Coefficients of S'(a)-S(a), S(a)=Sum (2a)^k/(k!)^2. *)
coefficientCheck=And@@Table[Simplify[2^(k+1)/((k+1)! k!)-2^k/(k!)^2+(k-1)/(k+1)*2^k/(k!)^2]==0,{k,1,20}];
tailBound=1/3+1/9+25/828;
(* H_k integral series and independent 80-digit numerical integral. *)
residuals=Table[With[{ell=Log[1+1/n]},Abs[NIntegrate[Exp[x] x^(k-1)/(k-1)!,{x,0,ell},WorkingPrecision->80,AccuracyGoal->65,PrecisionGoal->65]-N[ell^k/k! Sum[ell^j/j!*k/(k+j),{j,0,80}],80]]],{n,{2,4,16,256}},{k,1,8}];
result=<|"wolfram"->$Version,"precision"->100,"S"->ToString[s,InputForm],"S_over_e"->ToString[N[s/E,90],InputForm],"precision_doubling_difference"->ToString[Abs[s-s50],InputForm],"coefficient_check"->coefficientCheck,"exact_derivative_loss_bound"->ToString[tailBound,InputForm],"derivative_positive_bound"->TrueQ[tailBound<1],"max_integral_residual"->ToString[Max[residuals],InputForm]|>;
Export[FileNameJoin[{base,"rounds","r024","symbolic_checks.json"}],result,"RawJSON"];
Print[result];
