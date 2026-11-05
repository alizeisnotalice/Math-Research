root=DirectoryName[DirectoryName[$InputFileName]];
Clear[c,lam];
p=1/(1+Exp[c(1+lam)]);
m=1/(1+lam)+c p;
v=1/(1+lam)^2+c^2 p(1-p);
logMGF=-Log[1+lam]+Log[1+Exp[-c(1+lam)]]-Log[1+Exp[-c]];
entropy=-lam m-logMGF;
checks=<|"mean_derivative_plus_variance"->FullSimplify[D[m,lam]+v,Assumptions->{c>=0,lam>=0}],
 "entropy_derivative_minus_lambda_variance"->FullSimplify[D[entropy,lam]-lam v,Assumptions->{c>=0,lam>=0}],
 "variance_derivative_identity"->FullSimplify[D[v,lam]+2/(1+lam)^3+c^3 p(1-p)(1-2p),Assumptions->{c>=0,lam>=0}]|>;
result=<|"wolfram_version"->$Version,"identity_residuals"->Map[ToString[#,InputForm]&,checks],
 "variance_upper_constant"->ToString[1+4/E^2,InputForm],
 "variance_derivative_upper_constant"->ToString[2+27/E^3,InputForm]|>;
Export[FileNameJoin[{root,"rounds","r012","symbolic_checks.json"}],result,"RawJSON"];
Print[ExportString[result,"RawJSON"]];
