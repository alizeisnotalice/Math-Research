root=DirectoryName[DirectoryName[$InputFileName]];
Clear[mu,lam,v,a,n,k];
K=lam mu-Log[1+lam]+Log[1+Exp[-mu(1+lam)]]-Log[1+Exp[-mu]];
p=1/(1+Exp[mu(1+lam)]);
third=FullSimplify[D[K,{lam,3}]+2/(1+lam)^3+mu^3 p(1-p)(1-2p),
 Assumptions->{mu>0,lam>=0}];
b=Exp[1+ProductLog[1/E]]; mm=Log[b];
derivatives=FullSimplify[{mu-1-mu rr/(1+rr),1+mu^2 rr/(1+rr)^2-mu}/.rr->mu-1];
integralChecks=Table[Module[{lhs,rhs},
 lhs=NIntegrate[Exp[-(1-aa/nn)vv](1+aa vv/nn)^(nn-1),{vv,0,Infinity},WorkingPrecision->50];
 rhs=Sum[Factorial[nn-1]/Factorial[nn-1-kk](aa/nn)^kk/(1-aa/nn)^(kk+1),{kk,0,nn-1}];
 <|"n"->nn,"a"->ToString[aa,InputForm],"integral50"->ToString[N[lhs,25],InputForm],
 "finite_sum_difference"->ToString[N[Abs[lhs-rhs],15],InputForm]|>],
 {nn,{4,16,64}},{aa,{1/2,1}}];
result=<|"wolfram_version"->$Version,"third_derivative_identity_residual"->ToString[third,InputForm],
 "Kprime0_and_Ksecond0_minus_mu_residuals"->(ToString[#,InputForm]&/@derivatives),
 "mu40"->ToString[N[mm,40],InputForm],
 "integral_checks"->Flatten[integralChecks]|>;
Export[FileNameJoin[{root,"rounds","r003","symbolic_checks.json"}],result,"RawJSON"];
Print[ExportString[result,"RawJSON"]];
