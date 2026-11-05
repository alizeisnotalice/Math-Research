Clear[c,t,r,b,k];
(* On each knot interval kc<t<(k+1)c; full tail uses k=b. *)
res=Flatten[Table[With[{poly=Sum[Binomial[bb,j](t-j c)^(bb-1),{j,0,kk}],shift=Sum[Binomial[bb-1,j](t-(j+1)c)^(bb-2),{j,0,kk-1}]},FullSimplify[D[poly,c]+bb(bb-1)shift]],{bb,2,8},{kk,0,bb}]];
(* Scalar Orlicz derivative for h>1, using y=Log[h]. *)
phiDerivative=FullSimplify[D[h(1+Log[h])^(n-1),h]/.h->Exp[y],Assumptions->{y>0,n>=2}];
profileResidual=FullSimplify[phiDerivative-(1+y)^(n-2)(n+y),Assumptions->{y>0,n>=2}];
result=<|"wolfram"->$Version,"knot_polynomial_cases"->Length[res],"all_gradient_residuals_zero"->And@@(#===0&/@res),"soft_profile_derivative_residual"->ToString[profileResidual,InputForm]|>;
kernelResult=result;
Clear[a,h];rr=Exp[-a];ss=Exp[-a-h];
dPlus=Log[(1+rr)/(1+ss)]-Log[2](rr-ss)/(1+ss);
sm=Exp[-a+h];dMinus=Log[(1+rr)/(1+sm)]+2Log[2](sm-rr)/(1+sm);
p=rr/(1+rr);
slopes=FullSimplify[{Limit[dPlus/h,h->0],Limit[dMinus/h,h->0]}/p,Assumptions->a>0];
tvPlus=(rr-ss)(2/(1+rr)-1/(1+ss));tvSlope=FullSimplify[Limit[tvPlus/h,h->0]/p,Assumptions->a>0];
mu=1+ProductLog[Exp[-1]];
rows=Table[<|"delta"->ToString[hh,InputForm],"D_plus_over_pdelta"->ToString[N[(dPlus/(p h))/.{a->mu,h->hh},50],InputForm],"D_minus_over_pdelta"->ToString[N[(dMinus/(p h))/.{a->mu,h->hh},50],InputForm]|>,{hh,{1/1000,1/10^6,1/10^9,1/10^12}}];
result=<|"wolfram"->$Version,"precision"->"exact symbolic inputs; N[...,50]","KL_slopes_divided_by_p"->ToString[slopes,InputForm],"TV_slope_divided_by_p"->tvSlope,"rows"->rows|>;

entropyResult=result; result=<|"kernel"->kernelResult,"entropy"->entropyResult|>;
base=DirectoryName[DirectoryName[$InputFileName]]; Export[FileNameJoin[{base,"rounds","r030","symbolic_checks.json"}],result,"RawJSON"]; Print[result];
