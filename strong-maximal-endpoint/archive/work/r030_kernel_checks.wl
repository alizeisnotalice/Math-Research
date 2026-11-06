Clear[c,t,r,b,k];
(* On each knot interval kc<t<(k+1)c; full tail uses k=b. *)
res=Flatten[Table[With[{poly=Sum[Binomial[bb,j](t-j c)^(bb-1),{j,0,kk}],shift=Sum[Binomial[bb-1,j](t-(j+1)c)^(bb-2),{j,0,kk-1}]},FullSimplify[D[poly,c]+bb(bb-1)shift]],{bb,2,8},{kk,0,bb}]];
(* Scalar Orlicz derivative for h>1, using y=Log[h]. *)
phiDerivative=FullSimplify[D[h(1+Log[h])^(n-1),h]/.h->Exp[y],Assumptions->{y>0,n>=2}];
profileResidual=FullSimplify[phiDerivative-(1+y)^(n-2)(n+y),Assumptions->{y>0,n>=2}];
result=<|"wolfram"->$Version,"knot_polynomial_cases"->Length[res],"all_gradient_residuals_zero"->And@@(#===0&/@res),"soft_profile_derivative_residual"->ToString[profileResidual,InputForm]|>;
Export["/Users/zhengzhihao/Desktop/T/work/r030_kernel_symbolic.json",result,"RawJSON"];Print[result];
