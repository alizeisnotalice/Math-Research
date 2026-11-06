Clear[c,b,x];
smallPref=(1+Exp[-c])^-b;tailPref=(2/(1+Exp[-c]))^b;
resSmall=FullSimplify[D[smallPref,c]-b/(1+Exp[c])smallPref,Assumptions->{b>=1,c>0}];
resTail=FullSimplify[D[tailPref,c]-b/(1+Exp[c])tailPref,Assumptions->{b>=1,c>0}];
surv[cc_,xx_]:=If[xx<cc,(Exp[-cc]+Exp[-xx])/(1+Exp[-cc]),2Exp[-xx]/(1+Exp[-cc])];
rows=Table[<|"z"->zz,"survival_c1"->N[surv[1,zz],30],"survival_c2"->N[surv[6/5,zz],30],"difference_c2_minus_c1"->N[surv[6/5,zz]-surv[1,zz],30]|>,{zz,{1/2,2}}];
result=<|"wolfram"->$Version,"small_CDF_prefactor_derivative_residual"->resSmall,"large_tail_prefactor_derivative_residual"->resTail,"single_step_survival_crossing"->rows|>;
Export["/Users/zhengzhihao/Desktop/T/work/r032_order_checks.json",result,"RawJSON"];Print[result];
