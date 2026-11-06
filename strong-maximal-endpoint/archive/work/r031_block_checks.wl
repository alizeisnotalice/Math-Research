Clear[t,c];mu=1+ProductLog[Exp[-1]];
i2[prec_]:=Module[{cc=N[mu,prec+30],pp,rr,B,pieces},rr=Exp[-cc];pp=rr/(1+rr);B=1+Exp[cc];pieces={
 NIntegrate[Exp[-t]t(2pp)^2/(1+rr)^2,{t,0,cc},WorkingPrecision->prec,AccuracyGoal->prec-25,PrecisionGoal->prec-25],
 NIntegrate[Exp[-t](3t-2cc)(2pp(1-B/(3t-2cc)))^2/(1+rr)^2,{t,cc,2cc},WorkingPrecision->prec,AccuracyGoal->prec-25,PrecisionGoal->prec-25],
 NIntegrate[Exp[-t](4t-4cc)(2pp(1-2B/(4t-4cc)))^2/(1+rr)^2,{t,2cc,Infinity},WorkingPrecision->prec,AccuracyGoal->prec-25,PrecisionGoal->prec-25]};Total[pieces]];
p[c_]:=1/(1+Exp[c]);m[c_]:=1+c p[c];v[c_]:=1+c^2 p[c](1-p[c]);
mp=Together[D[m[c],c]];vp=Together[D[v[c],c]];
mpResidual=FullSimplify[(mp/.c->x)/.Exp[x]->1/(x-1)];vpResidual=FullSimplify[((vp/.c->x)/.Exp[x]->1/(x-1))-(x-1)];
mppResidual=FullSimplify[((Together[D[m[c],{c,2}]]/.c->x)/.Exp[x]->1/(x-1))+(x-1)/x];
normalMoment=Integrate[(z^2-1)^2 Exp[-z^2/2]/Sqrt[2Pi],{z,-Infinity,Infinity}];
k3=2+c^3 r(1-r)(1-2r);k4=6+c^4 r(1-r)(1-6r(1-r));
k3Residual=FullSimplify[(k3/.r->(c-1)/c)-c(3-c)];
k4Residual=FullSimplify[(k4/.r->(c-1)/c)-c(3-c)(4-c)];
gramResidual=FullSimplify[(c(3-c)(4-c)-(c(3-c))^2/c)-c(3-c)];
result=<|"wolfram"->$Version,"I2_60digits"->ToString[N[i2[60],35],InputForm],"I2_90digits"->ToString[N[i2[90],65],InputForm],"mean_derivative_critical_residual"->ToString[mpResidual,InputForm],"second_mean_derivative_critical_residual"->ToString[mppResidual,InputForm],"variance_derivative_critical_residual"->ToString[vpResidual,InputForm],"normal_centered_square_moment"->normalMoment,"critical_cumulant_3_residual"->k3Residual,"critical_cumulant_4_residual"->k4Residual,"critical_Gram_residual"->gramResidual,"stationary_fisher_limit_70digits"->ToString[N[p[mu]^2/2,70],InputForm]|>;
Export["/Users/zhengzhihao/Desktop/T/work/r031_block_symbolic.json",result,"RawJSON"];Print[result];
