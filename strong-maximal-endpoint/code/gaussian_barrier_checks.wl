Clear[a,b,u,cm,r,tau,xx,v];
normal[x_]:=Erfc[-x/Sqrt[2]]/2;
conv=normal[a/b]-Exp[-a+b^2/2]normal[a/b-b];
ode=FullSimplify[D[conv,a]-(normal[a/b]-conv),Assumptions->{Element[{a,b},Reals],b>0}];
coercivity=FullSimplify[u^2/2-cm r u/2+cm^2 r^2/8-(u-cm r/2)^2/2];
ratio=FullSimplify[Log[Sqrt[tau/(tau+1)]]+xx^2/(2v tau)-(xx-u)^2/(2v(tau+1))-
 (-Log[1+1/tau]/2+xx^2/(2v tau(tau+1))+(2xx u-u^2)/(2v(tau+1))),Assumptions->{tau>0,v>0}];
rows=Table[Module[{bb=20,aa=20zz,quad,formula},quad=NIntegrate[Exp[-y]normal[(aa-y)/bb],{y,0,Infinity},WorkingPrecision->80,AccuracyGoal->55,PrecisionGoal->55];formula=N[conv/.{a->aa,b->bb},80];
<|"z"->zz,"quadrature_formula_difference"->ToString[N[Abs[quad-formula],20],InputForm],"agrees_45_digits"->TrueQ[Abs[quad-formula]<10^-45]|>],{zz,{-3,0,3}}];
result=<|"wolfram_version"->$Version,"Exp_convolution_ODE_residual"->ToString[ode,InputForm],"mean_variance_square_residual"->ToString[coercivity,InputForm],"Gaussian_density_ratio_residual"->ToString[ratio,InputForm],"rows"->rows|>;
root=DirectoryName[DirectoryName[$InputFileName]];
Export[FileNameJoin[{root,"rounds","r020","symbolic_checks.json"}],result,"RawJSON"];
Print[ExportString[result,"RawJSON"]];
