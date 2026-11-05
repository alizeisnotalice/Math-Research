Clear[h,nn,c,t,z,gg];
phi=h*(1+Log[h])^(nn-1);
phiDerivative=FullSimplify[D[phi,h]-(1+Log[h])^(nn-2)*(nn+Log[h]),Assumptions->{h>1,nn>=2}];
f2[t_]:=Exp[-t]*(t+2 Piecewise[{{t-c,t>c}},0]+Piecewise[{{t-2c,t>2c}},0])/(1+Exp[-c])^2;
transform=FullSimplify[Integrate[Exp[-z*t]*f2[t],{t,0,Infinity},Assumptions->{c>0,z>0}]-(1+Exp[-(1+z)*c])^2/((1+z)^2*(1+Exp[-c])^2),Assumptions->{c>0,z>0}];
mu=1+ProductLog[Exp[-1]];
(* For r=0, b+b Gamma blocks, conditional law is s*Beta(b,b). *)
conditionalRows=Table[{2b,ToString[N[Sqrt[2b]*Exp[LogGamma[2b]-2LogGamma[b]-(2b-2)Log[2]-Log[2b]],80],InputForm]},{b,{1,2,4,8,16,32,64,128}}];
scaleIdentity=FullSimplify[(1+(cc-1)aa)/(aa)-(1/aa+cc-1)];
result=<|"wolfram"->$Version,"precision"->80,"Phi_derivative_residual"->ToString[phiDerivative,InputForm],"two_increment_density_Laplace_residual"->ToString[transform,InputForm],"Gamma_conditional_columns"->{"n","sqrt_n_times_sup_conditional_density_at_s_n"},"Gamma_conditional_rows"->conditionalRows,"Gamma_conditional_limit"->ToString[N[Sqrt[2/Pi],80],InputForm],"tilted_mean_rate_range"->{"1/3","2 mu"},"scale_algebra_residual"->ToString[scaleIdentity,InputForm]|>;
base=DirectoryName[DirectoryName[$InputFileName]];Export[FileNameJoin[{base,"rounds","r026","symbolic_checks.json"}],result,"RawJSON"];Print[result];
