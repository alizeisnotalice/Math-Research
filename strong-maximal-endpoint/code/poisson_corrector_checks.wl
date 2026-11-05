Clear[c,n,z,x];p[c_]:=1/(1+Exp[c]);m[c_]:=1+c*p[c];
drift=FullSimplify[D[m[c],c]+p[c]*(c-m[c])];
smooth[u_]:=35u^4-84u^5+70u^6-20u^7;
cc=1/4;big=6;theta[z_]:=Piecewise[{{smooth[(z-cc)/cc],cc<z<2cc},{1,2cc<=z<=big/2},{1-smooth[(z-big/2)/(big/2)],big/2<z<big}},0];
mu=1+ProductLog[Exp[-1]];gauss[z_]:=Exp[-z^2/(2mu)]/Sqrt[2Pi mu];
(* Explicit pieces avoid evaluation ambiguity in symbolic iterators. *)
cm[prec_]:=2 Total[Table[NIntegrate[z^2 theta[z](2-theta[z])gauss[z],{z,pair[[1]],pair[[2]]},WorkingPrecision->prec,AccuracyGoal->prec-25,PrecisionGoal->prec-25],{pair,{{cc,2cc},{2cc,big/2},{big/2,big}}}]];
harmonicLimit=Limit[(Log[n/Log[n]^2]-Log[Log[n]^4])/Log[n],n->Infinity];
ramp=Exp[x/n]-1; rampDerivative=D[ramp,x];
initialDerivative=rampDerivative/.x->0; finalDerivative=rampDerivative/.x->n Log[2];
softSup=FullSimplify[n finalDerivative,Assumptions->n>0];
softVariation=FullSimplify[n(initialDerivative+(finalDerivative-initialDerivative)+finalDerivative),Assumptions->n>0];
g=(mu-1)/8;coeff=FullSimplify[(9-mu)/24+g-(3+mu)/12];
result=<|"wolfram"->$Version,"drift_cancellation_residual"->ToString[drift,InputForm],"harmonic_window_limit"->harmonicLimit,"soft_derivative_sup_times_N"->softSup,"soft_derivative_total_variation_times_N"->softVariation,"R_log_coefficient_identity_residual"->ToString[coeff,InputForm],"cutoff_Gaussian_second_moment_60digits"->ToString[N[cm[60],40],InputForm],"cutoff_Gaussian_second_moment_90digits"->ToString[N[cm[90],70],InputForm],"mu_80digits"->ToString[N[mu,80],InputForm],"g_80digits"->ToString[N[g,80],InputForm],"R_lower_coefficient_80digits"->ToString[N[(3+mu)/12,80],InputForm]|>;
base=DirectoryName[DirectoryName[$InputFileName]];Export[FileNameJoin[{base,"rounds","r029","symbolic_checks.json"}],result,"RawJSON"];Print[result];
