(* Analytic Chebyshev/tilting bounds for the genuine budget stripe w=4.
   Formula bounds are rigorous; decimal evaluations here are not intervals. *)
mgB46[c_,e_]:=-Log[1-e]+Log[(1+Exp[(e-1)c])/(1+Exp[-c])];
yB46[c_,e_]:=1/(1-e)+c/(1+Exp[(1-e)c]);
runB46[n_,kap_,factor_:1]:=Module[{wp=factor Max[100,2IntegerLength[n,10]+80],mu,ll,zz,dd,tau,rad,ld,ud,th,thU,psi,psiU,q,loglo,logup,start=AbsoluteTime[]},
 mu=N[1+ProductLog[1/E],wp+80];ll=N[Log[n],wp+30];zz=z/.FindRoot[Erfc[z/Sqrt[2]]/2==Sqrt[mu/n]z,{z,Sqrt[ll]},WorkingPrecision->wp+20,AccuracyGoal->wp-20,PrecisionGoal->wp-20];dd=Sqrt[mu n]zz;tau=n mu-dd;rad=tau+Log[2Sqrt[n]];ld=kap Sqrt[n ll];ud=ld+4Sqrt[n];
 th=t/.FindRoot[(n/2)(yB46[mu,t]-yB46[mu,-t])==ld+2Sqrt[n],{t,(ld+2Sqrt[n])/(n mu)},WorkingPrecision->wp+20,AccuracyGoal->wp-20,PrecisionGoal->wp-20];
 thU=kap/mu Sqrt[ll/n];psi=(n/2)(mgB46[mu,th]+mgB46[mu,-th]);psiU=(n/2)(mgB46[mu,thU]+mgB46[mu,-thU]);q=1-9/16-(9/4)n/(dd-Log[2Sqrt[n]])^2;
 If[!(TrueQ[q>0]&&TrueQ[th<=1/4]&&TrueQ[thU<1]&&TrueQ[dd>Log[2Sqrt[n]]]&&TrueQ[rad>ud]),Return[<|"N_exact"->ToString[n,InputForm],"kappa"->N[kap],"status"->"finite sufficient conditions fail"|>]];
 loglo=psi-th ud+Log[q];logup=Min[0,psiU-thU ld];
 <|"N_exact"->ToString[n,InputForm],"kappa"->N[kap],"width"->4,"working_precision"->wp,"theta_mid"->ToString[N[th,23],InputForm],"root_relative_residual"->ToString[((n/2)(yB46[mu,th]-yB46[mu,-th])-ld-2Sqrt[n])/(ld+2Sqrt[n]),InputForm],"Q_event_lower_Chebyshev"->N[q,23],"negative_log_probability_lower_bound_rate"->N[-logup/ll,23],"negative_log_probability_upper_bound_rate"->N[-loglo/ll,23],"limiting_rate"->N[kap^2/(2mu),23],"log_probability_lower_formula"->N[loglo,23],"log_probability_upper_formula"->N[logup,23],"status"->"completed: analytic formulas evaluated numerically","seconds"->N[AbsoluteTime[]-start]|>];
If[!TrueQ[$stripeBounds46LibraryMode],argsB46=$ScriptCommandLine;outB46=If[Length[argsB46]>=2,argsB46[[2]],"/tmp/stripe_bounds46.json"];rowsB46={};Do[row=runB46[n,kap];AppendTo[rowsB46,row];Print[KeyTake[row,{"N_exact","kappa","negative_log_probability_lower_bound_rate","negative_log_probability_upper_bound_rate","seconds"}]];Export[outB46,<|"wolfram"->$Version,"status"->"analytic stripe probability bounds; high precision decimal evaluation, no directed interval claim","rows"->rowsB46|>,"RawJSON"],{n,2^{16,64,256,1024,4096}},{kap,{1/2,1,2}}];
];
