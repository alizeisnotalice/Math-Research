(* Actual Orlicz-budget indicator stripes: continuous finite-mixture probabilities. *)
mu46=N[1+ProductLog[1/E],360];wp46=60;
mg46[e_]:=-Log[1-e]+Log[(1+Exp[(e-1)mu46])/(1+Exp[-mu46])];
z46[n_,prec_]:=z46[n,prec]=z/.FindRoot[Erfc[z/Sqrt[2]]/2==Sqrt[mu46/n]z,{z,Sqrt[Log[n]]},WorkingPrecision->prec,AccuracyGoal->prec-20,PrecisionGoal->prec-20];
case46[n_,kap_,prec_:60,width_:4]:=Module[{b=n/2,ll=Log[n],tau,ss,ld,ud,ks,bi,ww,den,pdf,interval,x0,x1,hi,cent,sd,breaks,terms,pr,start=AbsoluteTime[]},
 tau=n mu46-Sqrt[mu46 n]z46[n,prec+40];ss=tau+Log[width Sqrt[n]/2];ld=N[kap Sqrt[n ll],prec+20];ud=ld+width Sqrt[n];
 If[!TrueQ[ss>ud],Return[<|"N"->n,"kappa"->N[kap],"status"->"radial threshold below stripe; simple exact budget formula not used"|>]];
 ks=Range[0,b];bi=N[Binomial[b,#]&/@ks,prec+20];ww=bi Exp[-mu46 ks]/(1+Exp[-mu46])^b;den=Gamma[b](1+Exp[-mu46])^b;
 pdf[t_?NumericQ]:=If[t<=0,0,Exp[-t]Dot[bi,Clip[t-mu46 ks,{0,Infinity}]^(b-1)]/den];
 interval[lo_?NumericQ,up_?NumericQ]:=If[up<=lo,0,Dot[ww,GammaRegularized[b,Clip[lo-mu46 ks,{0,Infinity}],Clip[up-mu46 ks,{0,Infinity}]]]];
 x0=(ss+ld)/2;x1=(ss+ud)/2;hi=(b mg46[1/4]+50Log[10])/(1/4);cent=b mu46;sd=Sqrt[b mu46];breaks=Sort[DeleteDuplicates@Join[{x1,hi},Select[cent+sd Range[-8,12,2],x1<#<hi&]]];
 pr=NIntegrate[pdf[x]interval[ss-x,x-ld],{x,x0,x1},WorkingPrecision->prec,AccuracyGoal->(Floor[prec/2]-5),PrecisionGoal->(Floor[prec/2]-8),MaxRecursion->16,Method->{"GlobalAdaptive","SymbolicProcessing"->0}];
 terms=Table[NIntegrate[pdf[x]interval[x-ud,x-ld],{x,breaks[[j]],breaks[[j+1]]},WorkingPrecision->prec,AccuracyGoal->(Floor[prec/2]-5),PrecisionGoal->(Floor[prec/2]-8),MaxRecursion->16,Method->{"GlobalAdaptive","SymbolicProcessing"->0}],{j,Length[breaks]-1}];pr+=Total[terms];
 <|"N"->n,"kappa"->N[kap],"radial_threshold"->N[ss,22],"difference_lower"->N[ld,22],"difference_upper"->N[ud,22],"width"->width,"normalized_actual_Orlicz_budget"->N[Exp[tau]width Sqrt[n]/2 Exp[-ss],22],"actual_probability"->N[pr,22],"probability_internal"->ToString[pr,InputForm],"observed_negative_log_probability_over_logN"->N[-Log[pr]/ll,22],"asymptotic_exponent"->N[kap^2/(2mu46),22],"analytic_omitted_probability_bound"->N[Exp[b mg46[1/4]-hi/4],22],"working_precision"->prec,"seconds"->N[AbsoluteTime[]-start]|>];
If[!TrueQ[$stripe46LibraryMode],args46=$ScriptCommandLine;out46=If[Length[args46]>=2,args46[[2]],"/tmp/stripe46.json"];rows46={};Do[row=case46[n,kap];AppendTo[rows46,row];Print[KeyTake[row,{"N","kappa","actual_probability","observed_negative_log_probability_over_logN","seconds"}]];Export[out46,<|"wolfram"->$Version,"status"->"actual indicator-stripe budget and finite-mixture probability; floating quadrature not interval certificate","rows"->rows46|>,"RawJSON"],{n,{64,256}},{kap,{1/2,1,2}}];];
