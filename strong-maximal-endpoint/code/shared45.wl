(* Genuine nonnegative 2D input, restricted common horizontal side vs
   a within-vertical-average relaxation. No full-maximal-operator claim. *)
Clear[r,s,eps,yy];r0=1/2;beta=1/4;f0[t_]:=1-beta(t-r0)^2/2;gg[t_]:=t-r0;
fplus[ee_,y_]:=f0[0]+ee y gg[0];
fminus[t_,ee_,y_]:=1-beta(t-r0)^2/2-beta(1+t)(t-r0)+ee y(1+2t-r0);
avgExact=FullSimplify[(fplus[eps,yy]+Integrate[fminus[s,eps,yy],{s,0,r}])/(1+r)-(f0[r]+eps yy gg[r])];
minCheck=Resolve[ForAll[{s,eps},0<=s<=1&&0<=eps<=1/8,And[fminus[s,eps,1]>=13/32,fminus[s,eps,-1]>=13/32]],Reals];
plusCheck=Resolve[ForAll[eps,0<=eps<=1/8,And[fplus[eps,1]>=29/32,fplus[eps,-1]>=29/32]],Reals];
maxShared=FullSimplify[f0[r0]];rOpt[y_]:=r0+eps y/beta;
maxAdaptive=FullSimplify[(f0[rOpt[1]]+eps gg[rOpt[1]]+f0[rOpt[-1]]-eps gg[rOpt[-1]])/2];
phi2[x_?NumericQ]:=x(1+Max[0,Log[x]]);
cost45[e_,prec_:60]:=Module[{extra=0,roots,ed,poly,lo,hi,mid,z,y},
 Do[poly=fminus[s,e,y];roots=N[s/.Solve[poly==1,s,Reals],prec+20];ed=Join[{0},Sort[Select[roots,0<#<1&]],{1}];
 Do[{lo,hi}=ed[[j;;j+1]];mid=(lo+hi)/2;If[(poly/.s->mid)>1,extra+=NIntegrate[Evaluate[poly Log[poly]],{s,lo,hi},WorkingPrecision->prec,AccuracyGoal->Floor[3prec/5],PrecisionGoal->Floor[prec/2],MaxRecursion->12]/2],{j,Length[ed]-1}];
 z=fplus[e,y];If[z>1,extra+=z Log[z]/2],{y,{-1,1}}];N[31/16+extra,prec]];
rows45=Table[<|"epsilon"->ToString[e,InputForm],"shared_value"->1,"within_average_relaxation_value"->N[1+2e^2,25],"gap_exact"->ToString[2e^2,InputForm],"actual_Phi2_denominator"->N[cost45[e],25],"input_mass"->31/16|>,{e,{1/128,1/64,1/32,1/16,1/8}}];
precCheck45=<|"epsilon"->"1/8","precision60"->N[cost45[1/8,60],25],"precision120"->N[cost45[1/8,120],25],"absolute_difference"->ToString[cost45[1/8,60]-cost45[1/8,120],InputForm]|>;
checks45=<|"fiber_average_identity_residual"->avgExact,"negative_fiber_nonnegative_certificate"->minCheck,"positive_fiber_nonnegative_certificate"->plusCheck,"shared_maximum"->maxShared,"adaptive_maximum"->ToString[maxAdaptive,InputForm],"nonproduct_determinant"->ToString[FullSimplify[fplus[eps,1]fminus[1/2,eps,-1]-fplus[eps,-1]fminus[1/2,eps,1]],InputForm],"gap_coefficient_residual"->FullSimplify[maxAdaptive-maxShared-eps^2/(2beta)]|>;
out45=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/shared45.json"];Export[out45,<|"wolfram"->$Version,"scope"->"restricted rectangles [-r,1] x [0,1] vs horizontal r selected inside y-average; not full M_S","checks"->checks45,"working_precision"->60,"method"->"Phi2 split at exact f=1 roots; only smooth positive log pieces integrated","precision_check"->precCheck45,"rows"->rows45|>,"RawJSON"];Print[checks45];
