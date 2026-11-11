a90=Rest[$ScriptCommandLine];out90=a90[[1]];workers90=If[Length[a90]>=2,ToExpression[a90[[2]]],12];mem90=If[Length[a90]>=3,ToExpression[a90[[3]]],256];If[!DirectoryQ[out90],CreateDirectory[out90,CreateIntermediateDirectories->True]];
ClearAll[price90,H90,q90,str90,run90];price90[x_,nn_]:=Log[1+x]+(nn-1)Log[1+Log[1+x]];H90[u_]:=(2u^2+7u+7)/(1+u)^3;q90[x_]:=Erfc[x/Sqrt[2]]/2;str90[x_]:=ToString[N[x,65],InputForm];
symbolic90=<|
"price_first_derivative"->TrueQ[FullSimplify[D[price90[x,nn],x]/.x->0]==nn],
"price_second_derivative"->TrueQ[FullSimplify[D[price90[x,nn],{x,2}]/.x->0]==1-2nn],
"price_third_derivative"->TrueQ[FullSimplify[D[price90[x,nn],{x,3}]/.x->0]==7nn-5],
"third_derivative_formula"->TrueQ[FullSimplify[D[price90[x,nn],{x,3}]==(2+(nn-1)H90[Log[1+x]])/(1+x)^3,Assumptions->x>=0&&nn>=1]],
"decreasing_H_identity"->TrueQ[FullSimplify[D[H90[u],u]==-(2u^2+10u+14)/(1+u)^4,Assumptions->u>=0]],
"positive_H"->TrueQ[FullSimplify[H90[u]>0,Assumptions->u>=0]],
"H_bound_polynomial"->TrueQ[FullSimplify[7(1+u)^3-(2u^2+7u+7)==7u^3+19u^2+14u]],
"H_le_seven"->TrueQ[FullSimplify[H90[u]<=7,Assumptions->u>=0]],
"critical_cubic_scale"->TrueQ[FullSimplify[nn (2ss Sqrt[ll/nn])^3==8ss^3 ll^(3/2)/Sqrt[nn],Assumptions->nn>0&&ll>0&&ss>0]],
"quadratic_marked_price"->TrueQ[FullSimplify[nn (Sqrt[ll/nn]u)-(nn-1/2)(Sqrt[ll/nn]u)^2==Sqrt[nn ll]u-(1-1/(2nn))ll u^2,Assumptions->nn>0&&ll>0]],
"CDF_grid_balance"->TrueQ[FullSimplify[(ll^2/(hh Sqrt[nn]))/hh/.hh->ll nn^(-1/4),Assumptions->nn>0&&ll>0]==1],
"ramp_perturbation_slope"->TrueQ[FullSimplify[(hh/(2Sqrt[nn]))Sqrt[nn]/hh,Assumptions->nn>0&&hh>0]==1/2]
|>;Print["symbolic ",symbolic90];
run90[k_Integer]:=Module[{rows={},wp,prec,nn=2^k,ll,mu,ss,x,exactprice,rem,cubic,zz,hh,de,rat,res},Do[
prec=wp+Ceiling[k Log[2]/Log[10]]+80;ll=N[k Log[2],prec];mu=u/.FindRoot[u-1-Exp[-u]==0,{u,N[128/100,prec]},WorkingPrecision->prec,AccuracyGoal->prec-20,PrecisionGoal->prec-20];ss=Sqrt[mu];x=2Sqrt[mu ll/N[nn,prec]];
(* Independent direct high-guard subtraction, no integral remainder or Python interval path. *)
exactprice=price90[x,nn];rem=exactprice-(nn x-(nn-1/2)x^2);cubic=(7nn-5)x^3/6;
zz=u/.FindRoot[Log[q90[u]]-Log[ss u]+ll/2==0,{u,N[4Sqrt[ll]/5,prec],N[11Sqrt[ll]/10,prec]},WorkingPrecision->prec,AccuracyGoal->prec-20,PrecisionGoal->prec-20,MaxIterations->100];res=Abs[Log[q90[zz]]-Log[ss zz]+ll/2];hh=ll/N[nn,prec]^(1/4);de=ss zz/Sqrt[N[nn,prec]];rat=(ll^2/(hh Sqrt[N[nn,prec]]))/hh;
AppendTo[rows,<|"working_precision"->wp,"guarded_working_precision"->prec,"mu"->str90[mu],"sigma"->str90[ss],"L"->str90[ll],"x_max"->str90[x],"log_remainder"->str90[Log[rem]],"log_cubic_bound"->str90[Log[cubic]],"remainder_to_cubic_bound"->str90[rem/cubic],"z"->str90[zz],"h"->str90[hh],"delta0"->str90[de],"delta0_over_h"->str90[de/hh],"normalized_identity_value"->str90[rat],"log_tail_residual"->str90[res],"pass"->TrueQ[0<x<1/4&&rem>0&&rem<=cubic&&Abs[rat-1]<10^-wp&&res<10^-wp]|>],{wp,{80,160}}];<|"k"->k,"kernel_id"->$KernelID,"precision_rows"->rows|>];
ks90={16,24,32,48,64,96,128,192,256,384,512,768,1024};LaunchKernels[workers90];ParallelEvaluate[$HistoryLength=0];DistributeDefinitions[price90,H90,q90,str90,run90,mem90];tm90=AbsoluteTime[];rows90=ParallelMap[MemoryConstrained[TimeConstrained[run90[#],120,"TIMEOUT"],mem90 2^20,"MEMORY_LIMIT"]&,ks90,Method->"FinestGrained"];CloseKernels[];pass90=And@@Values[symbolic90]&&AllTrue[rows90,AssociationQ[#]&&AllTrue[#["precision_rows"],TrueQ[#["pass"]]&]&];Export[FileNameJoin[{out90,"verification90.json"}],<|"version"->$Version,"pass"->TrueQ[pass90],"configured_workers"->workers90,"per_case_memory_limit_MiB"->mem90,"kernels_closed"->True,"seconds"->N[AbsoluteTime[]-tm90],"symbolic"->symbolic90,"rows"->rows90,"scope"->"Scalar Taylor-price and normal-grid scales only; no numerical proof of the physical theorem and no finite-N physical onset certificate."|>,"RawJSON"];Print[If[TrueQ[pass90],"PASS ","FAIL "],Length[rows90]," price/scale cases"];If[!TrueQ[pass90],Exit[3]];
