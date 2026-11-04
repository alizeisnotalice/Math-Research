Print["Wolfram version: ",$Version];
LaunchKernels[15];
Clear[calc];calc[n_,lam_,kap_]:=Module[{delta,q,q0,eq,cc,ee,dd,series,slackcoef,fullcoef},
delta=u/(1+2u)-cc u^3+ee u^4;eq=(2delta)^n+kap lam u^(n+2)-(2u(1-2delta))^n;
series=Normal[Series[eq,{u,0,n+3}]];cc=FullSimplify[cc/.First[Solve[Coefficient[series,u,n+2]==0,cc]]];
ee=FullSimplify[ee/.First[Solve[(Coefficient[series,u,n+3]/.cc->cc)==0,ee]]];
dd=u/(1+2u)-cc u^3+ee u^4;q=(2dd)^n+lam u^(n+2);q0=(2u/(1+2u))^n;
fullcoef=FullSimplify[Coefficient[Normal[Series[q/(1-q)-q0/(1-q0),{u,0,n+2}]],u,n+2]];
{n,lam,kap,"c"->cc,"next"->ee,"hypotheticalFullGainCoefficient"->fullcoef,"expected"->(1-kap)lam,"check"->FullSimplify[fullcoef==(1-kap)lam],"scope"->"full capture requires separate geometric proof; kappa=2 is candidate guarantee under explicit gap assumptions"}];
specs=Join[Flatten[Table[{2,lam,kap},{lam,{1/16,1/4,1}},{kap,{1,3/2,2}}],1],Flatten[Table[{3,lam,kap},{lam,{1/16,1/4,1}},{kap,{3/2,2}}],1]];

DistributeDefinitions[calc];
results=ParallelMap[Function[sp,calc@@sp],specs,Method->"FinestGrained"];
Print[InputForm[results]];
ok=And@@(TrueQ["check"/.Cases[#,_Rule]]& /@ results);
Print["case_count=",Length[results],"; all_checks=",ok];
CloseKernels[];
Exit[If[Length[results]===15 && TrueQ[ok],0,1]];
