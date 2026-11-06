(* Symbolic identities plus precision and mesh refinement, no interval claim. *)
dirV46=DirectoryName[$InputFileName];
$quadrant46LibraryMode=True;Get[FileNameJoin[{dirV46,"quadrant46.wl"}]];$quadrant46LibraryMode=False;
$stripe46LibraryMode=True;Get[FileNameJoin[{dirV46,"stripe46.wl"}]];$stripe46LibraryMode=False;
$stripeBounds46LibraryMode=True;Get[FileNameJoin[{dirV46,"stripe_bounds46.wl"}]];$stripeBounds46LibraryMode=False;
idsV46=<|
 "Phi_derivative"->TrueQ[FullSimplify[D[h(1+Log[h])^(nn-1),h]==(1+Log[h])^(nn-2)(nn+Log[h]),Assumptions->{h>1,nn>=2}]],
 "stripe_Jacobian"->TrueQ[Abs[Det[{{1/2,1/2},{1/2,-1/2}}]]==1/2],
 "stripe_budget"->TrueQ[FullSimplify[(ww Sqrt[n]/2)Exp[-(tau+Log[ww Sqrt[n]/2])]==Exp[-tau],Assumptions->{ww>0,n>0,tau>0}]],
 "variance_bound"->TrueQ[16/9+(13/10)^2/4==7921/3600&&7921/3600<9/4],
 "MGF_mean"->TrueQ[FullSimplify[D[-Log[1-e]+Log[(1+Exp[(e-1)c])/(1+Exp[-c])],e]/.e->0]==1+c/(1+Exp[c])],
 "rate_ratio_derivative"->TrueQ[FullSimplify[D[1/d+1/(1+Exp[d]),d]<0,Assumptions->d>0]]|>;
aV46=case46[256,1,60];bV46=case46[256,1,100];sdV46=Abs[ToExpression[aV46["probability_internal"]]-ToExpression[bV46["probability_internal"]]];
refV46=Table[Module[{wp=120,mu,ll,z,t,r1,r2,q1,q2,dd},mu=N[1+ProductLog[1/E],320];ll=N[Log[n],260];z=zz/.FindRoot[Erfc[zz/Sqrt[2]]/2==Sqrt[mu/n]zz,{zz,Sqrt[ll]},WorkingPrecision->260,AccuracyGoal->220,PrecisionGoal->220];t=n mu-Sqrt[mu n]z;r1=n mu/2+Sqrt[(n/2)ll]/2;r2=t-r1;
 q1=logTailQ46[n/2,mu,r1,True,120,1/16,32];q2=logTailQ46[n/2,mu,r1,True,240,1/32,40];dd=Abs[q1["log_probability"]-q2["log_probability"]];
 <|"N_exact"->ToString[n,InputForm],"B"->1/2,"log_upper_tail_difference"->ToString[dd,InputForm],"less_than_1e-30"->TrueQ[dd<10^-30],"both_roots_valid"->TrueQ[q1["root_valid"]&&q2["root_valid"]]|>],{n,{256,2^64,2^128}}];
bd1V46=runB46[2^4096,1];bd2V46=runB46[2^4096,1,2];
resV46=<|"wolfram"->$Version,"identities"->idsV46,"stripe_probability_WP60_to_100_difference"->ToString[sdV46,InputForm],"stripe_difference_less_than_1e-24"->TrueQ[sdV46<10^-24],"tail_WP_mesh_cutoff_refinements"->refV46,"stripe_bound_WP2548_to_5096_stored_23digit_rates_agree"->TrueQ[bd1V46["negative_log_probability_upper_bound_rate"]==bd2V46["negative_log_probability_upper_bound_rate"]],"status"->"completed: exact identities and floating refinements; no total interval certificate"|>;
outV46=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/verify46.json"];Export[outV46,resV46,"RawJSON"];Print[resV46];
