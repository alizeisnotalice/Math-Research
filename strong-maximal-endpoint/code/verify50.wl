$hessian50LibraryMode=True;Get[FileNameJoin[{DirectoryName[$InputFileName],"hessian50.wl"}]];$hessian50LibraryMode=False;
logPbeta50=-Log[cf[bb]]/bb;dd150=D[logPbeta50,bb];dd250=D[logPbeta50,{bb,2}];
rep50={cf[bb]->cv,Derivative[1][cf][bb]->-ii1,Derivative[2][cf][bb]->ii2};
F50=FullSimplify[((dd150^2+dd250)(9/4)^2+dd150(27/4))/.rep50/.bb->1/2,Assumptions->cv>1];UU50=2Log[cv]+ii1/cv;
aa50=Log[cv];r150=ii1/cv;r250=ii2/cv;mm50=3-3aa50;vv50=3/2-9UU50/4;TT50=(9/4)(4aa50^2+4aa50 r150+r250);BB50=(mm50^2-3mm50 r150+(9/4)r250)/2;
ids50=<|"Phi2_first_derivative"->TrueQ[FullSimplify[D[h(1+Log[h]),h]-(2+Log[h])]==0],"Phi2_second_derivative"->TrueQ[FullSimplify[D[h(1+Log[h]),{h,2}]-1/h]==0],"Schur_single_constraint_sign"->TrueQ[FullSimplify[(-lam(TT+BB)-(-MM-BC)^2/(-B00/lam))-lam(-TT-BB+(MM+BC)^2/B00)]==0],"theta_one_third_polynomial_formula"->TrueQ[FullSimplify[F50-(81/8)(2UU50^2-(8/3)UU50-ii2/cv+(ii1/cv)^2),Assumptions->cv>1]==0],"closed_T_B_Schur_matches_derivative"->TrueQ[FullSimplify[F50-3(-TT50-BB50+2vv50^2),Assumptions->cv>1]==0],"normalization_a0_power_C"->TrueQ[FullSimplify[(Exp[-Log[cv]/bb])^bb cv==1,Assumptions->{cv>1,bb>0}]]|>;
value50[th0_,wp_:110]:=Module[{th=N[th0,wp+60],be,ac=2+Log[2],cc},be=th/(1-th);cc=1+NIntegrate[Exp[a-2]a^-be,{a,2,ac},WorkingPrecision->wp+20,AccuracyGoal->wp-10,PrecisionGoal->wp-10];Exp[-Log[cc]/be]];
fd50=Table[Module[{rr=case50[th,100],true,vals,errs},true=rr["budget_second_derivative_analytic"];vals=Table[(value50[th+hh]-2value50[th]+value50[th-hh])/hh^2,{hh,{1/100,1/1000,1/10000}}];errs=Abs[vals-true];<|"theta"->ToString[th,InputForm],"steps"->{"1/100","1/1000","1/10000"},"finite_difference_errors"->(ToString[#,InputForm]&/@errs),"last_error_lt_1e-4"->TrueQ[Last[errs]<10^-4],"second_order_shrink_at_least_50"->TrueQ[And@@Thread[Most[errs]/Rest[errs]>50]]|>],{th,{1/3,1/2}}];
lo50=case50[1/3,80];hi50=case50[1/3,160];prec50=Abs[lo50["budget_second_derivative_full_Hessian"]-hi50["budget_second_derivative_full_Hessian"]];
certPath50=If[Length[$ScriptCommandLine]>=3,$ScriptCommandLine[[3]],FileNameJoin[{DirectoryName[DirectoryName[$InputFileName]],"rounds","r050","certificate50.json"}]];cert50=Import[certPath50,"RawJSON"]["certificate"];
bounds50[key_]:=ToExpression/@Lookup[cert50[key],{"lower","upper"}];
nb50=Association@Table[key->bounds50[key],{key,{"actual_P_second","omit_kernel_prediction","omit_soft_prediction"}}];
certCheck50=<|"exact_rational_signs"->TrueQ[nb50["actual_P_second"][[2]]<0&&nb50["omit_kernel_prediction"][[1]]>0&&nb50["omit_soft_prediction"][[1]]>0]|>;
(* Match each rational enclosure with an independently evaluated value. *)
pairs50={{"actual_P_second","budget_second_derivative_analytic"},{"omit_kernel_prediction","omit_kernel_T_prediction"},{"omit_soft_prediction","omit_soft_prediction"}};
certCheck50["high_precision_values_inside_rational_bounds"]=And@@Table[With[{bd=nb50[pair[[1]]],vv=hi50[pair[[2]]]},TrueQ[bd[[1]]<vv<bd[[2]]]],{pair,pairs50}];
resV50=<|"wolfram"->$Version,"exact_symbolic_checks"->ids50,"finite_differences"->fd50,"precision80_to_160_stored30digit_difference"->ToString[prec50,InputForm],"precision_difference_lt_1e-25"->TrueQ[prec50<10^-25],"certificate_rational_checks"->certCheck50,"scope"->"exact symbolic identities, finite-difference and precision checks; independent native reading of rational sign bounds"|>;
outV50=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/verify50.json"];Export[outV50,resV50,"RawJSON"];Print[resV50];
