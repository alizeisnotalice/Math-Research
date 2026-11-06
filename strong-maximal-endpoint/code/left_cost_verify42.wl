$left42LibraryMode=True;Get[FileNameJoin[{DirectoryName[$InputFileName],"left_cost_probe42.wl"}]];$left42LibraryMode=False;
argsV42=$ScriptCommandLine;outV42=If[Length[argsV42]>=2,argsV42[[2]],"/tmp/left_cost_verify42.json"];
fracRowsV42={};Do[ll=N[(1/(aa mu42(mu42-1)/8 tt))^(1/3),60];x=case42[tt,ll,aa,64,1/16];AppendTo[fracRowsV42,x];Print[KeyTake[x,{"t","alpha","residual_scaled_lambda2","min_action","seconds"}]],{tt,{512,8192}},{aa,{101/100,11/10}}];
contactV42={};Do[ll=N[(1/(aa mu42(mu42-1)/8 tt))^(1/3),60];qfun[l_?NumericQ]:=Module[{dj=setup42[tt,l,1/8]},aa dj[[2]]/dj[[1]]];lc=l/.FindRoot[qfun[l]==999/1000,{l,ll},WorkingPrecision->40,AccuracyGoal->28,PrecisionGoal->28];x=case42[tt,lc,aa,64,1/16];AppendTo[contactV42,x];Print[KeyTake[x,{"t","alpha","phase_fraction_alpha_J_over_D","residual_scaled_lambda2","min_action","seconds"}]],{tt,{512,8192,32768}},{aa,{101/100,11/10}}];
(* Independent finite binomial--Gamma mixture at the 512 time point. *)
tt=512;aa=101/100;ll=N[(1/(aa mu42(mu42-1)/8 tt))^(1/3),60];setup42[tt,ll,1/16];rr=rr42;ks=Range[0,tt];xx=Clip[rr-mu42 ks,{0,Infinity}];bi=N[Binomial[tt,#]&/@ks,65];ww=bi Exp[-mu42 ks]/(1+Exp[-mu42])^tt;
fExact=Exp[-rr]Dot[bi,xx^(tt-1)]/(Gamma[tt](1+Exp[-mu42])^tt);dExact=Dot[ww,GammaRegularized[tt,0,xx]]-Exp[-rr]Dot[bi,xx^tt]/(Gamma[tt+1](1+Exp[-mu42])^tt);
indV42=<|"t"->tt,"lambda"->N[ll,18],"density_relative_difference"->ToString[(fc42-fExact)/fExact,InputForm],"left_cost_over_density_difference"->ToString[dn42-dExact/fExact,InputForm]|>;
Export[outV42,<|"wolfram"->$Version,"status"->"quadrature refinement, exact phase-fraction roots and finite-mixture check; floating not certified","precision"->wp42,"refinements"->fracRowsV42,"near_contact"->contactV42,"independent_check"->indV42|>,"RawJSON"];Print[indV42];
