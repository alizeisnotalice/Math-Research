(* Exact finite checks for the median and aspect-ratio bounds. *)
args=Rest[$ScriptCommandLine]; out=First[args];
Clear[x,k,j];
residues=Table[FullSimplify[Floor[(8k+j)/2]-k+1-3(8k+j)/8,k\[Element]Integers],{j,0,7}];
medianOK=And@@Thread[residues>=0];
muLower=1278/1000;
(* e^(1278/1000)<1/(muLower-1) follows by summing with a geometric upper tail. *)
expUpper=Sum[muLower^j/j!,{j,0,12}]+muLower^13/13!/(1-muLower/14);
muOK=expUpper (muLower-1)<1;
(* tanh-series bounds with a positive remainder. *)
logBounds[t_]:=Module[{u=(t-1)/(t+1),v},v=Sum[2u^(2j+1)/(2j+1),{j,0,99}];{v,v+2u^201/(201(1-u^2))}];
eLower=Sum[1/j!,{j,0,9}];
rateUpper=2 Last[logBounds[2]]-First[logBounds[1+15 eLower^2/4]]/30;
rateOK=rateUpper<1275/1000<muLower;
marginOK=Last[logBounds[8/3]]<muLower;
(* Use an independent positive symbol for u^n when checking the rational inequality. *)
scalePoly=Resolve[ForAll[{b,t,w,v},0<b<=t && w>=1 && v>(t/b)^2 w \[Implies] t/(b v)<b/(t w)],Reals];
pass=TrueQ[medianOK&&muOK&&rateOK&&marginOK&&scalePoly];
Export[out,<|"pass_all"->pass,"median_residues"->(ToString[#,InputForm]&/@residues),"mu_lower"->"1278/1000","block_entropy_upper_less_than"->"1275/1000","exact_mu_bracket_pass"->TrueQ[muOK],"exact_entropy_pass"->TrueQ[rateOK],"selection_penalty_pass"->TrueQ[marginOK],"aspect_interval_separation"->TrueQ[scalePoly],"kernel"->$Version|>,"RawJSON"];Print[pass];If[!pass,Exit[3]];
