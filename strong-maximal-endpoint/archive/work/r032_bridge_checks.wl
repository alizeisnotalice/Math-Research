Clear[s];
posteriorC[tt_,cc_,qq_]:=Module[{weights,inds},inds=Select[Range[0,tt],qq-cc #>0&];weights=(Binomial[tt,#](qq-cc #)^(tt-1)&/@inds);Total[weights*((((qq-cc #)^2/(tt+1)+cc^2 # (tt-#)/(tt-1))&/@inds))]/Total[weights]];
exactCase[tt_,qq_,aa_]:=Module[{pt,conv,ct,formula},pt=Sum[Binomial[tt,k](qq-k)^(tt-1)Boole[qq>k],{k,0,tt}]/Factorial[tt-1];
conv=Total[Table[Sum[If[qq>k+l,Binomial[j,k]Binomial[tt-j,l]Integrate[(s-k)^(j-1)(qq-s-l)^(tt-j-1)(s-j aa)^2/j^2,{s,k,qq-l}]/(Factorial[j-1]Factorial[tt-j-1]),0],{k,0,j},{l,0,tt-j}],{j,1,tt-1}]]/pt+(qq-tt aa)^2/tt^2;
ct=posteriorC[tt,1,qq];formula=tt(qq/tt-aa)^2+ct/tt(HarmonicNumber[tt]-1);
<|"t"->tt,"q"->ToString[qq,InputForm],"center"->ToString[aa,InputForm],"direct_convolution_over_density_divided_by_p_over_8"->ToString[conv,InputForm],"posterior_formula"->ToString[formula,InputForm],"exact_residual"->Together[conv-formula]|>];
rows={exactCase[3,7/2,5/4],exactCase[4,5,5/4]};
Clear[j,t,a,v];harmonicResidual=FullSimplify[Sum[(t-j)/(j t^2),{j,1,t}]-(HarmonicNumber[t]-1)/t,Assumptions->Element[t,Integers]&&t>=1];
result=<|"wolfram"->$Version,"status"->"exact rational convolution checks with c=1; center arbitrary, not an endpoint certificate","rows"->rows,"harmonic_sum_residual"->harmonicResidual|>;
Export["/Users/zhengzhihao/Desktop/T/work/r032_bridge_checks.json",result,"RawJSON"];Print[result];
(* Target-matched bias cancellation, representative Gaussian thresholds. *)
mu=1+ProductLog[Exp[-1]];
matchRows=Table[Module[{z,q,cc,ct,pc,gg},z=zv/.FindRoot[Erfc[zv/Sqrt[2]]/2-Sqrt[mu]zv/Sqrt[nn]==0,{zv,1},WorkingPrecision->50,AccuracyGoal->35,PrecisionGoal->35];q=nn mu-Sqrt[nn mu]z;cc=cv/.FindRoot[2cv-(1+cv/(1+Exp[cv]))==q/nn,{cv,N[mu,50]},WorkingPrecision->50,AccuracyGoal->35,PrecisionGoal->35];ct=posteriorC[nn,cc,q];pc=1/(1+Exp[cc]);gg=pc ct/(8nn);<|"N"->nn,"q_over_N"->N[q/nn,18],"c"->N[cc,18],"bias_equation_residual"->ToString[N[2cc-(1+cc/(1+Exp[cc]))-q/nn,25],InputForm],"C_over_N"->N[ct/nn,18],"source_log_coefficient"->N[gg,18],"Poisson_shift"->N[gg(HarmonicNumber[nn]-1),18]|>],{nn,{2,4,8,16,32,64,128,256,512}}];
result=Join[result,<|"target_matched_rows"->matchRows,"limiting_source_log_coefficient"->N[(mu-1)/8,25]|>];Export["/Users/zhengzhihao/Desktop/T/work/r032_bridge_checks.json",result,"RawJSON"];Print[matchRows];
