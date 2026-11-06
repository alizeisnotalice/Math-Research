base=DirectoryName[DirectoryName[$InputFileName]];
Clear[s];
posteriorC[tt_,cc_,qq_]:=Module[{weights,inds},inds=Select[Range[0,tt],qq-cc #>0&];weights=(Binomial[tt,#](qq-cc #)^(tt-1)&/@inds);Total[weights*((((qq-cc #)^2/(tt+1)+cc^2 # (tt-#)/(tt-1))&/@inds))]/Total[weights]];
exactCase[tt_,qq_,aa_]:=Module[{pt,conv,ct,formula},pt=Sum[Binomial[tt,k](qq-k)^(tt-1)Boole[qq>k],{k,0,tt}]/Factorial[tt-1];
conv=Total[Table[Sum[If[qq>k+l,Binomial[j,k]Binomial[tt-j,l]Integrate[(s-k)^(j-1)(qq-s-l)^(tt-j-1)(s-j aa)^2/j^2,{s,k,qq-l}]/(Factorial[j-1]Factorial[tt-j-1]),0],{k,0,j},{l,0,tt-j}],{j,1,tt-1}]]/pt+(qq-tt aa)^2/tt^2;
ct=posteriorC[tt,1,qq];formula=tt(qq/tt-aa)^2+ct/tt(HarmonicNumber[tt]-1);
<|"t"->tt,"q"->ToString[qq,InputForm],"center"->ToString[aa,InputForm],"direct_convolution_over_density_divided_by_p_over_8"->ToString[conv,InputForm],"posterior_formula"->ToString[formula,InputForm],"exact_residual"->Together[conv-formula]|>];
rows={exactCase[3,7/2,5/4],exactCase[4,5,5/4]};
Clear[j,t,a,v];harmonicResidual=FullSimplify[Sum[(t-j)/(j t^2),{j,1,t}]-(HarmonicNumber[t]-1)/t,Assumptions->Element[t,Integers]&&t>=1];
result=<|"wolfram"->$Version,"status"->"exact rational convolution checks with c=1; center arbitrary, not an endpoint certificate","rows"->rows,"harmonic_sum_residual"->harmonicResidual|>;

(* Target-matched bias cancellation, representative Gaussian thresholds. *)
mu=1+ProductLog[Exp[-1]];
matchRows=Table[Module[{z,q,cc,ct,pc,gg},z=zv/.FindRoot[Erfc[zv/Sqrt[2]]/2-Sqrt[mu]zv/Sqrt[nn]==0,{zv,1},WorkingPrecision->50,AccuracyGoal->35,PrecisionGoal->35];q=nn mu-Sqrt[nn mu]z;cc=cv/.FindRoot[2cv-(1+cv/(1+Exp[cv]))==q/nn,{cv,N[mu,50]},WorkingPrecision->50,AccuracyGoal->35,PrecisionGoal->35];ct=posteriorC[nn,cc,q];pc=1/(1+Exp[cc]);gg=pc ct/(8nn);<|"N"->nn,"q_over_N"->N[q/nn,18],"c"->N[cc,18],"bias_equation_residual"->ToString[N[2cc-(1+cc/(1+Exp[cc]))-q/nn,25],InputForm],"C_over_N"->N[ct/nn,18],"source_log_coefficient"->N[gg,18],"Poisson_shift"->N[gg(HarmonicNumber[nn]-1),18]|>],{nn,{2,4,8,16,32,64,128,256,512}}];
result=Join[result,<|"target_matched_rows"->matchRows,"limiting_source_log_coefficient"->N[(mu-1)/8,25]|>];

bridgeResult=result;
Clear[t,s,n,mu];
inv2[u_]:=Exp[u]/ProductLog[Exp[u+1]];
invAsym=Limit[Exp[-t]t inv2[t-s],t->Infinity,Assumptions->Element[s,Reals]];
(* Fenchel coefficient with representative critical scale and subexponential epsilon. *)
coefLog=Log[n]+(n Log[2]-LogGamma[n]-(mu n-Sqrt[n Log[n]])+Log[n]^2)/(n-1);
coefLimit=Limit[coefLog,n->Infinity,Assumptions->mu>0];
mu0=1+ProductLog[Exp[-1]];
result=<|"wolfram"->$Version,"inverse_envelope_Gamma2_integrand_limit"->ToString[invAsym,InputForm],"Fenchel_log_coefficient_limit"->ToString[coefLimit,InputForm],"Fenchel_coefficient_limit_at_mu_50digits"->ToString[N[2Exp[1-mu0],50],InputForm],"fixed_N"->2,"s_for_numeric_check"->3,"Gamma2_inverse_integrand_samples"->Table[<|"t"->tt,"value_40digits"->ToString[N[Exp[-tt]tt inv2[tt-3],40],InputForm]|>,{tt,{20,100,1000}}]|>;


tailResult=result;
Clear[c,b,x];
smallPref=(1+Exp[-c])^-b;tailPref=(2/(1+Exp[-c]))^b;
resSmall=FullSimplify[D[smallPref,c]-b/(1+Exp[c])smallPref,Assumptions->{b>=1,c>0}];
resTail=FullSimplify[D[tailPref,c]-b/(1+Exp[c])tailPref,Assumptions->{b>=1,c>0}];
surv[cc_,xx_]:=If[xx<cc,(Exp[-cc]+Exp[-xx])/(1+Exp[-cc]),2Exp[-xx]/(1+Exp[-cc])];
rows=Table[<|"z"->zz,"survival_c1"->N[surv[1,zz],30],"survival_c2"->N[surv[6/5,zz],30],"difference_c2_minus_c1"->N[surv[6/5,zz]-surv[1,zz],30]|>,{zz,{1/2,2}}];
result=<|"wolfram"->$Version,"small_CDF_prefactor_derivative_residual"->resSmall,"large_tail_prefactor_derivative_residual"->resTail,"single_step_survival_crossing"->rows|>;


orderResult=result;

s=Import[FileNameJoin[{base,"rounds","r030","receiver_seed.json"}],"RawJSON"];r=Import[FileNameJoin[{base,"rounds","r031","checks.json"}],"RawJSON"]["receiver_repair"];
w=r["uniform_scale_with_extra_margin"] s["weights"];av=r["final_N256_inner64"]["averages_by_group"];co=Join[{1},ConstantArray[2,Length[w]-1]];q=1+10^-5;
slack=co w(av-q);rows=MapThread[<|"group"->#1,"v"->#2,"weight"->#3,"average"->#4,"complementarity_slack"->#5|>&,{Range[Length[w]],s["v_nonnegative"],w,av,slack}];
js=r["final_N256_inner64"]["scaled_cost"];gap=Total[slack];result=<|"wolfram"->$Version,"status"->"floating continuous KKT diagnostic; no certificate","primal_scaled"->js,"dual_scaled"->js-gap,"relative_gap"->gap/js,"largest_active_slacks"->TakeLargestBy[Select[rows,#["weight"]>10^-10&],#["complementarity_slack"]&,8],"active_count"->Length[Select[w,#>10^-10&]]|>;


kktResult=result;

result=<|"bridge"->bridgeResult,"tail"->tailResult,"stochastic_order"->orderResult,"fixed_seed_KKT"->kktResult|>;Export[FileNameJoin[{base,"rounds","r032","proof_checks.json"}],result,"RawJSON"];Print["r032 proof and diagnostic checks completed"];
