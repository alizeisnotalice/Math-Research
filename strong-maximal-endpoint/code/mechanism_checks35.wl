(* Exact mechanism identities and high-precision right-tail calibration.
   The rational certificates below are exact; decimal values are not intervals. *)
Clear[uu,xx,nn,jj,ss,cc0,theta];
alpha=15/16;mu=1+ProductLog[Exp[-1]];
q2=1/3;log2lo=2Sum[q2^(2jj+1)/(2jj+1),{jj,0,3}];
log2hi=log2lo+2q2^9/(9(1-q2^2));
q7=5/9;log7lo=2Sum[q7^(2jj+1)/(2jj+1),{jj,0,3}];
fexpr=(nn-1)Log[1+alpha xx]-(4-alpha)xx;
haloD=FullSimplify[D[Exp[4uu]-1-4Exp[uu]uu,uu]];
concave=FullSimplify[D[fexpr,{xx,2}]<0,Assumptions->{nn>=8,xx>=0}];
tensor=Table[Module[{aa=15/16,bb=1+Log[ss],expr},
 expr=Sum[Binomial[kk,jj]aa^jj(1-aa)^(kk-jj)(1+uu)^jj(bb+jj Log[1+uu])^(nn-1),{jj,0,kk}];
 FullSimplify[(D[expr,uu]/.uu->0)-bb^(nn-1)kk aa(1+(nn-1)/bb),Assumptions->{ss>=1,nn>=2}]],{kk,1,6}];
product=<|"alpha"->"15/16","log2_lower"->ToString[log2lo,InputForm],"log2_upper"->ToString[log2hi,InputForm],
 "log2_bracket_verifies_2over3_7over10"->TrueQ[log2lo>2/3&&log2hi<7/10],
 "log_7over2_lower_minus_5over4"->ToString[log7lo-5/4,InputForm],
 "endpoint_n8_rational_lower"->ToString[7(5/4)-8(343/320),InputForm],
 "n_derivative_rational_lower"->ToString[1+35/58-343/320,InputForm],
 "F_concavity"->concave,"halo_growth_derivative"->ToString[haloD,InputForm],
 "tensor_first_derivative_residuals"->(ToString[#,InputForm]&/@tensor)|>;
M[c_,th_]:=(1+Exp[-(1-th)c])/((1-th)(1+Exp[-c]));
ma[c_]:=1+c/(1+Exp[c]);
calibration=Table[Module[{ct,yy,ratio},
 ct=cv/.FindRoot[Exp[(1-th)cv]-(1-th)Exp[cv]+th==0,{cv,N[mu,80]},WorkingPrecision->80,AccuracyGoal->60,PrecisionGoal->60];
 yy=1/(1-th)+mu/(1+Exp[(1-th)mu]);ratio=M[ct,th]/M[mu,th];
 <|"eta"->ToString[th,InputForm],"c_star"->ToString[N[ct,50],InputForm],"y"->ToString[N[yy,50],InputForm],
 "Mstar_over_Mmu"->ToString[N[ratio,50],InputForm],
 "stationarity_residual"->ToString[N[Exp[(1-th)ct]-(1-th)Exp[ct]+th,50],InputForm]|>],{th,{1/20,1/4}}];
(* Exact critical-parameter identities after e^u=1/(u-1). *)
m2=FullSimplify[(D[ma[cc0],{cc0,2}]/.cc0->uu)/.{Exp[uu]->1/(uu-1),Exp[2uu]->1/(uu-1)^2,Exp[3uu]->1/(uu-1)^3}];
identity=FullSimplify[m2+(uu-1)/uu];
(* Fixed eta comparison: F(mu,eta)>0 follows from g(z)=e^(mu z)-z e^mu+1-z,
   g(1)=g'(1)=0 and g''>0.  Both identities are exact below. *)
gg=Exp[uu xx]-xx Exp[uu]+1-xx;
ggAt1=FullSimplify[gg/.xx->1];ggPrimeAt1=FullSimplify[(D[gg,xx]/.xx->1)/.Exp[uu]->1/(uu-1)];
cells={{{0,1},{1,2}},{{1,2},{2,3}},{{2,3},{1,2}},{{3,4},{2,3}}};
singleMax=Function[qq,16/((qq[[1,2]]+7/2)(qq[[2,2]]+7/2))]/@cells;
shape=Table[Boole[MemberQ[cells,{{ii-1,ii},{jj-1,jj}}]],{ii,1,4},{jj,1,4}];
geometry=<|"input_cells"->cells,"amplitude"->16,"full_mass"->64,"receiver_box"->{{"-4","-7/2"},{"-4","-7/2"}},
 "area_exact"->"1/4","full_witness"->"[x,4] x [y,3]","full_average_lower_exact"->ToString[64/56,InputForm],
 "single_component_maxima_exact"->(ToString[#,InputForm]&/@singleMax),"full_strict_margin"->"1/7",
 "single_component_uniform_gap"->ToString[1-Max[singleMax],InputForm],"shape_rank"->MatrixRank[shape],
 "checks_pass"->TrueQ[64/56>1&&And@@Thread[singleMax<1]&&MatrixRank[shape]==2]|>;
result=<|"wolfram"->$Version,"status"->"exact rational and symbolic checks plus 80-digit numerical calibration; decimal calibrations not interval certified",
 "product_obstruction"->product,"critical_m_second_derivative_residual"->ToString[identity,InputForm],
 "mgf_root_comparison_g_at1"->ToString[ggAt1,InputForm],"mgf_root_comparison_gprime_at1"->ToString[ggPrimeAt1,InputForm],
 "rectangle_baseline_at_A4"->ToString[N[(1/4+8Log[2])/(1+2Log[2]),50],InputForm],
 "grid_best_witness_below_rectangle_baseline"->TrueQ[8log2lo-11/4>0],
 "right_LD_calibration"->calibration,"continuous_cross_component_certificate"->geometry|>;
out=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],FileNameJoin[{$TemporaryDirectory,"mechanism_checks35.json"}]];
Export[out,result,"RawJSON"];Print[<|"output"->out,"tensor_residuals"->tensor,"critical_m2_residual"->identity|>];
