(* r034: exact affine/Gamma supersolution and independent finite-mixture checks.
   Wolfram 15; deterministic 70-digit arithmetic, no interval certificate. *)
ClearAll[mu,rho,bfun,bclosed,bint,left,vref,phi];
mu=1+ProductLog[Exp[-1]];
rho[nn_]:=nn+ProductLog[-nn Exp[-nn]];
phi[nn_,rr_]:=If[rr>0,Exp[-rr],1-rr/nn];
bfun[nn_,rr_]:=If[rr>rho[nn],Exp[-rr]+rr/nn-1,0];
bclosed[nn_,tt_,rr_]:=Module[{xx=rr-rho[nn]},If[xx<=0,0,
 Exp[-rr]xx^tt/Gamma[tt+1]+(rr/nn-1)GammaRegularized[tt,0,xx]
 -tt/nn GammaRegularized[tt+1,0,xx]]];
bint[nn_,tt_,rr_]:=Module[{xx=N[rr-rho[nn],70]},If[xx<=0,0,
 NIntegrate[(Exp[-rr+y]+(rr-y)/nn-1)y^(tt-1)Exp[-y]/Gamma[tt],{y,0,xx},
 WorkingPrecision->70,AccuracyGoal->48,PrecisionGoal->45]]];
left[nn_,tt_,rr_]:=1-rr/nn+tt mu/nn+bclosed[nn,tt,rr];
vref[nn_,tt_,rr_]:=Module[{ks=Range[0,tt],x,ww,qq,call,rs},
 If[rr<=0,Return[1-rr/nn+tt mu/nn]];
 x=rr-mu ks;ww=N[Binomial[tt,#]&/@ks,70]Exp[-mu ks]/(1+Exp[-mu])^tt;
 qq=Total[ww GammaRegularized[tt,Clip[x,{0,Infinity}]]];
 call=Total[ww(tt GammaRegularized[tt+1,Clip[x,{0,Infinity}]]-x GammaRegularized[tt,Clip[x,{0,Infinity}]])];
 rs=Exp[-rr]Total[N[Binomial[tt,#]&/@ks,70]Clip[x,{0,Infinity}]^tt]/((1+Exp[-mu])^tt Gamma[tt+1]);
 rs+qq+call/nn];
ints=Flatten[Table[With[{rr=N[rho[nn]+xx,70]},<|"N"->nn,"t"->tt,"x"->xx,
 "closed"->ToString[N[bclosed[nn,tt,rr],40],InputForm],
 "independent_integral_residual"->ToString[N[bclosed[nn,tt,rr]-bint[nn,tt,rr],35],InputForm]|>],
 {nn,{2,16,64}},{tt,{1,4,16}},{xx,{1/2,4,32}}]];
rows=Flatten[Table[Module[{rr=N[tt mu+sgn 2Sqrt[tt Log[tt]],70],vv,ll},
 vv=vref[nn,tt,rr];ll=left[nn,tt,rr];<|"N"->nn,"t"->tt,"side"->sgn,
 "r"->N[rr,20],"reference_payoff"->N[vv,20],"left_barrier"->N[ll,20],
 "left_minus_reference"->ToString[N[ll-vv,30],InputForm],
 "Gamma_correction"->ToString[N[bclosed[nn,tt,rr],30],InputForm]|>],
 {nn,{16,64,256}},{tt,{8,16}},{sgn,{-1,1}}]];
diagonal=Flatten[Table[Module[{tt=nn,rr,vv,ll,th,cc,cb},
 rr=N[tt mu+sgn 2Sqrt[tt Log[tt]],70];vv=vref[nn,tt,rr];ll=left[nn,tt,rr];
 th=N[If[sgn<0,1/nn,2/mu Sqrt[Log[tt]/tt]],70];
 cc=cv/.FindRoot[Exp[(1-th)cv]-(1-th)Exp[cv]+th==0,{cv,N[mu,70]},WorkingPrecision->70,AccuracyGoal->48,PrecisionGoal->48];
 cb=Exp[-th rr]((1+Exp[-(1-th)cc])/((1-th)(1+Exp[-cc])))^tt;
 <|"N"->nn,"t"->tt,"side"->sgn,"reference_payoff"->N[vv,20],"left_barrier"->N[ll,20],
 "theta"->N[th,20],"positive_mgf_max_c"->N[cc,20],"exponential_barrier"->N[cb,20],
 "combined_barrier"->N[Min[ll,cb],20],"combined_minus_reference"->N[Min[ll,cb]-vv,20]|>],
 {nn,{16,64,256}},{sgn,{-1,1}}]];
Clear[uu,aaa,bbb,zz];
ff=Exp[(1-zz)(uu+aaa zz+bbb zz^2)]-(1-zz)Exp[uu+aaa zz+bbb zz^2]+zz;
implicitCoefficients=Table[FullSimplify[SeriesCoefficient[ff,{zz,0,j}]/.Exp[uu]->1/(uu-1)],{j,1,3}];
implicitResiduals=FullSimplify[implicitCoefficients/.{aaa->uu/2,bbb->uu(9-uu)/24}];
Clear[zz,rrr,hhh];
weightedCrossingResidual=FullSimplify[Integrate[Exp[-(rrr-zz)]Exp[-zz],{zz,hhh,rrr},
 Assumptions->0<hhh<rrr]-(rrr-hhh)Exp[-rrr],Assumptions->0<hhh<rrr];
crossings=Table[Module[{aa=1,dd,kk,tt,rr,lo,hi,pr,zn},
 zn=zv/.FindRoot[Erfc[zv/Sqrt[2]]/2==Sqrt[mu]zv/Sqrt[nn],{zv,N[Sqrt[Log[nn]],70]},WorkingPrecision->70,AccuracyGoal->48,PrecisionGoal->48];
 dd=N[Sqrt[mu nn]zn,70];kk=Floor[N[(dd+2aa Sqrt[nn Log[nn]])/(mu-1),70]];tt=nn-kk;
 rr=N[aa Sqrt[tt Log[tt]],70];lo=N[nn mu-dd-tt mu-3rr,70];hi=N[nn mu-dd-tt mu-rr,70];
 pr=GammaRegularized[kk,0,Max[0,lo]]+GammaRegularized[kk,hi];
 <|"N"->nn,"elapsed_k"->kk,"remaining_t"->tt,"A"->aa,"implicit_d_N"->N[dd,20],
 "event"->"A sqrt(t log t) <= q-Gamma_k-t mu <= 3 A sqrt(t log t)",
 "outside_probability"->ToString[N[pr,40],InputForm],"inside_probability"->ToString[N[1-pr,40],InputForm]|>],
 {nn,{65536,1048576,16777216}}];
(* Verify the convolution differential identity using exact integer shapes,
   keeping rho formal and substituting exp(-rho)=1-rho/N at the end. *)
Clear[rr,aa,nn];
bc[0]=Exp[-rr]+rr/nn-1;
Do[bc[tt]=Exp[-rr](rr-aa)^tt/Factorial[tt]+(rr/nn-1)(1-Exp[-(rr-aa)]Sum[(rr-aa)^j/Factorial[j],{j,0,tt-1}])
 -tt/nn (1-Exp[-(rr-aa)]Sum[(rr-aa)^j/Factorial[j],{j,0,tt}]),{tt,1,5}];
ode=Table[FullSimplify[D[bc[tt],rr]-(bc[tt-1]-bc[tt])],{tt,1,5}];
(* Here B_1'=b-B_1 also holds without the root relation: the upper endpoint
   b(rho)=0 is used only for monotonicity and nonnegativity. *)
result=<|"wolfram"->$Version,"working_precision"->70,"status"->"exact identities plus high-precision quadrature evidence; not interval certification",
 "convolution_ODE_residuals"->(ToString[#,InputForm]&/@ode),
 "mgf_implicit_coefficients"->(ToString[#,InputForm]&/@implicitCoefficients),
 "mgf_action_series_residuals"->(ToString[#,InputForm]&/@implicitResiduals),
 "mgf_action_series"->"c(theta)=mu+mu theta/2+mu (9-mu) theta^2/24+O(theta^3)",
 "weighted_crossing_integral_residual"->ToString[weightedCrossingResidual,InputForm],
 "closed_vs_integral"->ints,"tube_rows"->rows,"diagonal_two_screens"->diagonal,
 "pure_Exp_control_annulus_counterexample"->crossings|>;
out=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],FileNameJoin[{$TemporaryDirectory,"barrier_checks34.json"}]];
Export[out,result,"RawJSON"];Print[<|"output"->out,"ODE"->ode,"integral_rows"->Length[ints]|>];
