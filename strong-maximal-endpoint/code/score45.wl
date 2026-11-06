(* Global score minimum: exact truncated-power sums at finite b, tilted
   center Fourier at large b. Floating checks do not certify the supremum. *)
mu45=N[1+ProductLog[1/E],180];p45=1/(1+Exp[mu45]);wp45=90;
u45[c_]:=c/(1+Exp[c]);m45[c_]:=1+u45[c];nu45[c_]:=1+c^2 Exp[c]/(1+Exp[c])^2;
third45[c_]:=2+c^3 Exp[c](Exp[c]-1)/(1+Exp[c])^3;
setupP45[b_Integer]:=(bin45=N[Binomial[b,#]&/@Range[0,b],wp45+20];k45=Range[0,b];);
powers45[b_,x_?NumericQ]:=Module[{z=Clip[x-k45,{0,Infinity}]},Table[If[b-1-j>=0,FactorialPower[b-1,j]Dot[bin45,z^(b-1-j)],0],{j,0,2}]];
dVal45[b_,x_?NumericQ]:=Module[{v=powers45[b,x]},x v[[2]]/v[[1]]-(b-1)];
dDer45[b_,x_?NumericQ]:=Module[{v=powers45[b,x]},v[[2]]/v[[1]]+x(v[[3]]/v[[1]]-(v[[2]]/v[[1]])^2)];
finite45[b_Integer]:=Module[{start=AbsoluteTime[],xs,vals,ii,lo,hi,xm,opt,dmax,gridmax},setupP45[b];xs=N[Subdivide[1/100,5b,256],wp45];vals=dVal45[b,#]&/@xs;gridmax=Max[vals];lo=N[b-2,wp45];hi=N[b+1,wp45];
 opt=x/.FindRoot[dDer45[b,x]==0,{x,N[b-(3-mu45)/(2mu45),wp45]},WorkingPrecision->70,AccuracyGoal->45,PrecisionGoal->45];dmax=dVal45[b,opt];
 <|"b_exact"->ToString[b,InputForm],"method"->"exact finite truncated-power sums and stationary root; grid check only","x_peak"->N[opt,20],"x_minus_b"->N[opt-b,20],"D_candidate"->N[dmax,20],"D_minus_b_u_mu"->N[dmax-b u45[mu45],20],"peak_root_residual"->ToString[dDer45[b,opt],InputForm],"grid_max_le_candidate"->TrueQ[gridmax<=dmax],"mu_score_inf_candidate"->N[b p45-dmax/mu45,20],"seconds"->N[AbsoluteTime[]-start]|>];
Clear[bridge45];bridge45[b_,d0_?NumericQ,step_:1/8,cut_:24]:=Module[{dd=N[d0,wp45+30],vv,ps0,ps1,ms,xx,fac},vv=N[Range[-cut/step,cut/step]step,wp45];ms=m45[dd];xx=ms-dd;
 ps0=Exp[b(-Log[1-I vv/Sqrt[b]]+Log[(1+Exp[-dd+I dd vv/Sqrt[b]])/(1+Exp[-dd])])-I ms vv Sqrt[b]];
 ps1=Exp[(b-1)(-Log[1-I vv/Sqrt[b-1]]+Log[(1+Exp[-dd+I dd vv/Sqrt[b-1]])/(1+Exp[-dd])])-I ms vv Sqrt[b-1]-I vv xx/Sqrt[b-1]];
 Sqrt[b/(b-1)]Re[Total[ps1]]/Re[Total[ps0]]];
large45[b_]:=Module[{start=AbsoluteTime[],d0,bv,k,dv,score},wp45=Max[100,2IntegerLength[b,10]+80];mu45=N[1+ProductLog[1/E],wp45+80];p45=1/(1+Exp[mu45]);d0=mu45+(3-mu45)/(2b);bv=bridge45[b,d0,1/8,If[b>=2^256,32,24]];k=u45[d0]bv;dv=b k;
 <|"b_exact"->ToString[b,InputForm],"method"->"Fourier evaluation at asymptotic maximizing d; not optimizer certificate","working_precision"->wp45,"d_minus_mu_times_b"->N[b(d0-mu45),20],"D_minus_b_u_mu"->N[dv-b u45[mu45],20],"bridge_correction_b_times_B_minus1"->N[b(bv-1),20],"mu_score_value"->N[b p45-dv/mu45,20],"local_scale_scores"->Table[Module[{cc=mu45+kk/Sqrt[b]},<|"kappa"->kk,"score_value"->N[b/(1+Exp[cc])-dv/cc,20],"predicted_limit"->N[-p45/2-p45 kk^2/(2mu45),20]|>],{kk,{-2,-1,0,1,2}}],"seconds"->N[AbsoluteTime[]-start]|>];
If[!TrueQ[$score45LibraryMode],args45=$ScriptCommandLine;out45=If[Length[args45]>=2,args45[[2]],"/tmp/score45.json"];fr45={};lr45={};Do[row=finite45[b];AppendTo[fr45,row];Print[KeyTake[row,{"b_exact","x_minus_b","D_minus_b_u_mu","mu_score_inf_candidate","seconds"}]],{b,{4,8,16,32,64,128,256,1024}}];
 Do[row=large45[b];AppendTo[lr45,row];Print[KeyTake[row,{"b_exact","D_minus_b_u_mu","mu_score_value","seconds"}]],{b,{4096,65536,2^24,2^48,2^96,2^192,2^384}}];
 Export[out45,<|"wolfram"->$Version,"status"->"score-minimum diagnostics; exact formulas but floating stationary roots/Fourier, not interval-certified","finite_rows"->fr45,"large_rows"->lr45,"limits"-><|"D_minus_b_u_mu"->N[(mu45-1)/2,20],"mu_score_min"->N[-p45/2,20],"d_shift_times_b"->N[(3-mu45)/2,20],"x_shift"->N[-(3-mu45)/(2mu45),20]|>|>,"RawJSON"];
];
