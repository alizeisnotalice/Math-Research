(* Reference Mills derivatives via continuous-time tilted Fourier extension.
   Numerical diagnostics, not interval certificates. *)
mu41=N[1+ProductLog[1/E],80];wp41=50;args41=$ScriptCommandLine;mode41=If[Length[args41]>=3,args41[[3]],"infinite"];
hn41[z_,nt_]:=1/(z(1-z))+If[nt===Infinity,0,1/(nt z^2)];
mf41[a_]:=1+a/(1+Exp[a]);
matchedA41[nt_]:=matchedA41[nt]=Module[{zz,target},zz=z/.FindRoot[Erfc[z/Sqrt[2]]/2==Sqrt[mu41]z/Sqrt[nt],{z,2},WorkingPrecision->60,AccuracyGoal->45,PrecisionGoal->45];target=mu41-Sqrt[mu41/nt]zz;a/.FindRoot[2a-mf41[a]==target,{a,mu41-.05},WorkingPrecision->60,AccuracyGoal->45,PrecisionGoal->45]];
y41[a_,e_]:=1/(1-e)+a Exp[(e-1)a]/(1+Exp[(e-1)a]);
Clear[delta41];
delta41[t_?NumericQ,e_?NumericQ,a_?NumericQ]:=delta41[t,e,a,24,Infinity];
delta41[t_?NumericQ,e_?NumericQ,a_?NumericQ,cut_]:=delta41[t,e,a,cut,Infinity];
delta41[t0_?NumericQ,e0_?NumericQ,a0_?NumericQ,cut_,nt_]:=delta41[t0,e0,a0,cut,nt]=Module[{tt=N[t0,70],ee=N[e0,70],aa=N[a0,70],yy,logf,nf,df,nv,dv},
 yy=y41[aa,ee];
 logf[v_?NumericQ]:=tt(Log[(1-ee)/(1-ee-I v/Sqrt[tt])]+Log[(1+Exp[(ee-1)aa+I aa v/Sqrt[tt]])/(1+Exp[(ee-1)aa])]-I yy v/Sqrt[tt]);
 df[v_?NumericQ]:=Re[Exp[logf[v]]];
 nf[v_?NumericQ]:=Re[Exp[logf[v]](hn41[ee+I v/Sqrt[tt],nt]/hn41[ee,nt]-1)];
 dv=NIntegrate[df[v],{v,0,cut},WorkingPrecision->wp41,AccuracyGoal->32,PrecisionGoal->28,MaxRecursion->18,Method->{"GlobalAdaptive","SymbolicProcessing"->0}];
 nv=NIntegrate[nf[v],{v,0,cut},WorkingPrecision->wp41,AccuracyGoal->32,PrecisionGoal->28,MaxRecursion->18,Method->{"GlobalAdaptive","SymbolicProcessing"->0}];nv/dv];
ks41=Range[-3,3];weights41[j_]:=weights41[j]=LinearSolve[Table[If[m==0,ConstantArray[1,Length[ks41]],ks41^m],{m,0,6}],Table[If[m==j,j!,0],{m,0,6}]];
fd41[f_,x_,h_,j_]:=Dot[weights41[j],f/@(x+h ks41)]/h^j;
mix41[tt_Integer,ee_,aa_,nt_]:=Module[{rr=tt y41[aa,ee],kk=Range[0,tt],xx,bb,ff,rs,qq,ww,call},xx=Clip[rr-aa kk,{0,Infinity}];bb=N[Binomial[tt,#]&/@kk,60];ww=bb Exp[-aa kk]/(1+Exp[-aa])^tt;
 ff=Exp[-rr]Dot[bb,xx^(tt-1)]/(Gamma[tt](1+Exp[-aa])^tt);rs=Exp[-rr]Dot[bb,xx^tt]/(Gamma[tt+1](1+Exp[-aa])^tt);qq=Dot[ww,GammaRegularized[tt,xx]];call=Dot[ww,tt GammaRegularized[tt+1,xx]-(rr-aa kk)GammaRegularized[tt,xx]];(rs+qq+If[nt===Infinity,0,call/nt])/(ff hn41[ee,nt])-1];
row41[tt_,A_]:=Module[{ee=N[A Sqrt[Log[tt]/tt],60],dd,d1,d2,dt1,dt2,eta1,eta2,tail,ac,bound,im,aa,nt},
 nt=If[mode41==="matched",tt,Infinity];aa=If[mode41==="matched",matchedA41[tt],mu41];dd=delta41[tt,ee,aa,24,nt];eta1=Table[fd41[delta41[tt,#,aa,24,nt]&,ee,ee/100,j],{j,1,3}];eta2=Table[fd41[delta41[tt,#,aa,24,nt]&,ee,ee/200,j],{j,1,3}];
 dt1=fd41[delta41[#,ee,aa,24,nt]&,tt,tt/100,1];dt2=fd41[delta41[#,ee,aa,24,nt]&,tt,tt/200,1];
 ac=(1-ee)^2 tt;bound=24(ac+24^2)/( (tt-1)24^2-ac)(1+24^2/ac)^(-tt/2);
 im=If[tt<=1024,mix41[tt,ee,aa,nt],delta41[tt,ee,aa,32,nt]];
 <|"t"->tt,"horizon_N"->If[nt===Infinity,"infinity",nt],"a"->N[aa,18],"A"->A,"eta"->N[ee,18],"delta"->N[dd,24],"scaled_delta_t_eta2"->N[tt ee^2 dd,20],"scaled_time_derivative"->N[tt^2 ee^2 dt2,20],
 "scaled_eta_derivatives"->N[Table[tt ee^(2+j)eta2[[j]],{j,1,3}],20],"FD_time_step_halving_absolute_difference"->ToString[N[Abs[dt2-dt1],8],InputForm],
 "FD_eta_scaled_halving_differences"->N[Table[tt ee^(2+j)Abs[eta2[[j]]-eta1[[j]]],{j,1,3}],16],"FD_eta_step_halving_absolute_differences"->(ToString[N[#,8],InputForm]&/@Abs[eta2-eta1]),"independent_check"->If[tt<=1024,"exact finite Gamma-binomial mixture","Fourier cutoff raised 24 to 32"],"independent_absolute_difference"->ToString[N[Abs[im-dd],8],InputForm],
 "Fourier_denominator_tail_absolute_bound"->N[bound,18]|>];
If[!TrueQ[$mills41LibraryMode],
args41=$ScriptCommandLine;out41=If[Length[args41]>=2,args41[[2]],"/tmp/mills_derivative_probe41.json"];
rows41=Flatten[Table[row41[t,A],{t,{256,1024,4096}},{A,{1,2}}],1];
Export[out41,<|"wolfram"->$Version,"precision"->wp41,"mode"->mode41,"integration_accuracy_goal"->32,"integration_precision_goal"->28,"Fourier_cutoff"->24,"FD_stencil"->ks41,"FD_moment_residuals"->Table[Table[Dot[weights41[j],If[m==0,ConstantArray[1,Length[ks41]],ks41^m]]-If[m==j,j!,0],{m,0,6}],{j,1,3}],"rows"->rows41,"status"->"finite differences + high-precision integrals; no interval guarantee"|>,"RawJSON"];
Print[<|"output"->out41,"cases"->Length[rows41]|>];

];
