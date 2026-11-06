(* Exact series coefficients and Doob identities, plus high-precision rate
   reconnaissance. No decimal scan is an interval certificate. *)
Clear[c,p,u,x,theta,nn,dd,kk,z];
op[f_]:=D[f,c]-p(1-p)D[f,p];critical[f_]:=Factor[f/.{c->u,p->(u-1)/u}];
m=1+c p;var=1+c^2 p(1-p);third=2+c^3 p(1-p)(1-2p);
A=u(u-1)/8;
Blambda=Factor[critical[Nest[op,m,3]](u/2)^3/6+critical[Nest[op,var,2]](u/2)^2/4+critical[op[third]](u/2)/6];
Bgap=Factor[Blambda-9A^2/(2u)];BD=Factor[A(2-u)];K4=Factor[Bgap-BD];
leftLogD2=FullSimplify[D[Log[1-x/nn],{x,2}],Assumptions->{nn>=1,x<0}];
moment=Integrate[2z^2 Exp[-(1-dd-2kk)z],{z,0,Infinity},Assumptions->{dd>=0,kk>=0,dd+2kk<1}];
binomialResiduals=Table[Table[Binomial[tt,j]j(tt-j)-tt(tt-1)Binomial[tt-2,j-1],{j,1,tt-1}],{tt,2,12}];
exact=<|"A"->ToString[A,InputForm],"B_lambda"->ToString[Blambda,InputForm],"B_gap"->ToString[Bgap,InputForm],"B_D"->ToString[BD,InputForm],
 "K_alpha1_quartic"->ToString[K4,InputForm],"K_positive_for_mu_gt1"->FullSimplify[K4>0,Assumptions->u>1],
 "phi_left_log_second_derivative"->ToString[leftLogD2,InputForm],"phi_log_slope_jump_allowed"->FullSimplify[-1/nn>=-1,Assumptions->nn>=1],
 "Doob_weighted_second_moment_upper"->ToString[moment,InputForm],"binomial_auxiliary_density_identity_residuals"->Flatten[binomialResiduals]|>;
mu=1+ProductLog[Exp[-1]];mp=N[mu,80];pp=1/(1+Exp[mp]);
lm[c_,th_]:=Log[(1+Exp[-(1-th)c])/((1-th)(1+Exp[-c]))];
cs[th_?NumericQ]:=cs[th]=cv/.FindRoot[Exp[(1-th)cv]-(1-th)Exp[cv]+th==0,{cv,mp+mp th/2},WorkingPrecision->60,AccuracyGoal->42,PrecisionGoal->42];
lst[th_?NumericQ]:=lm[cs[th],th];
ys[th_?NumericQ]:=1/(1-th)+cs[th]/(1+Exp[(1-th)cs[th]]);
ym[x_]:=1/(1-x)+mp/(1+Exp[(1-x)mp]);dexpr=pp/8(ym[x]-mp)^2 x(1-x);
d1=D[dexpr,x];d2=D[dexpr,{x,2}];y1=D[ym[x],x];y2=D[ym[x],{x,2}];
rateRow[eta_,al_]:=Module[{e=N[eta,70],yy,iv,dv,de,dee,ve,vee,jp,jpp,h,lo=0,hi,th,mid,is,gap},
 yy=ym[e];iv=e yy-lm[mp,e];dv=dexpr/.x->e;de=d1/.x->e;dee=d2/.x->e;ve=y1/.x->e;vee=y2/.x->e;
 jp=e-al de/ve;jpp=(1-al(dee ve-de vee)/ve^2)/ve;h=iv-al dv-yy jp+lst[jp];
 hi=N[1-10^-20,70];Do[mid=(lo+hi)/2;If[ys[mid]<yy,lo=mid,hi=mid],{85}];th=(lo+hi)/2;is=th yy-lst[th];gap=iv-is;
 <|"eta"->ToString[eta,InputForm],"alpha"->ToString[al,InputForm],"Jprime"->ToString[N[jp,35],InputForm],"Jsecond"->ToString[N[jpp,35],InputForm],
 "H_one_step_rate"->ToString[N[h,35],InputForm],"rate_gap_over_D"->ToString[N[gap/dv,35],InputForm]|>];
rows=Flatten[Table[rateRow[e,al],{e,{1/10000,1/10,2/5,9/10,99/100}},{al,{1,101/100,21/20}}]];
out=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],FileNameJoin[{$TemporaryDirectory,"rate_checks36.json"}]];
Export[out,<|"wolfram"->$Version,"status"->"exact algebra plus 70-digit inputs and 60-digit root computations; no interval certificate","exact"->exact,"rows"->rows|>,"RawJSON"];
Print[<|"output"->out,"K4"->K4,"auxiliary_residual_max"->Max[Abs[Flatten[binomialResiduals]]],"rows"->Length[rows]|>];
