(* Actual clipped nonlinear cost; negative tilted Fourier quadrature. Floating evidence. *)
mu43=N[1+ProductLog[1/E],120];a43=mu43;pa43=1/(1+Exp[a43]);wp43=50;cut43=18;N43=Infinity;source43=0;kappa43=0;ell43=mu43;mean43=mu43;
lam43[e_]:=-Log[1-e]+Log[(1+Exp[(e-1)a43])/(1+Exp[-a43])];
y43[e_]:=1/(1-e)+a43 Exp[(e-1)a43]/(1+Exp[(e-1)a43]);
h43[e_]:=-1/(e(1-e))-If[N43===Infinity,0,1/(N43 e^2)];g143[e_]:=1/(1-e)^2;g243[e_]:=Exp[a43 e]/(1-pa43+pa43 Exp[a43 e])^2;
exprel43[x_]:=If[Abs[x]<10^-12,1+x/2+x^2/6+x^3/24+x^4/120,2Exp[x/2]Sinh[x/2]/x];
std43[n_]:=std43[n]=Module[{b,j,ev,vec,ix},b=N[Range[n-1]/Sqrt[4Range[n-1]^2-1],wp43];j=DiagonalMatrix[b,1]+DiagonalMatrix[b,-1];{ev,vec}=Eigensystem[j];ix=Ordering[ev];{ev[[ix]],2vec[[ix,1]]^2}];
quad43[n_,aa_,bb_]:=Module[{v=std43[n]}, {(bb-aa)v[[1]]/2+(aa+bb)/2,(bb-aa)v[[2]]/2}];
Clear[setup43];
setup43[tt_,lambda_,step_:1/8]:=Module[{e=-N[lambda,wp43+15],vv,ss,ps,cur,pre,n=tt,fac},
 e43=e;t43=tt;yy43=y43[e];rr43=tt yy43;vv=N[Range[-cut43/step,cut43/step]step,wp43];ss=e+I vv/Sqrt[tt];
 ps=Exp[tt(Log[(1-e)/(1-ss)]+Log[(1+Exp[(ss-1)a43])/(1+Exp[(e-1)a43])]-I yy43 vv/Sqrt[tt])];
 score43=-e+Re[Total[ps(-I vv/Sqrt[tt])]]/Re[Total[ps]];cur=Re[Total[ps]] step/(2Pi Sqrt[tt]);fc43=Exp[tt lam43[e]-e rr43]cur;
 dn43=Re[Total[ps(h43/@ss)]] step/(2Pi Sqrt[tt])/cur;
 cn43=Re[Total[ps((g143/@ss)+a43^2 pa43(1-pa43)(g243/@ss))]]step/(2Pi Sqrt[tt])/cur;
 jn43=pa43/8(tt(yy43-ell43)^2+(HarmonicNumber[tt]-1)cn43);
 (* same eta, previous time; normalized by current physical density *)
 vv43=vv;ss43=e+I vv/Sqrt[tt-1];pre=Exp[(tt-1)(Log[(1-e)/(1-ss43)]+Log[(1+Exp[(ss43-1)a43])/(1+Exp[(e-1)a43])]-I yy43 vv/Sqrt[tt-1])];
 pre43=pre step/(2Pi Sqrt[tt-1])/cur;
 dpre43=pre43(h43/@ss43);cpre43=pre43((g143/@ss43)+a43^2 pa43(1-pa43)(g243/@ss43));
 {dn43,jn43}];
raw43[z_?NumericQ,alpha_]:=Module[{ph,fr,dr,cr,jr,x},ph=Exp[-I vv43(yy43-z)/Sqrt[t43-1]];fr=Re[Total[pre43 ph]];dr=Re[Total[dpre43 ph]];cr=Re[Total[cpre43 ph]];jr=pa43/8(fr(rr43-z-(t43-1)ell43)^2/(t43-1)+(HarmonicNumber[t43-1]-1)cr);{fr,dr,jr}=Exp[-lam43[e43]+e43 z]{fr,dr,jr};x=fc43(alpha jr+source43 fr)/(1-(rr43-z)/N43+(t43-1)mean43/N43-fc43 dr);dr-(alpha jr+source43 fr)exprel43[x]-kappa43(t43-1)/(N43 fc43)];
rawCurrent43[alpha_]:=Module[{x=fc43(alpha jn43+source43)/(1-rr43/N43+t43 mean43/N43-fc43 dn43)},dn43-(alpha jn43+source43)exprel43[x]-kappa43 t43/(N43 fc43)];
case43[tt_,lam_,alpha_,ord_:32,step_:1/8,clipped_:True]:=Module[{d,j,front,lo,hi,mid,q,qn,zn,wn,b0,bb,cost,grid,best,opt,cur,left,right,der,ca,start=AbsoluteTime[]},
 {d,j}=setup43[tt,lam,step];cur=rawCurrent43[alpha];front=80;
 If[clipped,If[raw43[0,alpha]<=0,front=0,If[raw43[80,alpha]<0,lo=0;hi=80;Do[mid=(lo+hi)/2;If[raw43[mid,alpha]>0,lo=mid,hi=mid],{55}];front=(lo+hi)/2]]];
 bb[cc_?NumericQ]:=If[cc>=front,0,Module[{x,w},{x,w}=quad43[ord,cc,front];Dot[w,Exp[-x](If[clipped,Max[0,raw43[#,alpha]],raw43[#,alpha]]&/@x)]]];
 b0=bb[0];cost[cc_?NumericQ]:=(b0+bb[cc])/(1+Exp[-cc]);
 grid=Table[{N[cc],cost[N[cc,wp43]]},{cc,0,6,1/4}];best=First[MinimalBy[grid,Last]];
 der[cc_?NumericQ]:=cost[cc]-If[clipped,Max[0,raw43[cc,alpha]],raw43[cc,alpha]];left=N[0,wp43];right=N[6,wp43];If[der[left]<=0&&der[right]>=0,Do[mid=(left+right)/2;If[der[mid]<0,left=mid,right=mid],{55}];ca=(left+right)/2,ca=N[best[[1]],wp43]];opt={cost[ca],{cc->ca}};
 <|"N_exact"->ToString[N43,InputForm],"N"->N[N43,18],"matched_a"->N[a43,18],"b"->N[mean43-a43,18],"source"->N[source43],"kappa"->N[kappa43,18],"active"->TrueQ[cur>0],"floor_over_density"->N[kappa43 t43/(N43 fc43),18],"working_precision"->wp43,"clipped"->clipped,"t"->tt,"lambda"->N[lam,18],"alpha"->N[alpha],"phase_fraction_alpha_J_over_D"->N[alpha j/d,18],"phase_fraction_total_U_over_D"->N[(alpha j+source43)/d,18],"score"->N[score43,18],"current_cost_over_f"->N[cur,18],"front_z"->N[front,18],"min_action"->N[cc/.Last[opt],18],"min_cost_over_f"->N[First[opt],18],"residual_over_f"->N[First[opt]-Max[0,cur],18],"residual_scaled_lambda2"->N[(First[opt]-Max[0,cur])/lam^2,18],"reference_action_residual_over_f"->N[cost[a43]-Max[0,cur],18],"grid"->N[grid,18],"density"->N[fc43,18],"outer_order"->ord,"Fourier_step"->N[step],"seconds"->N[AbsoluteTime[]-start]|>];

matched43[nn_]:=Module[{z,target},z=zv/.FindRoot[Erfc[zv/Sqrt[2]]/2==Sqrt[mu43/nn]zv,{zv,3},WorkingPrecision->(wp43+30),AccuracyGoal->(wp43+10),PrecisionGoal->(wp43+10)];target=mu43-Sqrt[mu43/nn]z;a/.FindRoot[2a-(1+a/(1+Exp[a]))==target,{a,mu43-.01},WorkingPrecision->(wp43+30),AccuracyGoal->(wp43+10),PrecisionGoal->(wp43+10)]];
configure43[nn_,ss_]:=(N43=nn;a43=matched43[nn];pa43=1/(1+Exp[a43]);mean43=1+a43 pa43;ell43=2a43-mean43;source43=ss;kappa43=(mean43-a43)^2-(mu43-mean43));
If[!TrueQ[$finite43LibraryMode],args43=$ScriptCommandLine;out43=If[Length[args43]>=2,args43[[2]],"/tmp/finite_left43.json"];rows43={};Do[configure43[tt^power,ss];ll=N[(fraction/(aa mu43(mu43-1)/8 tt))^(1/3),60];row=case43[tt,ll,aa,32,1/8,True];AppendTo[rows43,row];Print[KeyTake[row,{"t","N","source","active","phase_fraction_alpha_J_over_D","floor_over_density","residual_scaled_lambda2","seconds"}]];Export[out43,<|"wolfram"->$Version,"status"->"finite matched left cost; floating sufficient-condition scan, not certified","floor_Lambda"->1,"working_precision"->wp43,"Fourier_cutoff"->cut43,"physical_jump_cutoff"->80,"action_grid"->"0 to 6 in increments of 1/4; stationary root, not a global numerical certificate","root_bisection_steps"->55,"rows"->rows43|>,"RawJSON"],{tt,{512,8192}},{power,{2,8}},{ss,{0,1/10}},{aa,{11/10}},{fraction,{1/4,1}}]];
