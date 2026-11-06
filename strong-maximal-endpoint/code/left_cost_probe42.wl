(* Actual clipped nonlinear cost; negative tilted Fourier quadrature. Floating evidence. *)
mu42=N[1+ProductLog[1/E],70];pa42=1/(1+Exp[mu42]);wp42=50;
lam42[e_]:=-Log[1-e]+Log[(1+Exp[(e-1)mu42])/(1+Exp[-mu42])];
y42[e_]:=1/(1-e)+mu42 Exp[(e-1)mu42]/(1+Exp[(e-1)mu42]);
h42[e_]:=-1/(e(1-e));g142[e_]:=1/(1-e)^2;g242[e_]:=Exp[mu42 e]/(1-pa42+pa42 Exp[mu42 e])^2;
exprel42[x_]:=If[Abs[x]<10^-12,1+x/2+x^2/6+x^3/24+x^4/120,2Exp[x/2]Sinh[x/2]/x];
std42[n_]:=std42[n]=Module[{b,j,ev,vec,ix},b=N[Range[n-1]/Sqrt[4Range[n-1]^2-1],wp42];j=DiagonalMatrix[b,1]+DiagonalMatrix[b,-1];{ev,vec}=Eigensystem[j];ix=Ordering[ev];{ev[[ix]],2vec[[ix,1]]^2}];
quad42[n_,aa_,bb_]:=Module[{v=std42[n]}, {(bb-aa)v[[1]]/2+(aa+bb)/2,(bb-aa)v[[2]]/2}];
Clear[setup42];
setup42[tt_,lambda_,step_:1/8]:=Module[{e=-N[lambda,wp42+15],vv,ss,ps,cur,pre,n=tt,fac},
 e42=e;t42=tt;yy42=y42[e];rr42=tt yy42;vv=N[Range[-18/step,18/step]step,wp42];ss=e+I vv/Sqrt[tt];
 ps=Exp[tt(Log[(1-e)/(1-ss)]+Log[(1+Exp[(ss-1)mu42])/(1+Exp[(e-1)mu42])]-I yy42 vv/Sqrt[tt])];
 cur=Re[Total[ps]] step/(2Pi Sqrt[tt]);fc42=Exp[tt lam42[e]-e rr42]cur;
 dn42=Re[Total[ps(h42/@ss)]] step/(2Pi Sqrt[tt])/cur;
 cn42=Re[Total[ps((g142/@ss)+mu42^2 pa42(1-pa42)(g242/@ss))]]step/(2Pi Sqrt[tt])/cur;
 jn42=pa42/8(tt(yy42-mu42)^2+(HarmonicNumber[tt]-1)cn42);
 (* same eta, previous time; normalized by current physical density *)
 vv42=vv;ss42=e+I vv/Sqrt[tt-1];pre=Exp[(tt-1)(Log[(1-e)/(1-ss42)]+Log[(1+Exp[(ss42-1)mu42])/(1+Exp[(e-1)mu42])]-I yy42 vv/Sqrt[tt-1])];
 pre42=pre step/(2Pi Sqrt[tt-1])/cur;
 dpre42=pre42(h42/@ss42);cpre42=pre42((g142/@ss42)+mu42^2 pa42(1-pa42)(g242/@ss42));
 {dn42,jn42}];
raw42[z_?NumericQ,alpha_]:=Module[{ph,fr,dr,cr,jr,x},ph=Exp[-I vv42(yy42-z)/Sqrt[t42-1]];fr=Re[Total[pre42 ph]];dr=Re[Total[dpre42 ph]];cr=Re[Total[cpre42 ph]];jr=pa42/8(fr(rr42-z-(t42-1)mu42)^2/(t42-1)+(HarmonicNumber[t42-1]-1)cr);{fr,dr,jr}=Exp[-lam42[e42]+e42 z]{fr,dr,jr};x=alpha fc42 jr/(1-fc42 dr);dr-alpha jr exprel42[x]];
rawCurrent42[alpha_]:=Module[{x=alpha fc42 jn42/(1-fc42 dn42)},dn42-alpha jn42 exprel42[x]];
case42[tt_,lam_,alpha_,ord_:32,step_:1/8,clipped_:True]:=Module[{d,j,front,lo,hi,mid,q,qn,zn,wn,b0,bb,cost,grid,best,opt,cur,left,right,der,ca,start=AbsoluteTime[]},
 {d,j}=setup42[tt,lam,step];cur=rawCurrent42[alpha];front=80;
 If[clipped,If[raw42[0,alpha]<=0,front=0,If[raw42[80,alpha]<0,lo=0;hi=80;Do[mid=(lo+hi)/2;If[raw42[mid,alpha]>0,lo=mid,hi=mid],{55}];front=(lo+hi)/2]]];
 bb[cc_?NumericQ]:=If[cc>=front,0,Module[{x,w},{x,w}=quad42[ord,cc,front];Dot[w,Exp[-x](If[clipped,Max[0,raw42[#,alpha]],raw42[#,alpha]]&/@x)]]];
 b0=bb[0];cost[cc_?NumericQ]:=(b0+bb[cc])/(1+Exp[-cc]);
 grid=Table[{N[cc],cost[N[cc,wp42]]},{cc,0,6,1/4}];best=First[MinimalBy[grid,Last]];
 der[cc_?NumericQ]:=cost[cc]-If[clipped,Max[0,raw42[cc,alpha]],raw42[cc,alpha]];left=N[0,wp42];right=N[6,wp42];If[der[left]<=0&&der[right]>=0,Do[mid=(left+right)/2;If[der[mid]<0,left=mid,right=mid],{55}];ca=(left+right)/2,ca=N[best[[1]],wp42]];opt={cost[ca],{cc->ca}};
 <|"working_precision"->wp42,"clipped"->clipped,"t"->tt,"lambda"->N[lam,18],"alpha"->N[alpha],"phase_fraction_alpha_J_over_D"->N[alpha j/d,18],"current_cost_over_f"->N[cur,18],"front_z"->N[front,18],"min_action"->N[cc/.Last[opt],18],"min_cost_over_f"->N[First[opt],18],"residual_over_f"->N[First[opt]-Max[0,cur],18],"residual_scaled_lambda2"->N[(First[opt]-Max[0,cur])/lam^2,18],"reference_action_residual_over_f"->N[cost[mu42]-Max[0,cur],18],"grid"->N[grid,18],"density"->N[fc42,18],"outer_order"->ord,"Fourier_step"->N[step],"seconds"->N[AbsoluteTime[]-start]|>];
If[!TrueQ[$left42LibraryMode],
 args42=$ScriptCommandLine;out42=If[Length[args42]>=2,args42[[2]],"/tmp/left_cost_probe42.json"];
 rows42={};Do[lam=N[(frac/(alpha mu42(mu42-1)/8 tt))^(1/3),60];row=case42[tt,lam,alpha];AppendTo[rows42,row];Print[KeyTake[row,{"t","lambda","alpha","phase_fraction_alpha_J_over_D","residual_scaled_lambda2","min_action","seconds"}]];Export[out42,<|"wolfram"->$Version,"status"->"negative Fourier clipped cost scan; floating not certified","precision"->wp42,"Fourier_cutoff"->18,"physical_jump_cutoff"->80,"root_bisection_steps"->55,"action_grid"->"0 to 6 in increments of 1/4; stationary root only, not global certificate","rows"->rows42|>,"RawJSON"],{tt,{512,2048,8192}},{alpha,{101/100,11/10}},{frac,{1/4,3/4,1}}];
];
