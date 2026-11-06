(* Exact identities for centered reference/payoff bridge and the local gain. *)
Clear[x,nu,k3,A1,A2,A3,eta,m,u,v,nn,z];
expect41[p_]:=Expand[p]/.{x^4->3nu^2,x^3->-k3,x^2->nu,x->0};
B41=1/2-k3 x/(2nu^2)-x^2/(2nu);
H41=1+A1 x+A2(x^2-nu)/2+A3(x^3+k3)/6;
hN41=FullSimplify[Integrate[(1-z/nn)Exp[eta z],{z,-Infinity,0},Assumptions->nn>0&&0<eta<1]+Integrate[Exp[(eta-1)z],{z,0,Infinity},Assumptions->0<eta<1],nn>0&&0<eta<1];
weights41={1/(1-eta),1/eta,1/(nn eta^2)}/hN41;
lam41[w_]:=-Log[1-eta]+Log[1+(m-1)Exp[(m-w)eta+w eta^2]]-Log[1+(m-1)Exp[-w eta]];
y41[w_]:=1/(1-eta)+(m+w eta)(m-1)Exp[(m-w)eta+w eta^2]/(1+(m-1)Exp[(m-w)eta+w eta^2]);
pa41=(m-1)Exp[-u eta]/(1+(m-1)Exp[-u eta]);ma41=1+(m+u eta)pa41;ell41=2(m+u eta)-ma41;
D41=pa41/8 eta(1-eta)(y41[u]-ell41)^2;
gain41=FullSimplify[Coefficient[Normal[Series[lam41[v]-lam41[u]-D41,{eta,0,3}]],eta,3],m>1];
meanDiff41=FullSimplify[Normal[Series[y41[v]-y41[u],{eta,0,2}]],m>1];
checks41=<|"density_bridge_B_mean"->FullSimplify[expect41[B41]],"payoff_bridge_polynomial_mean"->FullSimplify[expect41[H41]],"finite_terminal_normalizer"->ToString[hN41,InputForm],"finite_terminal_mixture_weights_sum"->FullSimplify[Total[weights41],nn>0&&0<eta<1],"gain_minus_source_eta3"->ToString[gain41,InputForm],"negative_square_gain_identity"->TrueQ[FullSimplify[gain41==-(m-1)/(2m)(v-m/2)^2,m>1]],"tilted_mean_difference_through_eta2"->ToString[meanDiff41,InputForm]|>;

muNum41=N[1+ProductLog[1/E],100];
lamNum41[c_,e_]:=-Log[1-e]+Log[1+Exp[(e-1)c]]-Log[1+Exp[-c]];
meanNum41[c_,e_]:=1/(1-e)+c/(1+Exp[(1-e)c]);
varNum41[c_,e_]:=1/(1-e)^2+c^2 Exp[(1-e)c]/(1+Exp[(1-e)c])^2;
thirdNum41[c_,e_]:=2/(1-e)^3+c^3 Exp[(1-e)c](Exp[(1-e)c]-1)/(1+Exp[(1-e)c])^3;
opt41[e_]:=Module[{ee=N[e,100]},ct/.FindRoot[(Exp[(1-ee)ct]-(1-ee)Exp[ct]+ee)/ee==0,{ct,muNum41+muNum41 ee/2},WorkingPrecision->80,AccuracyGoal->55,PrecisionGoal->55]];
gainRows41=Flatten[Table[Module[{ee=N[ep,100],aa,co,yy,ell,dd,gg},aa=muNum41+hh ee;co=opt41[ee];yy=meanNum41[aa,ee];ell=2aa-meanNum41[aa,0];dd=ee(1-ee)(yy-ell)^2/(8(1+Exp[aa]));gg=lamNum41[co,ee]-lamNum41[aa,ee];<|"eta"->N[ee,18],"reference_offset_over_eta"->hh,"gain_minus_source_over_eta4"->N[(gg-dd)/ee^4,24],"gain_minus_1_1_source_over_eta3"->N[(gg-(11/10)dd)/ee^3,24]|>],{ep,{1/200,1/100,1/50,1/25}},{hh,{-2,-1,-1/2,0}}],1];
polyInt41[lo_,hi_,ll_]:=(Gamma[1,ll lo]-If[hi===Infinity,0,Gamma[1,ll hi]])/ll+(Gamma[5,ll lo]-If[hi===Infinity,0,Gamma[5,ll hi]])/ll^5;
weightedTV41[aa_,cc_,ee_]:=Module[{qa=1/(1+Exp[-aa]),qc=1/(1+Exp[-cc]),lo=Min[aa,cc],hi=Max[aa,cc],mid,ll=1-ee},mid=If[cc<aa,2qc-qa,qc-2qa];(Abs[qc-qa]polyInt41[0,lo,ll]+Abs[mid]polyInt41[lo,hi,ll]+2Abs[qc-qa]polyInt41[hi,Infinity,ll])/Exp[lamNum41[aa,ee]]];
momentRows41=Flatten[Table[Module[{ee=N[ep,100],aa,cc,yy,vv,kk,mc,vc,kc,rat,dx},aa=muNum41-ee;cc=aa+jj ee;yy=meanNum41[aa,ee];vv=varNum41[aa,ee];kk=thirdNum41[aa,ee];mc=meanNum41[cc,ee];vc=varNum41[cc,ee];kc=thirdNum41[cc,ee];rat=Exp[lamNum41[cc,ee]-lamNum41[aa,ee]];dx=yy-mc;<|"eta"->N[ee,18],"action_shift_over_eta"->jj,"mass_change_over_eta3"->N[(rat-1)/ee^3,24],"centered_first_moment_over_eta2"->N[rat dx/ee^2,24],"centered_second_moment_over_eta"->N[rat(vc+dx^2-vv)/ee,24],"centered_third_moment_over_eta"->N[rat(dx^3+3dx vc-kc+kk)/ee,24],"weighted_TV_over_abs_action_shift"->N[weightedTV41[aa,cc,ee]/Abs[cc-aa],24]|>],{ep,{1/200,1/100,1/50,1/25}},{jj,{-2,-1,1,2}}],1];
args41=$ScriptCommandLine;out41=If[Length[args41]>=2,args41[[2]],"/tmp/bridge_checks41.json"];
Export[out41,<|"wolfram"->$Version,"arithmetic"->"exact symbolic under m>1, exp(-m)=m-1","checks"->checks41,"numerical_working_precision"->80,"gain_grid"->gainRows41,"moving_jump_moment_grid"->momentRows41|>,"RawJSON"];Print[checks41];
