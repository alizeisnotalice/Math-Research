(* Exact common-translation bookkeeping and a genuine clipping counterexample. *)
Clear[a,a1,a2,t1,t2,u1,u2,tau,AA,ell,uu,tt];
amp49[t1_,t2_,u1_,u2_,tau_]:=P[t1-u1]Q[t2-u2]Exp[u1+u2-tau];
checksG49=<|"weighted_Jacobian"->TrueQ[FullSimplify[Exp[-(t1+a1)-(t2+a2)]==Exp[-a1-a2]Exp[-t1-t2]]],"sourceAmp_with_tau_shift"->TrueQ[FullSimplify[amp49[t1,t2,u1+a1,u2+a2,tau+a1+a2]-amp49[t1-a1,t2-a2,u1,u2,tau]]==0]|>;
bneg49=FullSimplify[Integrate[AA Exp[-tt],{tt,-uu,ell-uu}],Assumptions->{AA>1,ell>0,uu<=0}];
bpos49=FullSimplify[Integrate[AA Exp[-tt],{tt,0,ell-uu}],Assumptions->{AA>1,0<=uu<=ell}];
lo49=-Log[AA(1-Exp[-ell])];hi49=ell+Log[1-1/AA];
vol49=FullSimplify[Integrate[Exp[-uu],{uu,lo49,hi49}],Assumptions->{AA>1,ell>Log[AA/(AA-1)]}];
checksG49=Join[checksG49,<|"negative_u_receiver_average"->TrueQ[FullSimplify[bneg49-AA Exp[uu](1-Exp[-ell])]==0],"positive_u_receiver_average"->TrueQ[FullSimplify[bpos49-AA(1-Exp[uu-ell])]==0],"halo_volume_formula"->TrueQ[FullSimplify[vol49-(AA-AA^2 Exp[-ell]/(AA-1)),Assumptions->{AA>1,ell>0}]==0]|>];
ratio49[x_]:=FullSimplify[(2-4Exp[-x])/(2(1-Exp[-x]))];
outG49=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/gauge49.json"];
resG49=<|"wolfram"->$Version,"exact_checks"->checksG49,"gauge"->"H_a(t)=H(t-a), receiver labels u_j+a, tau'=tau+a1+a2; full-domain J and V both multiply exp(-a1-a2)","clipping_counterexample"-><|"model"->"m=1, T~Exp(1), Phi_1(H)=H, H=2 1_[0,l]","original_l"->2,"translation_a"->-1,"clipped_l"->1,"original_actual_V_over_J"->N[ratio49[2],40],"clipped_actual_V_over_J"->N[ratio49[1],40],"difference_exact"->ToString[FullSimplify[ratio49[2]-ratio49[1]],InputForm],"difference_positive"->TrueQ[FullSimplify[ratio49[2]-ratio49[1]>0]]|>,"scope"->"full weighted-domain translation gauge; clipping is a new shape, not an invariant coordinate change; no C_N improvement"|>;
Export[outG49,resG49,"RawJSON"];Print[checksG49];
