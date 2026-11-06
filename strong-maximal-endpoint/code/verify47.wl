(* Exact identities, independent finite mixture, and numerical refinements. *)
$plateau47LibraryMode=True;Get[FileNameJoin[{DirectoryName[$InputFileName],"plateau47.wl"}]];$plateau47LibraryMode=False;
check47=FullSimplify[D[h(1+Log[h])^(nn-1),{h,2}]/.h->1,Assumptions->nn>=2];
ids47=<|"Phi_second_derivative_at_1_plus"->TrueQ[check47==(nn-1)^2],"optimal_parameter_quadratic"->TrueQ[Expand[-pp tt^2/2-pp ss tt/2-pp ss^2/8+pp(tt+ss/2)^2/2]==0],"plateau_scaling"->TrueQ[FullSimplify[aa(1+Log[aa])^(nn-1)==Exp[ee](1+ee)^(nn-1)/.aa->Exp[ee],Assumptions->ee>=0]]|>;
a47=run47[256];q47=ToExpression[a47["q_internal"]];mu47=N[1+ProductLog[1/E],180];n47=256;w47=N[Log[n47]/2+Log[Log[n47]],150];
direct47[c_,q_,w_,wp_]:=Module[{rho=1-1/n47,ks=Range[0,n47],weights,aa,ff},weights=N[Binomial[n47,#]&/@Range[0,n47],wp]Exp[-c ks]/(1+Exp[-c])^n47;
 aa=Exp[-(q+w)/n47]rho^(-n47)Dot[weights Exp[c ks/n47],GammaRegularized[n47,rho Clip[q+w-c ks,{0,Infinity}]]];
 ff[t_]:=Dot[weights,GammaRegularized[n47,0,Clip[t-c ks,{0,Infinity}]]];aa+ff[q+w]-ff[q]-1];
mixture47=Table[Module[{c=mu47+tt Sqrt[Log[n47]/n47],ff,di,diff},ff=mean47[n47,c,q47,w47,120]["mean_minus_one"];di=direct47[c,q47,w47,120];diff=Abs[ff-di];<|"t_over_sqrtmu"->tt/Sqrt[mu47],"absolute_mean_difference"->ToString[diff,InputForm],"less_than_1e-30"->TrueQ[diff<10^-30]|>],{tt,{-Sqrt[mu47],0,Sqrt[mu47]/2}}];
refine47=Table[Module[{a=run47[n],b=run47[n,2,1/32,40],dif},dif=Max[Abs[a["rows"][[All,"fixed_input_mean_deficit_normalized"]]-b["rows"][[All,"fixed_input_mean_deficit_normalized"]]]];<|"N_exact"->ToString[n,InputForm],"maximum_normalized_deficit_difference"->ToString[dif,InputForm],"less_than_1e-20"->TrueQ[dif<10^-20],"both_root_residuals_valid"->TrueQ[a["root_valid"]&&b["root_valid"]]|>],{n,{256,2^48,2^256,2^1024}}];
res47=<|"wolfram"->$Version,"symbolic_checks"->ids47,"finite_original_mixture_check"->mixture47,"precision_mesh_cutoff_refinements"->refine47,"status"->"completed; floating refinements and exact identities, not interval certificate"|>;
outV47=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/verify47.json"];Export[outV47,res47,"RawJSON"];Print[res47];
