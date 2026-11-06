(* Exact identities and asymptotic scale checks for r037; not a Bellman certificate. *)
Clear[l,v,u,x,n,g,k,c,t];
bb=l/(1+(l-v)/v Exp[-u l/(v(l-v))]);
ass=0<v<l && Element[{l,v,u},Reals];
checks=<|
 "initial_value"->FullSimplify[(bb/.u->0)==v,ass],
 "initial_derivative"->FullSimplify[(D[bb,u]/.u->0)==1,ass],
 "second_derivative"->FullSimplify[(D[bb,{u,2}]/.u->0)==(l-2v)/(v(l-v)),ass],
 "complement_symmetry"->FullSimplify[l-bb==(bb/.{v->l-v,u->-u}),ass],
 "positive_derivative"->FullSimplify[D[bb,u]>0,ass],
 "bernoulli_call_area"->Integrate[(1-Abs[x])/2,{x,-1,1}],
 "tail_integral"->FullSimplify[Integrate[x Exp[-c x^2/t],{x,0,t}]+Integrate[x Exp[-c x],{x,t,Infinity}],c>0&&t>0],
 "late_tail_ratio_limit"->Limit[n^(2-g(k+1))/Log[n],n->Infinity,Assumptions->0<g<1&&k>2/g-1],
 "matched_reference_error_limit"->Limit[n^g Sqrt[Log[n]/n],n->Infinity,Assumptions->0<g<1/2],
 "critical_logistic_error_over_density"->Limit[Log[n]^(5/2)/Sqrt[n],n->Infinity]|>;
out=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"extension_checks37.json"];
Export[out,<|"wolfram"->$Version,"precision"->"exact symbolic","checks"->Map[ToString[#,InputForm]&,checks]|>,"RawJSON"];
Print[checks];
mu=1+ProductLog[1/Exp[1]];mean[cc_]:=1+cc/(1+Exp[cc]);
rows=Table[Module[{zz,dd,target,cc},
 zz=z/.FindRoot[Erfc[z/Sqrt[2]]/2==Sqrt[mu]z/Sqrt[nn],{z,2},WorkingPrecision->50];
 dd=Sqrt[mu nn]zz;target=mu-dd/nn;
 cc=c0/.FindRoot[mean[c0]==target,{c0,N[mu+.5,50]},WorkingPrecision->50];
 <|"N"->nn,"z"->N[zz,16],"action"->N[cc,16],"mean_residual"->ToString[N[mean[cc]-target,8],InputForm],
 "action_offset_over_fourth_scale"->N[(cc-mu)/(Log[nn]/nn)^(1/4),16],
 "central_over_reference_density_scale"->N[Sqrt[nn]/Log[nn],16]|>],{nn,{4096,65536,1048576,16777216}}];
Export[out,<|"wolfram"->$Version,"precision"->"exact symbolic; root checks 50 digits",
 "checks"->Map[ToString[#,InputForm]&,checks],"mean_shift_obstruction_checks"->rows|>,"RawJSON"];
Print[rows];
