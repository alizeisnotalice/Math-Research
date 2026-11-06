(* r038: exact identities and stable Gaussian comparison scales. *)
Clear[x,m0,nn,ga];mf[x_]:=1+x/(1+Exp[x]);vf[x_]:=1+x^2 Exp[x]/(1+Exp[x])^2;
red[e_]:=FullSimplify[e/.{Exp[2m0]->1/(m0-1)^2,Exp[m0]->1/(m0-1)},m0>1];
checks38=<|"affine_regret"->Simplify[(1-r/nn+t m0/nn)-(1-(r-mc)/nn+(t-1)m0/nn)],
 "convex_gap_identity"->Simplify[lam (ko-kw)-(w-kw)-(lam(o-w)-(1-lam)(w-kw)-lam(o-ko))],
 "mean_at_mu"->red[mf[m0]],"mean_second_at_mu"->red[D[mf[x],{x,2}]/.x->m0],"variance_at_mu"->red[vf[m0]],"variance_first_at_mu"->red[D[vf[x],x]/.x->m0],
 "mean_loss_scaled_limit"->Limit[nn^(ga-1)Log[nn],nn->Infinity,Assumptions->0<ga<1],
 "centered_variance_scaled_limit"->Limit[nn^((ga-1)/2)Sqrt[Log[nn]],nn->Infinity,Assumptions->0<ga<1]|>;
calc38[n_Integer,prec_Integer]:=Module[{mu=N[1+ProductLog[1/E],prec],zz,dd,aa,tt,delta,gauss,old},
 zz=z/.FindRoot[Erfc[z/Sqrt[2]]/2==Sqrt[mu]z/Sqrt[n],{z,3},WorkingPrecision->prec];dd=Sqrt[mu n]zz;
 aa=a/.FindRoot[2a-mf[a]==mu-dd/n,{a,mu-.01},WorkingPrecision->prec];
 tt=Floor[n^(3/4)];delta=mu-mf[aa];gauss=tt delta+Sqrt[2/Pi]Sqrt[tt]Abs[Sqrt[mu]-Sqrt[vf[aa]]];old=tt Abs[aa-mu];
 {aa,delta,gauss,old,tt}];
rows38=Table[Module[{u=calc38[n,60],v=calc38[n,100]},<|"N"->n,"T"->u[[5]],"a"->N[u[[1]],16],"mean_loss"->N[u[[2]],16],"Gaussian_W1_upper"->N[u[[3]],16],"old_stepwise_scale"->N[u[[4]],16],"precision_doubling_difference"->ToString[N[Max[Abs[Take[u,4]-Take[v,4]]],8],InputForm]|>],{n,{2^16,2^24,2^32,2^40}}];
args38=$ScriptCommandLine;out38=If[Length[args38]>=2,args38[[2]],"/tmp/area_checks38.json"];
Export[out38,<|"wolfram"->$Version,"symbolic_checks"->Map[ToString[#,InputForm]&,checks38],"root_precision"->{60,100},"gamma"->.75,"rows"->rows38,"scope"->"Gaussian comparison term only; unscaled iid CLT adds a uniform O(1) term"|>,"RawJSON"];
tailRows=Table[Module[{poly,cal,qq,lead},
 qq=Sum[Binomial[tt,j]Sum[(r-c j)^i/i!,{i,0,tt-1}],{j,0,tt}];
 cal=Sum[Binomial[tt,j](tt Sum[(r-c j)^i/i!,{i,0,tt}]-(r-c j)Sum[(r-c j)^i/i!,{i,0,tt-1}]),{j,0,tt}];
 poly=Sum[Binomial[tt,j](r-c j)^tt/tt!,{j,0,tt}]+qq+cal/nn;
 <|"T"->tt,"call_polynomial_degree"->Exponent[Expand[cal],r],"payoff_degree"->Exponent[Expand[poly],r],"payoff_leading_coefficient_check"->TrueQ[FullSimplify[Coefficient[Expand[poly],r,tt]==2^tt/tt!]]|>],{tt,{1,2,4,8}}];
result=Import[out38,"RawJSON"];AssociateTo[result,"right_tail_checks"->tailRows];Export[out38,result,"RawJSON"];
Print[checks38];Print[rows38];Print[tailRows];
