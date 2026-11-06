(* Exact polynomial checks for fixed-time initialization; no Bellman certificate. *)
Clear[k,r,s,nn,a,ell];
small39=FullSimplify[Integrate[(r-s)s^(k-1),{s,0,r},Assumptions->r>0&&k>0]/r^(k-1),r>0&&k>0];
poly39[kk_Integer]:=Module[{ff,fc,vv,call,qq,jj,pa=1/(1+Exp[a])},
 ff=Sum[Binomial[kk,j](r-a j)^(kk-1),{j,0,kk}]/Gamma[kk];
 fc=Sum[Binomial[kk,j]((r-a j)^(kk+1)/(kk+1)+a^2 j(kk-j)(r-a j)^(kk-1)/(kk-1)),{j,0,kk}]/Gamma[kk];
 qq=Sum[Binomial[kk,j]Sum[(r-a j)^i/i!,{i,0,kk-1}],{j,0,kk}];
 call=Sum[Binomial[kk,j](kk Sum[(r-a j)^i/i!,{i,0,kk}]-(r-a j)Sum[(r-a j)^i/i!,{i,0,kk-1}]),{j,0,kk}];
 vv=Sum[Binomial[kk,j](r-a j)^kk/kk!,{j,0,kk}]+qq+call/nn;
 jj=pa/8(ff(r-kk ell)^2/kk+(HarmonicNumber[kk]-1)fc/kk);
 <|"k"->kk,"tail_slope_identity"->TrueQ[FullSimplify[Coefficient[Expand[jj],r,kk+1]/Coefficient[Expand[vv],r,kk]==pa/8(kk+HarmonicNumber[kk])/(kk+1)]],
 "density_payoff_ratio_identity"->TrueQ[FullSimplify[Coefficient[Expand[ff],r,kk-1]/Coefficient[Expand[vv],r,kk]==kk]],
 "small_r_deficit_over_density_coefficient"->ToString[(1-1/nn)/(kk(kk+1)),InputForm],
 "tail_slope_at_mu"->N[(pa/8(kk+HarmonicNumber[kk])/(kk+1))/.a->1+ProductLog[1/E],17]|>];
rows39=poly39/@{4,8,16};
gaussianCorrection=FullSimplify[Integrate[Exp[-nu w^2/2](ka4 w^4/24-ka3^2 w^6/72+ka3 zz w^4/6-zz^2 w^2/2),{w,-Infinity,Infinity},Assumptions->nu>0]/Integrate[Exp[-nu w^2/2],{w,-Infinity,Infinity},Assumptions->nu>0],nu>0];
h0check=FullSimplify[Integrate[Exp[et z],{z,-Infinity,0},Assumptions->0<et<1]+Integrate[Exp[(et-1)z],{z,0,Infinity},Assumptions->0<et<1],0<et<1];
auxcheck=FullSimplify[D[Log[(q+p Exp[et a])/(1-et)],{et,2}]==1/(1-et)^2+a^2 p q Exp[et a]/(q+p Exp[et a])^2,p>0&&q>0&&a>0&&0<et<1];
args39=$ScriptCommandLine;out39=If[Length[args39]>=2,args39[[2]],"/tmp/initialization_checks39.json"];
Export[out39,<|"wolfram"->$Version,"precision"->"exact identities; 17-digit display of exact constants","beta_integral_identity"->ToString[small39,InputForm],"rows"->rows39,"Fourier_Gaussian_correction"->ToString[gaussianCorrection,InputForm],"terminal_h0"->ToString[h0check,InputForm],"auxiliary_variance_identity"->auxcheck|>,"RawJSON"];
Print[small39];Print[rows39];
