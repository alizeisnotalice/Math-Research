(* Exact formula guards plus numerical derivative-sign diagnostics. *)
ClearAll[guard,a,u,beta,b];guard[name_,test_]:=If[TrueQ[test],Print["PASS ",name],Print["FAIL ",name];Exit[1]];
M[a_,u_,beta_]:=a^beta+(u-a)^beta;Den[a_,beta_]:=a^beta+(1-a)^beta;
KK[a_,u_,beta_]:=D[M[x,u,beta],x]Den[a,beta]-M[a,u,beta]D[Den[x,beta],x]/.x->a;
KprimeExpected=beta(beta-1)((a^(beta-2)+(u-a)^(beta-2))Den[a,beta]-M[a,u,beta](a^(beta-2)+(1-a)^(beta-2)));
guard["Kprime_formula",TrueQ[FullSimplify[D[KK[a,u,beta],a]==KprimeExpected]]];
ah=u/(2(1+Sqrt[2u]+u));oh=(u+Sqrt[2u])/(1+Sqrt[2u]);guard["half_K_zero",TrueQ[FullSimplify[KK[ah,u,1/2]==0,0<u<1]]];guard["half_Omega_closed",TrueQ[FullSimplify[M[ah,u,1/2]/Den[ah,1/2]==oh,0<u<1]]];guard["half_u_1_8",TrueQ[FullSimplify[(ah/.u->1/8)==1/26&&(oh/.u->1/8)==5/12]]];
grid=Flatten[Table[{bb,uu},{bb,{1/10,1/4,1/2,3/4,9/10}},{uu,{2^-4,2^-16,2^-64,1-2^-4,1-2^-16,1-2^-64}}],1];checks=Table[Module[{bb=N[pair[[1]],100],uu=N[pair[[2]],100],vals},vals=Table[N[(D[KK[x,uu,bb],x]/.x->(uu j/400)),80],{j,1,199}];Max[vals]<0],{pair,grid}];guard["Kprime_sampled_negative",And@@checks];Print["PASS round32 symbolic audit; derivative-grid cases=",Length[grid]];
(* Independent algebra for the Cartesian-grid lower construction. *)
ClearAll[ell,eta,nn,gg,bb,kk,TT,hh];
TT=kk(1+ell)+1;hh=(1-eta)(2/ell)^(nn-1);
Jgrid=(kk ell/TT)^(nn-1) hh^bb;
expected=(1-eta)^bb*(2^bb ell^(1-bb)/(1+ell))^(nn-1);
(* Positive rational specialisations avoid branch ambiguity in fractional powers. *)
gridLimitChecks=Flatten[Table[TrueQ[FullSimplify[Limit[Jgrid/.{nn->n0,bb->b0,ell->l0,eta->e0},kk->Infinity]==expected/.{nn->n0,bb->b0,ell->l0,eta->e0}]],{n0,{2,3,4}},{b0,{1/5,1/2,4/5}},{l0,{1/3,2,7/3}},{e0,{1/10,1/4}}],Infinity];
guard["grid_J_limit_positive_specialisations",And@@gridLimitChecks];
shape=ell^gg/(1+ell);guard["grid_optimal_ell",TrueQ[FullSimplify[(D[shape,ell]/.ell->gg/(1-gg))==0,0<gg<1]]];
lambda=2^(1-gg)(1-gg)^(1-gg)gg^gg;
guard["Lambda_gamma_1_5",TrueQ[FullSimplify[(lambda/.gg->1/5)^5==4096/3125]]&&TrueQ[4096/3125>1]];
Jfinite=(2(3/38)/(2(1+3/38)+1))*25^(19/20);
guard["grid_finite_J_power",TrueQ[FullSimplify[Jfinite^20==5^18/2^40]]&&TrueQ[5^18/2^40>2]];
Print["PASS round32 grid-lower algebra audit"];
