(* Exact inputs; independent 60/100 digit full finite mixture roots. *)
mu=1+ProductLog[1/E];p[c_]:=1/(1+Exp[c]);m[c_]:=1+c p[c];v[c_]:=1+c^2 p[c](1-p[c]);
local=FullSimplify[{D[m[c],{c,2}],D[v[c],c]}/.c->mu];
Clear[calc];
calc[n_Integer,prec_Integer]:=Module[{u,z,h,eta,ce,pp,ww,rate,avg,q,qq,guess,mn,dn},
 u=N[mu,prec+30];
 z=zz/.FindRoot[Erfc[zz/Sqrt[2]]/2==Sqrt[u/n] zz,{zz,2},WorkingPrecision->prec+20,AccuracyGoal->prec];
 h=1/(1+Sqrt[u/n]/(Exp[-z^2/2]/Sqrt[2 Pi]));eta=z h/(2 Sqrt[u n]);
 ce=cc/.FindRoot[Evaluate[D[m[cc]-eta v[cc],cc]]==0,{cc,u},WorkingPrecision->prec+10,AccuracyGoal->prec-10];
 pp=1/(1+Exp[ce]);ww=Table[Binomial[n,k]pp^k(1-pp)^(n-k),{k,0,n}];rate=1-1/n;
 avg[q_?NumericQ]:=Sum[ww[[k+1]] Exp[-(q-ce k)/n-n Log[rate]] Gamma[n,rate Max[q-ce k,0]]/Gamma[n],{k,0,n}];
 mn=n(1+ce pp);dn=Sqrt[u n]z;guess=mn-dn+(3-u)z^2/6;
 qq=q/.FindRoot[avg[q]==1,{q,guess},WorkingPrecision->prec,AccuracyGoal->prec-20];
 <|"n"->n,"precision"->prec,"c_eta"->ToString[ce,InputForm],"z"->ToString[z,InputForm],
 "soft_centered"->ToString[qq-n u,InputForm],"root_residual"->ToString[Abs[avg[qq]-1],InputForm]|>
];
rows=Flatten[Table[calc[n,prec],{n,{64,256}},{prec,{60,100}}]];
large=Table[Module[{n=10^j,z,h},z=zz/.FindRoot[Log[Erfc[zz/Sqrt[2]]/2]==Log[Sqrt[mu/n]zz],{zz,Sqrt[Log[n]]},WorkingPrecision->100];h=1/(1+Sqrt[mu/n]/(Exp[-z^2/2]/Sqrt[2 Pi]));<|"log10_n"->j,"z"->ToString[z,InputForm],"H_constant_diagnostic"->ToString[z^2(h^2-1),InputForm]|>],{j,{10,50,200,1000}}];
out=<|"mu"->ToString[N[mu,80],InputForm],"B_Q"->ToString[N[(13mu-15)/12,80],InputForm],
 "B_soft"->ToString[N[(5mu-6)/6-(mu-1)/4,80],InputForm],
 "m_second_v_first"->ToString[local,InputForm],
 "m_second_identity"->ToString[FullSimplify[local[[1]]+(mu-1)/mu],InputForm],
 "v_first_identity"->ToString[FullSimplify[local[[2]]-(mu-1)],InputForm],
 "H_expansion"->ToString[Series[(1+x- x^2+3x^3)^(-2),{x,0,2}],InputForm],
 "constant_sum"->ToString[FullSimplify[(5u-6)/6-(u-1)/4+(u-1)/2-(13u-15)/12],InputForm],"large_dimension_Mills_checks"->large,"roots"->rows|>;
Export[FileNameJoin[{DirectoryName[$InputFileName],"wolfram.json"}],out,"RawJSON"];Print[out];Quit[];
