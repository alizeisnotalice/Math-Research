mu=1+ProductLog[1/E];q=(mu-1)/mu;t=1-2q;
kap[j_]:=(1-t^(2j))/2-2q(1-q)^2 j t^(2j-2);
rows=Flatten[Table[Module[{qq=N[q,prec+20],uu,solutions,values,lg},
 solutions=Table[pp/.FindRoot[Log[(1-pp)/pp]-1/(1-pp)==2 s eps (1-2pp)^(j-1),{pp,qq},WorkingPrecision->prec,AccuracyGoal->prec-20],{s,{1,-1}}];
 values=MapThread[j(1+#1 Log[(1-#1)/#1]-N[mu,prec+20])+#2 eps (1-2#1)^j&,{solutions,{1,-1}}];
 lg=Log[(Exp[values[[1]]]+Exp[values[[2]]])/2]-Log[Cosh[eps]];
 <|"order"->j,"epsilon"->ToString[eps,InputForm],"precision"->prec,"log_ratio"->ToString[lg,InputForm],
 "negative_quadratic_diagnostic"->ToString[-lg/eps^2,InputForm],"predicted_kappa"->ToString[N[kap[j],prec],InputForm],
 "global_concavity_margin"->ToString[27/4-4eps(j-1),InputForm]|>],{j,1,8},{eps,{1/100,1/1000}},{prec,{60,100}}]];
(* Symbolic check of the strictly positive spectral increment formula. *)
Clear[z,kk];tt=1-2z;aa=tt^2;bb=2z(1-z)^2;
positive=2z(1-z)(1-z-tt^3)+bb(kk-1)(1-aa);
identity=FullSimplify[aa(1-aa)/2+bb(kk(1-aa)-aa)-positive];
out=<|"version"->$Version,"spectral_increment_identity"->ToString[identity,InputForm],
 "kappa"->Table[<|"order"->j,"value"->ToString[N[kap[j],80],InputForm]|>,{j,1,8}],"rows"->rows|>;
Export[FileNameJoin[{DirectoryName[$InputFileName],"walsh_wolfram.json"}],out,"RawJSON"];Print[<|"identity"->identity,"rows"->Length[rows],"kappa1"->N[kap[1],30],"kappa4"->N[kap[4],30]|>];Quit[];
