mu=1+ProductLog[Exp[-1]];
root=DirectoryName[DirectoryName[$InputFileName]];old=Import[FileNameJoin[{root,"rounds","r019","checks.json"}],"RawJSON"];
rows=Table[Module[{initial,sols,softs,exactK,offsets},
 initial=SelectFirst[old["rows"],#["n"]==nn&]["plateau_root"];
 {sols,softs}=Transpose@Table[Module[{cc=N[mu,wp],pp,weights,w=N[Log[nn]/2+Log[Log[nn]],wp],aa,ff,bb},
 pp=1/(1+Exp[cc]);weights=Table[Binomial[nn,k]pp^k(1-pp)^(nn-k),{k,0,nn}];
 aa[q_?NumericQ]:=Exp[-q/nn]((nn-1)/nn)^(-nn) Sum[weights[[k+1]]Exp[cc k/nn]GammaRegularized[nn,(nn-1)/nn Max[q-cc k,0]],{k,0,nn}];
 ff[q_?NumericQ]:=Sum[weights[[k+1]]If[q<=cc k,0,GammaRegularized[nn,0,q-cc k]],{k,0,nn}];
 bb[q_?NumericQ]:=aa[q+w]+ff[q+w]-ff[q]-1;
 {q/.FindRoot[bb[q]==0,{q,SetPrecision[initial,wp]},WorkingPrecision->wp,AccuracyGoal->Floor[wp/2],PrecisionGoal->Floor[wp/2]],
 q/.FindRoot[aa[q]-1==0,{q,SetPrecision[initial+1,wp]},WorkingPrecision->wp,AccuracyGoal->Floor[wp/2],PrecisionGoal->Floor[wp/2]]}],{wp,{60,100}}];
 exactK=N[1+Exp[-(Log[nn]/2+Log[Log[nn]])](nn Exp[nn-1]Gamma[nn,nn-1]/(nn-1)^nn-1),100];
 offsets=sols[[2]]-softs[[2]]-Log[exactK];
 <|"n"->nn,"independent_plateau_root_80"->ToString[N[sols[[2]],80],InputForm],"root_precision_doubling_stable_25_digits"->TrueQ[Abs[sols[[1]]-sols[[2]]]<10^-25],"root_minus_soft"->ToString[N[sols[[2]]-softs[[2]],30],InputForm],"score_minus_soft"->ToString[N[offsets,30],InputForm],"log_modular_factor"->ToString[N[Log[exactK],30],InputForm]|>],{nn,{256,1024}}];
result=<|"wolfram_version"->$Version,"working_precisions"->{60,100},"method"->"full original-binomial gamma sums; exact modular factor","rows"->rows|>;
Export[FileNameJoin[{root,"rounds","r019","symbolic_checks.json"}],result,"RawJSON"];
Print[ExportString[result,"RawJSON"]];
