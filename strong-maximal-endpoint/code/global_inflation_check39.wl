(* r039 finite-rate checks for the proposed intermediate-time obstruction.
   WolframScript calculation only; no interval claim. *)
$MaxExtraPrecision = 100000;
args39=$ScriptCommandLine;out39=If[Length[args39]>=2,args39[[2]],FileNameJoin[{$TemporaryDirectory,"global_inflation_check39.json"}]];wp39=If[Length[args39]>=3,ToExpression[args39[[3]]],80];
mu = x /. FindRoot[x == 1 + Exp[-x], {x, 1.3}, WorkingPrecision -> (wp39+20)];
p = 1/(1+Exp[mu]);
lm[c_,e_] := Log[(1+Exp[(e-1)c])/((1-e)(1+Exp[-c]))];
yp[e_] := 1/(1-e)+mu Exp[(e-1)mu]/(1+Exp[(e-1)mu]);
cp[e_?NumericQ] := Module[{ee=SetPrecision[e,wp39+20]},c /. FindRoot[Exp[(1-ee)c]-(1-ee)Exp[c]+ee==0,{c,mu+mu ee/2},WorkingPrecision->wp39,AccuracyGoal->(wp39-25),PrecisionGoal->(wp39-25)]];
ls[e_?NumericQ] := lm[cp[e],e];
lsp[e_?NumericQ] := 1/(1-e)+cp[e] Exp[(e-1)cp[e]]/(1+Exp[(e-1)cp[e]]);
ii[e_] := e yp[e]-lm[mu,e];
dd[e_] := p/8 (yp[e]-mu)^2 e (1-e);
de[e_] := D[dd[x],x]/.x->e;
ye[e_] := D[yp[x],x]/.x->e;
rateRow[eps_,cc_] := Module[{e=N[cc eps,wp39+20], al=1+eps,y,th,is,j,jp,hh,del},
 y=yp[e]; th= z/.FindRoot[lsp[z]==y,{z,e-e^2/10},WorkingPrecision->(wp39-10),AccuracyGoal->(wp39-35),PrecisionGoal->(wp39-35)];
 is=th y-ls[th]; j=ii[e]-al dd[e]; jp=e-al de[e]/ye[e];
 hh=j-y jp+ls[jp]; del=j-is;
 <|"eps"->N[eps,18],"c"->cc,"eta"->N[e,18],"thetaStar"->N[th,18],"Jprime"->N[jp,18],"delta"->N[del,24],"H"->N[hh,24],
 "deltaOverEps4"->N[del/eps^4,18],"HOverEps4"->N[hh/eps^4,18]|>
];
If[!TrueQ[$rate39LibraryMode],
A3=mu (mu-1)/8; K4=mu (mu-1)(3+5mu)/384; crit=A3/K4;
rows=Flatten[Table[rateRow[ep,cc],{ep,{1/100,1/200,1/500,1/1000}},{cc,{5,6,7}}],1];
Export[out39,<|"wolfram"->$Version,"working_precision"->wp39,"A3"->N[A3,30],"K4"->N[K4,30],"criticalC"->N[crit,30],"rows"->rows|>,"RawJSON"];
Print[<|"output"->out39,"working_precision"->wp39,"cases"->Length[rows]|>];

];
