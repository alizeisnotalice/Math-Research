(* r039 finite-rate checks for the proposed intermediate-time obstruction.
   WolframScript calculation only; no interval claim. *)
$MaxExtraPrecision = 100000;
mu = x /. FindRoot[x == 1 + Exp[-x], {x, 1.3}, WorkingPrecision -> 100];
p = 1/(1+Exp[mu]);
lm[c_,e_] := Log[(1+Exp[(e-1)c])/((1-e)(1+Exp[-c]))];
yp[e_] := 1/(1-e)+mu Exp[(e-1)mu]/(1+Exp[(e-1)mu]);
cp[e_?NumericQ] := c /. FindRoot[Exp[(1-e)c]-(1-e)Exp[c]+e==0,{c,mu+mu e/2},WorkingPrecision->80,AccuracyGoal->55,PrecisionGoal->55];
ls[e_?NumericQ] := lm[cp[e],e];
lsp[e_?NumericQ] := 1/(1-e)+cp[e] Exp[(e-1)cp[e]]/(1+Exp[(e-1)cp[e]]);
ii[e_] := e yp[e]-lm[mu,e];
dd[e_] := p/8 (yp[e]-mu)^2 e (1-e);
de[e_] := D[dd[x],x]/.x->e;
ye[e_] := D[yp[x],x]/.x->e;
rateRow[eps_,cc_] := Module[{e=N[cc eps,70], al=1+eps,y,th,is,j,jp,hh,del},
 y=yp[e]; th= z/.FindRoot[lsp[z]==y,{z,e-e^2/10},WorkingPrecision->70,AccuracyGoal->45,PrecisionGoal->45];
 is=th y-ls[th]; j=ii[e]-al dd[e]; jp=e-al de[e]/ye[e];
 hh=j-y jp+ls[jp]; del=j-is;
 <|"eps"->N[eps,18],"c"->cc,"eta"->N[e,18],"thetaStar"->N[th,18],"Jprime"->N[jp,18],"delta"->N[del,24],"H"->N[hh,24],
 "deltaOverEps4"->N[del/eps^4,18],"HOverEps4"->N[hh/eps^4,18]|>
];
A3=mu (mu-1)/8; K4=mu (mu-1)(3+5mu)/384; crit=A3/K4;
rows=Flatten[Table[rateRow[SetPrecision[ep,60],cc],{ep,{1/100,1/200,1/500,1/1000}},{cc,{5,6,7}}],1];
Export["/tmp/r039_global_inflation_check.json",<|"A3"->N[A3,30],"K4"->N[K4,30],"criticalC"->N[crit,30],"rows"->rows|>,"RawJSON"];
Print[ExportString[<|"A3"->N[A3,20],"K4"->N[K4,20],"criticalC"->N[crit,20],"rows"->rows|>,"RawJSON"]];
