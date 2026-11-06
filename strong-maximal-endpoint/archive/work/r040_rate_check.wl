(* r040 nonlinear rate uplift check; numerical diagnostics only. *)
$MaxExtraPrecision=100000;
mu=x/.FindRoot[x==1+Exp[-x],{x,1.3},WorkingPrecision->100];
p=1/(1+Exp[mu]);
lm[c_,e_]:=Log[(1+Exp[(e-1)c])/((1-e)(1+Exp[-c]))];
yp[e_]:=1/(1-e)+mu Exp[(e-1)mu]/(1+Exp[(e-1)mu]);
cp[e_?NumericQ]:=c/.FindRoot[Exp[(1-e)c]-(1-e)Exp[c]+e==0,{c,mu+mu e/2},WorkingPrecision->80,AccuracyGoal->55,PrecisionGoal->55];
ls[e_?NumericQ]:=lm[cp[e],e];
lsp[e_?NumericQ]:=1/(1-e)+cp[e] Exp[(e-1)cp[e]]/(1+Exp[(e-1)cp[e]]);
imu[e_]:=e yp[e]-lm[mu,e];
thst[e_?NumericQ]:=z/.FindRoot[lsp[z]==yp[e],{z,e-e^2/10},WorkingPrecision->70,AccuracyGoal->45,PrecisionGoal->45];
istar[e_?NumericQ]:=Module[{y=yp[e],th=thst[e]},th y-ls[th]];
dd[e_]:=p/8 (yp[e]-mu)^2 e(1-e);
de[e_]:=D[dd[x],x]/.x->e;
ye[e_]:=D[yp[x],x]/.x->e;
row[eps_,e_]:=Module[{al=1+eps,y,j,jp,hh,ee,ff},
 y=yp[e]; ee=imu[e]-istar[e]; ff=ee-dd[e];
 j=istar[e]-eps dd[e]; jp=thst[e]-eps de[e]/ye[e];
 (* Directly use envelope derivative at J' in H. *)
 hh=j-y jp+ls[jp];
 <|"eps"->N[eps,14],"eta"->N[e,14],"E"->N[ee,20],"D"->N[dd[e],20],"F"->N[ff,20],"Jprime"->N[jp,16],"Hnew"->N[hh,20]|>
];
A3=mu (mu-1)/8;K4=mu (mu-1)(3+5mu)/384;
ans=Flatten[Table[row[ep,c ep],{ep,SetPrecision[{.01,.005,.002,.001},60]},{c,{3,6}}],1];
Export["/tmp/r040_rate_check.json",<|"A3"->N[A3,30],"K4"->N[K4,30],"rows"->ans|>,"RawJSON"];
Print[ExportString[<|"A3"->N[A3,20],"K4"->N[K4,20],"rows"->ans|>,"RawJSON"]];
