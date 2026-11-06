(* Rate comparison for the nonlinear r036 candidate. Wolfram 15 actual scan;
   no claim of interval certification. *)
$MaxExtraPrecision = 100000;
mu = x /. FindRoot[x == 1 + Exp[-x], {x, 1.3}, WorkingPrecision -> 90];
pp = 1/(1 + Exp[mu]);
lm[c_, e_] := Log[(1 + Exp[(e-1)c])/((1-e) (1+Exp[-c]))];
lmprime[c_, e_] := 1/(1-e) + c Exp[(e-1)c]/(1+Exp[(e-1)c]);
cstar[e_?NumericQ] := c /. FindRoot[Exp[(1-e)c]-(1-e) Exp[c]+e == 0,
  {c, mu + mu e/2}, WorkingPrecision -> 70, PrecisionGoal -> 45, AccuracyGoal -> 45];
lst[e_?NumericQ] := lm[cstar[e],e];
lstprime[e_?NumericQ] := lmprime[cstar[e],e];
lmu[e_?NumericQ] := lm[mu,e];
lmuprime[e_?NumericQ] := lmprime[mu,e];
ym[x_] := 1/(1-x) + mu Exp[(x-1)mu]/(1+Exp[(x-1)mu]);
row[e_?NumericQ, al_:1] := Module[{y, th, im, is, dd, gap, de, ye, dee, yee, jprime, jpp, hh},
 y = lmuprime[e];
 th = z /. FindRoot[lstprime[z] == y, {z, Max[e/3, 10^-12], Min[.999, e]},
   WorkingPrecision -> 60, PrecisionGoal -> 36, AccuracyGoal -> 36];
 im = e y-lmu[e]; is = th y-lst[th];
 gap = im-is;
 dd = pp/8 (y-mu)^2 e (1-e);
 de = (pp/8) D[(ym[x]-mu)^2 x (1-x),x]/.x->e;
 ye = D[ym[x],x]/.x->e;
 dee = (pp/8) D[(ym[x]-mu)^2 x (1-x),{x,2}]/.x->e;
 yee = D[ym[x],{x,2}]/.x->e;
 jprime = e-al de/ye;
 jpp = (1-al (dee ye-de yee)/ye^2)/ye;
 hh = (im-al dd)-y jprime+lst[jprime];
 <|"eta"->N[e,18],"thetaStar"->N[th,18],"y"->N[y,18],
   "gap"->N[gap,18],"D"->N[dd,18],"gapOverD"->N[gap/dd,18],
   "alpha"->al,"Jprime"->N[jprime,18],"Jpp"->N[jpp,18],"H"->N[hh,18]|>
];
etas = Join[SetPrecision[{.0001,.0003,.001,.003,.01,.03,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,.98,.99},50],
  SetPrecision[Range[1,99]/100,50]];
ans = row /@ DeleteDuplicates[etas];
Export["/tmp/r036_nonlinear_rate_scan.json",ans,"RawJSON"];
Print[ExportString[ans,"Table"]];
Print["alpha tests at selected eta:"];
Print[Table[row[SetPrecision[e,50],al],{al,{1,1.01,1.05,1.1,1.25,1.5}},{e,{.01,.1,.4,.7,.9,.95,.98,.99}}]];

(* Formal power-series algebra: c*(e) and the Legendre-matched theta(e). *)
ord = 8;
cc = mu + Sum[Symbol["a"<>ToString[k]] ee^k,{k,1,ord}];
eqc = Normal@Series[Exp[(1-ee)cc]-(1-ee)Exp[cc]+ee,{ee,0,ord+1}];
solc = First@Solve[Table[Coefficient[eqc,ee,k]==0,{k,2,ord+1}], Table[Symbol["a"<>ToString[k]],{k,1,ord}]];
cs = cc/.solc;
ls = Normal@Series[lm[cs,ee],{ee,0,ord}];
lms = Normal@Series[lm[mu,ee],{ee,0,ord}];
thser = ee + Sum[Symbol["d"<>ToString[k]] ee^k,{k,2,ord-1}];
eqth = Normal@Series[(D[ls,ee]/.ee->thser) - D[lms,ee],{ee,0,ord-2}];
sold = First@Solve[Table[Coefficient[eqth,ee,k]==0,{k,1,ord-2}],Table[Symbol["d"<>ToString[k]],{k,2,ord-1}]];
ths=thser/.sold;
ys=D[lms,ee];
imser=Normal@Series[ee ys-lms,{ee,0,ord}];
isser=Normal@Series[ths ys-(ls/.ee->ths),{ee,0,ord}];
gapser=Normal@Series[imser-isser,{ee,0,ord}];
dser=Normal@Series[(pp/8)(ys-mu)^2 ee(1-ee),{ee,0,ord}];
Print["cstar series: ",InputForm[cs]];
Print["theta series: ",InputForm[ths]];
Print["gap series: ",InputForm[gapser]];
Print["D series: ",InputForm[dser]];
Print["gap-D series: ",InputForm[Normal@Series[gapser-dser,{ee,0,ord}]]];
