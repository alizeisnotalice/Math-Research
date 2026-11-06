(* r052 exact symbolic corroboration; boundary/soft-support arguments are in proof.md. *)
Clear[r,t,d,ph,hh,lam,rr,xx];
translationResidual52=FullSimplify[(D[Exp[d]r[t-d],d]/.d->0)-(r[t]-r'[t])];
fenchelResidual52=FullSimplify[(D[ph[hh[xx]]-lam rr[xx]hh[xx],xx]/.ph'[hh[xx]]->lam rr[xx])+lam hh[xx]rr'[xx]];
Clear[t, c, th, kk, dens, f1, f2, g1, g2, h1, h2, p, qDirect, qFormula,
  s1, s2, k1, k2, bracket, samples];

dens[bb_, jj_, x_, cc_] := Exp[-x]/(Gamma[bb] (1 + Exp[-cc])^bb) *
  Sum[Binomial[bb, kk] (x - kk cc)^(bb - 1), {kk, 0, jj}];

(* On the open cell  c<t1<2c, 0<t2<c.  The cell indices remain fixed for
   sufficiently small contrast perturbations c+th,c-th. *)
f1 = dens[4, 1, t, c + th];
f2 = dens[4, 0, 3 t/10, c - th];
p = 1/(1 + Exp[c]);
g1 = dens[3, 0, t - c, c]; g2 = 0;
h1 = 0; h2 = 0;
s1 = 4 p (1 - g1/(f1 /. th -> 0));
s2 = 4 p (1 - g2/(f2 /. th -> 0));
k1 = FullSimplify[D[Log[Exp[t] (f1 /. th -> 0)], t]];
k2 = FullSimplify[D[Log[Exp[3 t/10] (f2 /. th -> 0)], t]/(3/10)];
qDirect = FullSimplify[D[f1 f2, {th, 2}]/(f1 f2) /. th -> 0];
qFormula = FullSimplify[4 p (-2 (1-p) + g1/(f1 /. th -> 0) +
    g2/(f2 /. th -> 0) + 3 p (h1/(f1 /. th -> 0) + h2/(f2 /. th -> 0)) -
    8 p (g1/(f1 /. th -> 0)) (g2/(f2 /. th -> 0)) +
    D[g1, t]/(f1 /. th -> 0))];


localPass52=TrueQ[FullSimplify[qDirect-qFormula,Assumptions->c>0&&c<t<2c]===0];
Clear[b, c, T, eps, kk, dens, tailRow, knotRow, rows];
dens[bb_, jj_, x_, cc_] := Exp[-x]/(Gamma[bb] (1 + Exp[-cc])^bb) *
  Sum[Binomial[bb, kk] (x - kk cc)^(bb - 1), {kk, 0, jj}];

tailRow[bb_] := Module[{p, f, g, r, k, coeff, target},
 p = 1/(1 + Exp[c]);
 f = dens[bb, bb, T, c];
 g = dens[bb - 1, bb - 1, T - c, c];
 r = FullSimplify[g/f];
 k = FullSimplify[D[Log[Exp[T] f], T]];
 coeff = FullSimplify[Limit[T^2 (bb p r)^2 - T^2 bb p k r,
   T -> Infinity, Assumptions -> c > 0]];
 target = bb (bb - 1)^2 (bb - 2)/4;
 <|"b" -> bb, "tail_positive_limit" -> TrueQ[FullSimplify[coeff == target,
      Assumptions -> c > 0]], "tail_coefficient" -> ToString[InputForm[coeff]]|>
];

knotRow[bb_] := Module[{p, f, g, r, k, coeff, target},
 p = 1/(1 + Exp[c]);
 f = dens[bb, 1, c + eps, c];
 g = dens[bb - 1, 0, eps, c];
 r = FullSimplify[g/f];
 k = FullSimplify[D[Log[Exp[T] dens[bb, 1, T, c]], T] /. T -> c + eps];
 coeff = FullSimplify[Limit[((bb p r)^2 - bb p k r)/eps^(bb - 2), eps -> 0,
   Direction -> "FromAbove", Assumptions -> c > 0]];
 target = -bb (bb - 1)^2/c^bb;
 <|"b" -> bb, "knot_negative_limit" -> TrueQ[FullSimplify[coeff == target,
      Assumptions -> c > 0]], "knot_coefficient" -> ToString[InputForm[coeff]]|>
];

rows = Table[Join[tailRow[b], KeyDrop[knotRow[b], "b"]], {b, 4, 8}];

result52=<|"wolfram"->$Version,"translation_derivative_pass"->TrueQ[translationResidual52===0],"Fenchel_derivative_pass"->TrueQ[fenchelResidual52===0],"contrast_density_b4_open_cell_pass"->localPass52,"two_sign_limit_checks"->rows,"all_pass"->(TrueQ[translationResidual52===0]&&TrueQ[fenchelResidual52===0]&&localPass52&&And@@Flatten[Lookup[rows,{"tail_positive_limit","knot_negative_limit"}]]),"scope"->"symbolic corroboration of exact algebra/limits; general boundary and actual-soft-support constructions are proved separately"|>;
out52=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/identities52.json"];Export[out52,result52,"RawJSON"];Print[KeyTake[result52,{"all_pass","translation_derivative_pass","Fenchel_derivative_pass","contrast_density_b4_open_cell_pass"}]];
