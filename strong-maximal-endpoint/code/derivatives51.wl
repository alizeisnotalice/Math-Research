(* r051 reproducible symbolic check: second c-derivatives on every open moving-knot cell. *)
(* Keep iterator symbols fresh and retain the explicit multiplication across
   line breaks in the two finite-sum formulas. *)
Clear[t, c, kk, rr, densCell, cdfCell, checkDensity, checkCDF];
densCell[bb_, jj_, xx_, cc_] := Exp[-xx]/(Gamma[bb] (1 + Exp[-cc])^bb) *
  Sum[Binomial[bb, kk] (xx-kk cc)^(bb-1), {kk, 0, jj}];
cdfCell[bb_, jj_, xx_, cc_] := 1/(1 + Exp[-cc])^bb *
  Sum[Binomial[bb, kk] (Exp[-kk cc] - Exp[-xx]
    Sum[(xx-kk cc)^rr/Factorial[rr], {rr, 0, bb-1}]), {kk, 0, jj}];

checkDensity[b_, j_] := Module[{f, g, h, rhs, residual, asm},
 f = densCell[b, j, t, c];
 g = If[j >= 1, densCell[b-1, j-1, t-c, c], 0];
 h = If[j >= 2, densCell[b-2, j-2, t-2 c, c], 0];
 rhs = b/(1 + Exp[c]) (((b+1)/(1 + Exp[c])-1) f-
   (2 b/(1 + Exp[c])-1) g+(b-1)/(1 + Exp[c]) h+D[g, t]);
 asm = c > 0 && t > j c && If[j < b, t < (j+1)c, True];
 residual = FullSimplify[D[f, {c, 2}]-rhs, Assumptions -> asm];
 TrueQ[residual === 0]
];
checkCDF[b_, j_] := Module[{f, g, h, gd, rhs, residual, asm},
 f = cdfCell[b, j, t, c];
 g = If[j >= 1, cdfCell[b-1, j-1, t-c, c], 0];
 h = If[j >= 2, cdfCell[b-2, j-2, t-2 c, c], 0];
 gd = If[j >= 1, densCell[b-1, j-1, t-c, c], 0];
 rhs = b/(1 + Exp[c]) (((b+1)/(1 + Exp[c])-1) f-
   (2 b/(1 + Exp[c])-1) g+(b-1)/(1 + Exp[c]) h+gd);
 asm = c > 0 && t > j c && If[j < b, t < (j+1)c, True];
 residual = FullSimplify[D[f, {c, 2}]-rhs, Assumptions -> asm];
 TrueQ[residual === 0]
];

rows = Flatten[Table[
  <|"b" -> b, "cell" -> j, "domain" -> If[j < b,
    ToString[j] <> " c < t < " <> ToString[j+1] <> " c",
    ToString[b] <> " c < t"],
   "density_second_c_pass" -> checkDensity[b, j],
   "cdf_second_c_pass" -> checkCDF[b, j]|>,
  {b, 4, 8}, {j, 0, b}], 1];
allPass = And @@ Flatten[Lookup[rows, {"density_second_c_pass", "cdf_second_c_pass"}]];
result = <|"wolfram" -> $Version, "scope" -> "exact symbolic open moving-knot cells; b=4..8; c>0",
 "formula" -> "f''=bp(((b+1)p-1)f-(2bp-1)g+(b-1)ph+g_t); CDF replaces g_t by g",
 "rows" -> rows, "all_pass" -> allPass|>;
out = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]],
  "/tmp/derivatives51.json"];
Export[out, result, "RawJSON"];
Print[<|"all_pass" -> allPass, "cells" -> Length[rows], "output" -> out|>];
