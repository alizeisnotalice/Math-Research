(* Actual numerical probe of the all-time, centered quadratic Poisson formula.
   It tests U_t'(r)/R_t(r); U' <= R is sufficient for (v+U)' <= 0,
   since -v'=R+Q/N.  No MATLAB is used. *)
$MaxExtraPrecision = 20000;
mu = N[x /. FindRoot[x == 1 + Exp[-x], {x, 1.3}, WorkingPrecision -> 80], 70];
p = N[1/(1 + Exp[mu]), 70];

mom[tt_Integer, rr_?NumericQ, pow_Integer, wt_: (1 &)] := Module[{kk, uu},
  kk = Range[0, Min[tt, Floor[rr/mu]]];
  If[Length[kk] == 0, Return[0``60]];
  uu = rr - mu kk;
  Total[SetPrecision[Binomial[tt, #] wt[#] uu[[# + 1]]^pow, 65] & /@ kk]
];

probe[tt_Integer, z_?NumericQ] := Module[
 {rr, x, aa, a1, cn, cnp, cc, ccp, hh, ffprime, ffactor, fprimeOverF,
  fprimeUOverF, roverf, uprimeOverf, ratio},
 rr = N[tt mu + z Sqrt[tt], 60]; x = rr - tt mu;
 aa = mom[tt, rr, tt - 1];
 a1 = If[tt >= 2, mom[tt, rr, tt - 2], 0];
 cn = mom[tt, rr, tt + 1]/(tt + 1) + mu^2/(tt - 1) mom[tt, rr, tt - 1, Function[k, k (tt-k)]];
 cnp = mom[tt, rr, tt] + mu^2 mom[tt, rr, tt - 2, Function[k, k (tt-k)]];
 cc = cn/aa; ccp = (cnp aa - cn (tt - 1) a1)/aa^2;
 hh = HarmonicNumber[tt] - 1;
 fprimeOverF = -1 + (tt - 1) a1/aa;
 ffactor = p/8 (x^2/tt + cc hh/tt);
 fprimeUOverF = p/8 (2 x/tt + ccp hh/tt);
 uprimeOverf = fprimeOverF ffactor + fprimeUOverF;
 roverf = mom[tt, rr, tt]/(tt aa);
 ratio = uprimeOverf/roverf;
 <|"t" -> tt, "z" -> z, "x" -> N[x, 12], "C/t" -> N[cc/tt, 12],
   "UprimeOverf" -> N[uprimeOverf, 12], "ROverf" -> N[roverf, 12],
   "UprimeOverR" -> N[ratio, 12]|>
];

tests = Flatten[Table[probe[tt, z], {tt, {64, 256, 1024}}, {z, {-6, -4, -2, 0, 2, 4, 6}}], 1];
Export["/tmp/r035_raw_poisson_monotonicity.json", tests, "RawJSON"];
Print[ExportString[tests, "Table"]];
