(* Exact Wolfram-language verification for the displayed cube formula. *)
Do[
  q = Rational[n + 2, 3 n];
  Do[
    s2 = 1 + delta^2 q;
    Print[{"n", n, "delta", delta, "q", q,
      "SindividualSquared", s2, "strictGap", FullSimplify[s2 > 1]}];
    Do[
      deff = delta (cp - cm)/(cp + cm);
      grouped = 1 + deff^2 q;
      If[!TrueQ[FullSimplify[1 <= grouped <= s2]], Abort[]],
      {cp, {1, 2, 3}}, {cm, {1, 2, 5}}],
    {delta, {1/2, 3/4}}],
  {n, 1, 3}];
Print["PASS: exact cube coefficients and grouping inequalities"];
