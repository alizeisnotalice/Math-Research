(* Exact symbolic companion check for round 10.  Run with wolframscript -file. *)
Clear[m, u, c];
equalGap = FullSimplify[
  ((1 - 1/m) + u) (1 - 1/m)^(m - 1) - (1 - 1/m)^m,
  Assumptions -> Element[m, Integers] && m >= 2 && u > 0
];
unequalGap = FullSimplify[(1/2 + u) (1/2) - 1/4, Assumptions -> u > 0];
waterFill = FullSimplify[
  (1 - 1/2)^1 (1 - 1/2)^1 (1 - 0)^2
];
Print[{"equal-gap", equalGap}];
Print[{"unequal-gap", unequalGap}];
Print[{"a=(1,1,2)-candidate-value", waterFill}];
