(* Symbolic checks; assumptions and mathematical proof are in the manuscript. *)
ClearAll[b, c, t, a, s, z];
Print["Wolfram: ", $Version];
forms = {
  (1 - (1 - t)^(1/b))^2,
  2 - 4 (1 - t)/(2 - b) + 2 b (1 - t)^(2/b)/(2 - b),
  4 (c (1 - t)^(-1/c)/(c + 1) - 1 + (1 - t)/(c + 1)),
  ((1 - t)^(-1/c) - 1)^2
};
expected = {
  t^2/b^2 - (1 - b) t^3/b^3,
  2 t^2/b - 4 (1 - b) t^3/(3 b^2),
  2 t^2/c + 2 (1 + 2 c) t^3/(3 c^2),
  t^2/c^2 + (1 + c) t^3/c^3
};
residuals = Table[FullSimplify[Normal[Series[forms[[j]], {t, 0, 3}]] - expected[[j]],
  Assumptions -> 0 < b < 2 && 0 < c < 2], {j, 4}];
Print["Cubic coefficient residuals: ", residuals];
If[residuals =!= {0, 0, 0, 0}, Exit[1]];
(* Ordered-chamber integrals for the B and C active branches. *)
resB = FullSimplify[8 Integrate[1 - 2 s - a (1 - 2 s)^(1 - b),
  {s, 0, (1 - a^(1/b))/2}, Assumptions -> 0 < a < 1 && 1/2 <= b < 1]
  - (2 - 4 a/(2 - b) + 2 b a^(2/b)/(2 - b)),
  Assumptions -> 0 < a < 1 && 1/2 <= b < 1];
resC = FullSimplify[8 Integrate[1 - a (1 + 2 s)^c,
  {s, 0, (a^(-1/c) - 1)/2}, Assumptions -> 0 < a < 1 && 1/2 < c <= 1]
  - 4 (c a^(-1/c)/(c + 1) - 1 + a/(c + 1)),
  Assumptions -> 0 < a < 1 && 1/2 < c <= 1];
Print["Area integral residuals: ", {resB, resC}];
If[{resB, resC} =!= {0, 0}, Exit[1]];
Print["All symbolic checks passed."];
