(* Mathematica 15 symbolic cross-check.  Its output is evidence for the
   scalar identities only, not a proof of the analytic covering claims. *)
Print["Wolfram version: ", $Version];
Do[
  lstar = ((m - 1)/m)^m;
  s = m (1 - lstar^(1/m));
  phi = 1 - s;
  Print[{"m", m, "lambdaStar", lstar, "s", FullSimplify[s],
    "phi", FullSimplify[phi], "check", FullSimplify[{s == 1, phi == 0}]}],
  {m, {2, 3, 4, 8}}
];

Print["Series at lambda=1 (lambda=1-t):"];
Print[Assuming[m > 0,
  Series[m (1 - (1 - t)^(1/m)), {t, 0, 3}]]];

Print["Exact star-model algebra:"];
Print[FullSimplify[
  1 + 2^n (1 - eta) (1 - alpha)/alpha -
    (1 + 2^n (1 - alpha)/alpha),
  Assumptions -> n \[Element] Integers && n >= 1 && 0 < alpha < 1 && 0 < eta < 1]];
Print[FullSimplify[
  1/(1 + (1 - eta) (1 - alpha)/alpha) > alpha,
  Assumptions -> 0 < alpha < 1 && 0 < eta < 1]];
