(* Symbolic check of the corner layer-cake formula; no differentiation at 0. *)
Print["Wolfram version: ", $Version];
Clear[cornerMass, n, tau, delta];
cornerMass[n_Integer, tau_, delta_] := Module[{mid},
  mid = Sum[
    Binomial[n, k] (1 - delta)^(n - k) *
      ((1 + delta)^(tau - n + k) - (1 - delta)^(tau - n + k)) /
        (tau - n + k),
    {k, 0, n}
  ];
  (n - tau) (2^n (1 - delta)^tau/tau + mid) +
    2^n (1 + delta)^(tau - n)
];

Do[
  Module[{dim = pair[[1]], exponent = pair[[2]], kappa, ratio, target},
    kappa = 2^dim dim/exponent;
    ratio = FullSimplify[cornerMass[dim, exponent, delta]/kappa];
    target = -exponent (dim - exponent)/2;
    Print["n=", dim, ", tau=", exponent, ", series: ", Series[ratio, {delta, 0, 3}]];
    Print["quadratic coefficient verified: ",
      FullSimplify[Coefficient[Normal@Series[ratio, {delta, 0, 3}], delta, 2] == target]]
  ],
  {pair, {{2, 1/2}, {3, 3/4}, {4, 1/2}}}
];

Clear[b, aa];
Print["cutoff compatibility:", Reduce[
  b/2 < aa < (1 - 2 b)/(n - n b) && 0 < b < 1/2 && n > 0,
  aa, Reals]];
