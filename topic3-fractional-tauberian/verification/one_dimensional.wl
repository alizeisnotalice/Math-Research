(* Symbolic checks used by the one-dimensional section.  Wolfram Language 15. *)
Clear[t, beta, lambda, ell, p1, p2, p3, p4];

inner[t_, beta_] := 1 - (1 - t)^(1/beta);
outer[t_, beta_] := (1 - t)^(-1/(1 - beta)) - 1;

Print["inner Taylor: ", Series[inner[t, beta], {t, 0, 4}]];
Print["outer Taylor: ", Series[outer[t, beta], {t, 0, 4}]];

(* Solve the transition equation in ell=-Log[1-t]. *)
b = 1/2 + p1 ell + p2 ell^2 + p3 ell^3 + p4 ell^4;
transition = Exp[-ell/b] + Exp[ell/(1 - b)] - 2;
series = Normal[Series[transition, {ell, 0, 6}]];
rules = First@Solve[
   Table[Coefficient[series, ell, j] == 0, {j, 2, 5}],
   {p1, p2, p3, p4}
];
Print["transition coefficients in ell: ", rules];
Print["transition coefficients in t: ",
  Series[(b /. rules) /. ell -> -Log[1 - t], {t, 0, 4}]
];

criticalBeta = 1/2 + lambda t;
Print["critical inner: ", Series[inner[t, criticalBeta], {t, 0, 5}]];
Print["critical outer: ", Series[outer[t, criticalBeta], {t, 0, 5}]];
Print["critical difference outer-inner: ",
  Factor@Normal@Series[outer[t, criticalBeta] - inner[t, criticalBeta], {t, 0, 5}]
];
