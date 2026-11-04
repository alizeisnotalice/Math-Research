(* Wolfram Language 15 cross-check for the first water-filling regime. *)
Print["Wolfram version: ", $Version];
Clear[a, b, k, s, eps];
rho = k a/(k - s);
Print["inactive threshold:", FullSimplify[rho <= b,
  Assumptions -> a > 0 && b > a && 0 <= s <= k (1 - a/b)]];
Print["first-regime product:", FullSimplify[
  (1 - s/k)^(k a), Assumptions -> a > 0 && k > 0 && 0 <= s < 1]];
Print["inverse series:", Series[k (1 - (1 - eps)^(1/(k a))), {eps, 0, 3}]];
