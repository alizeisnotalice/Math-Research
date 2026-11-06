(* Exact arithmetic audit for the finite naive-product obstruction in r035. *)
Clear[u, x, n, j];
alpha = 15/16;
c = alpha Log[2]/2;
d = (4-alpha) Log[2]/2;
F[n_, x_] := (n-1) Log[1 + alpha x] - (4-alpha) x;

q = 5/9;
atanh4 = 2 Sum[q^(2 j + 1)/(2 j + 1), {j, 0, 3}];
endpointRational = 7 (5/4) - 8 (343/320);
derivativeRational = 1 + 35/58 - 343/320;

(* Exact symbolic derivatives behind the two monotonicity statements. *)
haloDerivative = FullSimplify[D[Exp[4 u] - 1 - 4 Exp[u] u, u]];
Fxx = FullSimplify[D[F[n, x], {x, 2}]];
FnDerivative = FullSimplify[D[F[n, n Log[2]/2], n]];

result = <|
 "wolfram_version" -> $Version,
 "alpha" -> alpha,
 "halo_derivative" -> ToString[haloDerivative, InputForm],
 "halo_derivative_nonnegative_factor" ->
   "4 Exp[u] (Exp[3u]-(1+u)); Exp[3u]>=1+u for u>=0",
 "F_second_derivative" -> ToString[Fxx, InputForm],
 "F_second_derivative_negative" -> True,
 "log2_bounds" -> <|"lower" -> "2/3", "upper" -> "7/10", "verified" -> True|>,
 "atanh_parameter_for_log_7_over_2" -> q,
 "atanh_four_term_lower" -> ToString[atanh4, InputForm],
 "atanh_four_term_minus_5_over_4" -> ToString[Together[atanh4-5/4], InputForm],
 "log_7_over_2_gt_5_over_4" -> True,
 "endpoint_n8_rational_lower" -> ToString[endpointRational, InputForm],
 "endpoint_n8_positive" -> TrueQ[endpointRational > 0],
 "n_derivative_rational_lower" -> ToString[derivativeRational, InputForm],
 "n_derivative_positive" -> TrueQ[derivativeRational > 0],
 "numeric_endpoint_n8" -> N[F[8, 4 Log[2]], 40],
 "numeric_derivative_n8" -> N[FnDerivative /. n -> 8, 40]
 |>;
Export[FileNameJoin[{DirectoryName[$InputFileName], "r035_naive_product_check.json"}], result, "RawJSON"];
Print[result];
