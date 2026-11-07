(* r061/block61: exact coordinate accounting for the far angular-band test.
   This is a deterministic Wolfram diagnostic.  It does not claim a weak
   theorem. *)

ClearAll[lowerCDF, quantile, logPhi, row];
wp = 80;
mu = x /. FindRoot[x == 1 + Exp[-x], {x, 6/5}, WorkingPrecision -> wp,
    AccuracyGoal -> 60, PrecisionGoal -> 60];

(* S is the sum of N independent Z_mu increments.  The finite
   Gamma--Binomial mixture is exact. *)
lowerCDF[n_Integer, z_?NumericQ] := Module[{p = 1/(1 + Exp[mu])},
  Sum[PDF[BinomialDistribution[n, p], j]
    If[z - mu j <= 0, 0, CDF[GammaDistribution[n, 1], z - mu j]],
   {j, 0, n}]
  ];

quantile[n_Integer, alpha_?NumericQ, start_?NumericQ] :=
 z /. FindRoot[lowerCDF[n, z] == alpha, {z, start},
   WorkingPrecision -> wp, AccuracyGoal -> 48, PrecisionGoal -> 48,
   MaxIterations -> 100];

logPhi[n_Integer, a_?NumericQ] :=
 Log[a] + (n - 1) Log[1 + Log[a]];

row[n_Integer, ell_Integer] := Module[
 {a = SetPrecision[3/2, wp], q, w, r, phi, pTailCheb, avgLower,
  globalJ, receiverV, scoreFactor, crudeTranslatedUpper,
  centralTranslatedExact, crudeOverExact, cdfResidual},
 q = quantile[n, SetPrecision[1/4, wp], n mu - Sqrt[n]];
 (* Var(T_1-T_2)=n mu.  Hence P(|D|>W)<=n mu/W^2 exactly by Chebyshev. *)
 w = 8 Sqrt[n Log[n]];
 If[ell <= Max[w, q], Return[<|"N" -> n, "error" -> "Require L>Max[W,q] for the chosen exact local-budget branch."|>]];
 r = ell + q;
 phi = Exp[logPhi[n, a]];
 pTailCheb = n mu/w^2;
 avgLower = a (SetPrecision[3/4, wp] - pTailCheb);
 (* H=a 1_{t1+t2>=r, |t1-t2|<=ell}; r>=ell makes this exact. *)
 globalJ = ell phi Exp[-r];
 (* E={r-q<=u1+u2<=r-q+1, |u1-u2|<=ell-w}; its points are received
    whenever avgLower>1. *)
 receiverV = (ell - w) Exp[-(r - q)] (1 - Exp[-1]);
 scoreFactor = receiverV/globalJ;
 (* The usual global-to-local upper bound is exp(u1+u2) J(H).  At the
    central receiver its exact upper-cone modular budget is the integral
    over s>=0 of Phi(H(u+s)) exp(-s1-s2). *)
 crudeTranslatedUpper = ell phi Exp[-q];
 centralTranslatedExact = phi ((q + 1) Exp[-q] - Exp[-ell]);
 crudeOverExact = crudeTranslatedUpper/centralTranslatedExact;
 cdfResidual = lowerCDF[n, q] - SetPrecision[1/4, wp];
 <|"N" -> n, "mu" -> N[mu, 30], "amplitude" -> N[a, 30],
   "sum_quantile_q_P_S_ge_q_3_4" -> N[q, 30],
   "W" -> N[w, 30], "angular_halfwidth_L" -> ell,
   "radial_threshold_R" -> N[r, 30],
   "chebyshev_D_tail_upper" -> N[pTailCheb, 30],
   "received_average_lower_bound" -> N[avgLower, 30],
   "global_modular_J_exact" -> ToString[N[globalJ, 30], InputForm],
   "receiver_weight_exact" -> ToString[N[receiverV, 30], InputForm],
   "receiver_over_J_exact" -> N[scoreFactor, 30],
   "crude_exp_sum_u_times_global_J" -> ToString[N[crudeTranslatedUpper, 30], InputForm],
   "central_actual_translated_budget_exact" -> ToString[N[centralTranslatedExact, 30], InputForm],
   "crude_over_exact_ratio" -> N[crudeOverExact, 30],
   "quantile_CDF_residual" -> ToString[N[cdfResidual, 12], InputForm]|>
 ];

out = <|"wolfram" -> $Version,
 "status" -> "Exact Gamma--Binomial quantile plus closed-form band accounting; Chebyshev is a rigorous D-tail bound.",
 "scope" -> "Shows only that the bound J_u<=exp(sum u) J(H) is too coarse for a global r046 localization argument; it neither disproves nor proves a global weak theorem.",
 "rows" -> {row[64, 10^8], row[128, 10^12]}|>;
out61=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/block61.json"];
Export[out61, out, "RawJSON"];
Print[out];
