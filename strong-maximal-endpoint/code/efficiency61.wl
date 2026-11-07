(* r061: actual mixed Phi_N cost and Gamma--Bernoulli quantile diagnostic.
   The rows are floating diagnostics, not interval certificates. *)

ClearAll[logPhi, mesh, mixedRatio, lowerCDF, densityS, quantile, oneRow];
args61=$ScriptCommandLine;
wp=If[Length[args61]>=3,FromDigits[args61[[3]]],70];
dims61=If[Length[args61]>=4,{FromDigits[args61[[4]]]},{128,256,512,1024}];
kappa = 1/4;
mu = x /. FindRoot[x == 1 + Exp[-x], {x, 6/5}, WorkingPrecision -> wp,
    AccuracyGoal -> wp-15, PrecisionGoal -> wp-15];

logPhi[n_Integer, x_?NumericQ] := If[x <= 1, Log[x],
  Log[x] + (n - 1) Log[1 + Log[x]]];

mesh[n_Integer] := Module[{ll, qq, rr, mm, aa, qg, radii, hh, pp, dd, low, high,
   eta, bb, t},
 ll = Log[N[n, wp]]; qq = kappa ll^2; rr = N[qq/n, wp];
 mm = Max[3, Floor[N[n, wp]/ll^5]]; aa = 2 rr; qg = 1 + n rr/aa;
 radii = Join[{0}, Table[qg^(j - 1), {j, 1, mm}]];
 hh = Differences[radii];
 pp = Join[{1 - (mm - 1) aa}, ConstantArray[aa, mm - 1]];
 dd = rr radii[[2 ;;]]/pp;
 low = Table[1 - dd[[i]]/hh[[i]], {i, mm}];
 high = Table[1 + dd[[i]]/hh[[i + 1]], {i, mm - 1}];
 t = 1 - First[pp]; eta = rr t/(First[pp] - rr);
 bb = (1 + eta) (1 + 1/n^2);
 <|"N" -> n, "L" -> ll, "q" -> qq, "R" -> rr, "M" -> mm, "Q" -> qg,
   "r" -> radii, "h" -> hh, "p" -> pp, "low" -> low, "high" -> high,
   "eta" -> eta, "b" -> bb, "gamma" -> (1 + 1/n^2) (1 - rr)|>
 ];

(* Average over the full two-dimensional mesh support, after multiplying it
   by b.  It is divided by Phi_N(A), the cost of the bare core terminal. *)
mixedRatio[x_Association, aa_?NumericQ] := Module[{n = x["N"], m = x["M"], p = x["p"],
  h = x["h"], low = x["low"], high = x["high"], b = x["b"], rm, den},
 rm = Last[x["r"]]; den = logPhi[n, aa];
 Sum[p[[i]] Sum[h[[j]] Exp[logPhi[n, aa b Which[i == j, low[[i]],
       i < m && i + 1 == j, high[[i]], True, 1]] - den], {j, m}], {i, m}]/rm
 ];

(* S_n=Gamma(n,1)+mu Binomial(n,1/(1+Exp[mu])).  The formulas below are
   exact finite mixture identities; their numerical evaluation is Wolfram's. *)
lowerCDF[n_Integer, z_?NumericQ] := Module[{pp = 1/(1 + Exp[mu])},
 Sum[PDF[BinomialDistribution[n, pp], j]
   If[z - mu j <= 0, 0, CDF[GammaDistribution[n, 1], z - mu j]], {j, 0, n}]
 ];
densityS[n_Integer, z_?NumericQ] := Module[{pp = 1/(1 + Exp[mu])},
 Sum[PDF[BinomialDistribution[n, pp], j]
   If[z - mu j <= 0, 0, PDF[GammaDistribution[n, 1], z - mu j]], {j, 0, n}]
 ];
quantile[n_Integer, target_?NumericQ, start_?NumericQ] :=
 z /. FindRoot[lowerCDF[n, z] == target, {z, start}, WorkingPrecision -> wp,
   AccuracyGoal -> wp-32, PrecisionGoal -> wp-32, MaxIterations -> 100];

oneRow[n_Integer] := Module[{x, ll, a, gam, q0, qg, shift, dens, ratio, c, unchanged,
  upper, lower, res0, resg},
 x = mesh[n]; ll = x["L"]; a = 1 + Sqrt[ll/n]; gam = x["gamma"];
 q0 = quantile[n - 2, 1 - 1/a, (n - 2) mu - Sqrt[(n - 2) ll]];
 qg = quantile[n - 2, 1 - 1/(a gam), q0 - Sqrt[ll]];
 shift = q0 - qg; dens = densityS[n - 2, (q0 + qg)/2];
 ratio = mixedRatio[x, a]; c = x["b"] (1 + 1/n);
 unchanged = 1 - (First[x["p"]]/Last[x["r"]] + 2 x["R"] (1 - 1/Last[x["r"]])) -
   (First[x["p"]] (x["Q"] - 1)/x["Q"]^(x["M"] - 1) +
     2 x["R"] (1 - x["Q"]/x["Q"]^(x["M"] - 1)));
 upper = x["b"] (1 - x["R"])
   Exp[(n - 1) Log[(1 + Log[a c])/(1 + Log[a])]];
 lower = unchanged Exp[logPhi[n, a x["b"]] - logPhi[n, a]];
 res0 = lowerCDF[n - 2, q0] - (1 - 1/a);
 resg = lowerCDF[n - 2, qg] - (1 - 1/(a gam));
 <|"N_total" -> n, "core_dimension" -> n - 2, "mu" -> N[mu, 28],
   "L" -> N[ll, 28], "R" -> N[x["R"], 28], "M" -> x["M"],
   "eta" -> N[x["eta"], 20], "gamma_common_action" -> N[gam, 28],
   "A_core" -> N[a, 28], "q_critical" -> N[q0, 28],
   "q_with_common_mesh" -> N[qg, 28], "quantile_shift" -> N[shift, 28],
   "shift_over_logN" -> N[shift/ll, 28], "predicted_kappa_sqrt_mu" -> N[kappa Sqrt[mu], 28],
   "midpoint_density" -> N[dens, 28], "N_density_over_logN" -> N[n dens/ll, 28],
   "predicted_1_over_sqrt_mu" -> N[1/Sqrt[mu], 28],
   "mixed_Phi_ratio" -> N[ratio, 28], "unchanged_area" -> N[unchanged, 28],
   "analytic_ratio_lower" -> N[lower, 28], "analytic_ratio_upper" -> N[upper, 28],
   "critical_CDF_residual" -> ToString[N[res0, 12], InputForm],
   "mesh_CDF_residual" -> ToString[N[resg, 12], InputForm]|>
 ];

rows = oneRow /@ dims61;
out = <|"wolfram" -> $Version, "working_precision"->wp, "status" ->
 "actual finite Gamma--Binomial mixture and cellwise Phi_N evaluation; floating diagnostic, not an interval certificate",
 "formula_scope" -> "distinguished common mesh receiver only; no assertion about the full product halo or general C_N",
 "rows" -> rows|>;
out61=If[Length[args61]>=2,args61[[2]],"/tmp/efficiency61.json"];
If[!AllTrue[rows,#["analytic_ratio_lower"]<=#["mixed_Phi_ratio"]<=#["analytic_ratio_upper"]&],Exit[2]];
Export[out61, out, "RawJSON"];
Print[out];
