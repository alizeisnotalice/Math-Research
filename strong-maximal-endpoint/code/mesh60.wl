(* Exact finite checks for the non-uniform outer-width mesh in mesh60.md.
   No random sampling: every displayed prefix/suffix candidate is enumerated. *)

ClearAll[phi, data, fibrePrefix, suffixPrefix, modular, clean, oneRow];

phi[n_Integer][x_] := Piecewise[{{x, 0 <= x <= 1}},
  x (1 + Log[x])^(n - 1)];

data[n_Integer, m_Integer, rr_Rational] := Module[
  {aa = 2 rr, q, radii, h, p, dd, low, high},
  q = 1 + n rr/aa;
  radii = Join[{0}, Table[q^(j - 1), {j, 1, m}]];
  h = Differences[radii];
  p = Join[{1 - (m - 1) aa}, ConstantArray[aa, m - 1]];
  dd = rr/p;
  low = Table[1 - dd[[i]] radii[[i + 1]]/h[[i]], {i, m}];
  high = Table[1 + dd[[i]] radii[[i + 1]]/h[[i + 1]], {i, m - 1}];
  <|"N" -> n, "M" -> m, "R" -> rr, "alpha" -> aa, "Q" -> q,
    "r" -> radii, "h" -> h, "p" -> p, "d" -> dd,
    "low" -> low, "high" -> high|>
];

fibrePrefix[x_Association, i_Integer, j_Integer] := Module[
  {h = x["h"], r = x["r"][[j + 1]], m = x["M"], low = x["low"], high = x["high"]},
  Sum[h[[k]] Which[i == k, low[[i]], i < m && i + 1 == k, high[[i]],
      True, 1], {k, j}]/r
];

(* Average over the outer suffix l,...,M at inner grid prefix j. *)
suffixPrefix[x_Association, l_Integer, j_Integer] := Module[
  {p = x["p"]}, Sum[p[[i]] fibrePrefix[x, i, j], {i, l, x["M"]}]/
    Sum[p[[i]], {i, l, x["M"]}]
];

(* Normalized integral of Phi_n(scale*f), direct cell-by-cell evaluation. *)
modular[x_Association, scale_] := Module[
  {n = x["N"], m = x["M"], p = x["p"], h = x["h"], low = x["low"],
   high = x["high"], rm = Last[x["r"]]},
  Sum[p[[i]] Sum[h[[k]] phi[n][scale Which[
       i == k, low[[i]], i < m && i + 1 == k, high[[i]], True, 1]],
    {k, m}], {i, m}]/rm
];

oneRow[n_Integer, m_Integer, rr_Rational] := Module[
  {x = data[n, m, rr], pp, eta, eps, b, t, raw, scaled, prefixOK,
   fullOK, suffixOK, full, esc, adaptive, strong},
  pp = x["p"]; t = 1 - First[pp];
  eta = rr t/(First[pp] - rr); eps = 1/n^2; b = (1 + eta) (1 + eps);
  prefixOK = And @@ Flatten[Table[
     FullSimplify[fibrePrefix[x, i, j] == If[i == j, 1 - rr/pp[[i]], 1]],
     {i, m}, {j, m}]];
  fullOK = And @@ Table[FullSimplify[Sum[pp[[i]] fibrePrefix[x, i, j], {i, m}]
       == 1 - rr], {j, m}];
  suffixOK = And @@ Flatten[Table[
      FullSimplify[suffixPrefix[x, l, j] <= If[l == 1, 1 - rr, 1]], {l, m}, {j, m}]];
  full = (1 - rr)/(1 + eta); esc = t/(t + eta);
  adaptive = b/(1 + eta); strong = b Max[full, esc];
  raw = modular[x, 1]; scaled = modular[x, b];
  <|"N" -> n, "M" -> m, "R" -> ToString[rr, InputForm],
    "alpha" -> ToString[x["alpha"], InputForm], "Q" -> x["Q"],
    "p1" -> ToString[First[pp], InputForm], "tail" -> ToString[t, InputForm],
    "etaStar" -> ToString[eta, InputForm], "prefixTableOK" -> prefixOK,
    "weightedFullPrefixOK" -> fullOK, "allSuffixGridUpperOK" -> suffixOK,
    "crossingOK" -> FullSimplify[full == esc],
    "rawBudget" -> N[raw, 40], "scaledBudget" -> N[scaled, 40],
    "adaptiveScaled" -> N[adaptive, 40], "strongScaled" -> N[strong, 40],
    "strictSeparationOK" -> TrueQ[FullSimplify[adaptive > 1 && strong < 1]]|>
];

clean[x_Association] := AssociationMap[clean, x];
clean[x_List] := clean /@ x;
clean[x_Rational] := ToString[x, InputForm];
clean[x_?NumericQ] := N[x, 40];
clean[x_] := x;

rows = {oneRow[64, 4, 1/32], oneRow[256, 5, 1/64],
  oneRow[1024, 6, 1/128]};
out60=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/mesh60.json"];
If[!And@@Flatten[Lookup[rows,#]&/@{"prefixTableOK","weightedFullPrefixOK","allSuffixGridUpperOK","crossingOK","strictSeparationOK"}],Exit[2]];
Export[out60, clean[rows], "RawJSON"];
Print[rows];
