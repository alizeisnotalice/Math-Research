(* Exact checks for the r059 geometric-mesh construction.
   All parameters are rational in the executed rows. *)

ClearAll[phi, meshData, directCost, formulaCost, prefixAverage,
  directFibrePrefix, suffixValue, suffixMax, fullValue, runSmall, runLarge,
  jsonClean];

phi[n_Integer][t_] := Piecewise[{{t, 0 <= t <= 1}},
  t (1 + Log[t])^(n - 1)];

meshData[n_Integer, m_Integer, rr_Rational] := Module[
  {q, radii, hs, ds, low, high},
  q = 1 + n m rr;
  radii = Join[{0}, Table[q^(j - 1), {j, 1, m}]];
  hs = Differences[radii];
  ds = m rr Rest[radii];
  low = Table[1 - ds[[i]]/hs[[i]], {i, m}];
  high = Table[1 + ds[[i]]/hs[[i + 1]], {i, m - 1}];
  <|"N" -> n, "M" -> m, "R" -> rr, "Q" -> q, "r" -> radii,
    "h" -> hs, "D" -> ds, "low" -> low, "high" -> high|>
];

(* Integrate the actual cell values, rather than using a pre-summed formula. *)
directCost[d_Association] := Module[{n = d["N"], m = d["M"], h = d["h"],
    low = d["low"], high = d["high"], rm, total},
  rm = Last[d["r"]];
  total = Sum[h[[j]] Sum[
      phi[n][Which[i == j, low[[i]], i < m && i + 1 == j, high[[i]],
        True, 1]], {i, m}], {j, m}];
  total/(m rm)
];

formulaCost[d_Association] := Module[{n = d["N"], m = d["M"], rr = d["R"],
    rm = Last[d["r"]]},
  1 - rr + (1 - 1/rm)/m (phi[n][1 + 1/n] - (1 + 1/n))
];

(* Average of fibre i over the anchored inner prefix [0,r_j]. *)
prefixAverage[d_Association, i_Integer, j_Integer] := Module[
  {m = d["M"], rr = d["R"]}, If[i == j, 1 - m rr, 1]];

directFibrePrefix[d_Association, i_Integer, j_Integer] := Module[
  {m = d["M"], h = d["h"], low = d["low"], high = d["high"], r},
  r = d["r"][[j + 1]];
  Sum[h[[k]] Which[i == k, low[[i]], i < m && i + 1 == k,
      high[[i]], True, 1], {k, j}]/r
];

(* At outer receiver M+b, a suffix [a,M+b] which drops the deficient
   j-th fibre has its best value at a=j.  The all-fibre choice is a=0. *)
suffixValue[m_Integer, j_Integer, b_Rational] := (m - j)/(m - j + b);
fullValue[m_Integer, rr_Rational, b_Rational] := (m - m rr)/(m + b);
suffixMax[m_Integer, b_Rational] := (m - 1)/(m - 1 + b);

runSmall[] := Flatten[Table[Module[{d, dc, fc, m = mm, rr = 1/(4 mm),
      b = 1/3, prefOK, properOK, outerOK},
    d = meshData[nn, mm, rr]; dc = directCost[d]; fc = formulaCost[d];
    prefOK = And @@ Flatten[Table[
       FullSimplify[directFibrePrefix[d, i, j] == prefixAverage[d, i, j]],
       {i, m}, {j, m}]];
    properOK = And @@ Table[
       And @@ Table[directFibrePrefix[d, i, j] == 1,
         {i, DeleteCases[Range[m], j]}],
       {j, m}];
    outerOK = And @@ Table[
       suffixValue[m, j, b] <= suffixMax[m, b], {j, m}] &&
       fullValue[m, rr, b] <= Max[suffixMax[m, b], fullValue[m, rr, b]];
    <|"N" -> nn, "M" -> mm, "R" -> rr, "Q" -> d["Q"],
      "directEqualsFormula" -> (FullSimplify[dc == fc]),
      "allGridPrefixAverage" -> (1 - rr), "prefixTableOK" -> prefOK,
      "properSubsetWitnessOK" -> properOK, "outerEndpointChecksOK" -> outerOK,
      "cost" -> dc|>], {nn, {2, 4, 8, 16}}, {mm, {2, 3, 5}}], 1];

runLarge[] := Table[Module[{m, ell, rr, d, eta, b, lambda, lower, cellLower},
    m = Floor[nn^(1/4)]; ell = Ceiling[N[Log[nn]^2, 80]];
    rr = ell/nn; d = meshData[nn, m, rr];
    (* eta is in the outerly rescaled coordinate: receiver y=1+eta.
       b=M eta is its physical separation when the outer support is [0,M]. *)
    eta = (m - 1) rr/(1 - m rr); b = m eta;
    lambda = 1 + eta;
    (* The recovery cell of fibre M-1 on the largest inner cell alone is
       a lower bound for the scaled actual modular. *)
    lower = (1 - 1/d["Q"])/m phi[nn][lambda];
    cellLower = (1 - 1/d["Q"])/m phi[nn][lambda (1 + 1/nn)];
    <|"N" -> nn, "M" -> m, "L" -> ell, "R" -> rr, "Q" -> d["Q"],
      "MR" -> m rr, "nonnegative" -> (m rr + 1/nn < 1),
      "rawCost" -> N[formulaCost[d], 30], "b" -> b,
      "etaStar" -> eta, "adaptiveOuter" -> N[1/(1 + eta), 30],
      "escapeSuffix" -> N[suffixMax[m, b], 30],
      "fullOuter" -> N[fullValue[m, rr, b], 30],
      "branchesEqual" -> (suffixMax[m, b] == fullValue[m, rr, b]),
      "thresholdNormalizer" -> N[lambda, 30],
      "oneRecoveryCellLowerBound" -> N[lower, 30],
      "actualRecoveryCellContribution" -> N[cellLower, 30],
      "logCellLowerBound" -> N[Log[lower], 30]|>],
   {nn, {64, 256, 1024, 4096, 16384, 65536}}];

small = runSmall[];
large = runLarge[];
precisionCheck = Table[Module[{n=row["N"],m=row["M"],r=row["R"],d,z,exprs,delta},
 d=meshData[n,m,r];z=(1-r)/(1-m r);
 exprs={formulaCost[d],Log[(1-1/d["Q"])/m phi[n][z]]};
 delta=Abs[N[exprs,30]-N[exprs,60]];
 <|"N"->n,"precisions"->{30,60},"absoluteDifferences_WolframInputForm"->(ToString[#,InputForm]& /@ delta),
 "pass"->TrueQ[Max[delta]<10^-25]|>],{row,large}];
jsonClean[x_Association] := Association @@
  (Rule[#[[1]], jsonClean[#[[2]]]] & /@ Normal[x]);
jsonClean[x_List] := jsonClean /@ x;
jsonClean[x_Integer] := x;
jsonClean[x_Rational] := ToString[x, InputForm];
jsonClean[x_?BooleanQ] := x;
jsonClean[x_DirectedInfinity] := ToString[x, InputForm];
jsonClean[x_?NumericQ] := N[x, 30];
jsonClean[x_] := x;
out59=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/mesh59.json"];
If[!TrueQ[And@@Flatten[Lookup[small,#]&/@{"directEqualsFormula","prefixTableOK","properSubsetWitnessOK","outerEndpointChecksOK"}]] || !TrueQ[And@@Lookup[large,"branchesEqual"]] || !TrueQ[And@@Lookup[precisionCheck,"pass"]],Exit[2]];
Export[out59,
  jsonClean[<|"smallExact" -> small, "largeRationalDiagnostics" -> large,"precisionDoubling"->precisionCheck|>], "RawJSON"];
Print["small exact rows: ", Length[small], "; all direct formula checks: ",
  And @@ Lookup[small, "directEqualsFormula"]];
Print["large rational diagnostics:"];
Print[large];
