args95 = Rest[$ScriptCommandLine];
out95 = args95[[1]];
workers95 = If[Length[args95] >= 2, ToExpression[args95[[2]]], 1];
If[!DirectoryQ[out95], CreateDirectory[out95, CreateIntermediateDirectories -> True]];
If[workers95 < 1 || workers95 > 2, Print["workers must be 1 or 2"]; Exit[2]];

ClearAll[str95, overlapRow95, generalRow95, rowText95, summaryRow95];
str95[x_] := ToString[x, InputForm];
rowText95[row_List] := StringRiffle[Join[{If[StringQ[row[[1]]], row[[1]], str95[row[[1]]]]},
  str95 /@ Rest[row]], "|"];
summaryRow95[row_List] := <|
  "m" -> row[[1]], "k" -> row[[2]], "p" -> str95[row[[3]]],
  "raw" -> str95[row[[4]]], "surviving" -> str95[row[[5]]],
  "loss" -> str95[row[[6]]], "pair_bound" -> str95[row[[7]]]|>;

overlapRow95[bits_List] := Module[
  {p = 1/31, q, probability, raw, colors, surviving, loss, pairCount},
  q = 1 - p;
  probability = Times @@ (If[# == 1, p, q] & /@ bits);
  raw = Total[bits];
  colors = Table[Total[bits[[{2 i - 1, 2 i}]]] > 0, {i, 1, 3}];
  surviving = Total[Table[If[colors[[i]] && Not[Or @@ Delete[colors, i]], 1, 0], {i, 1, 3}]];
  loss = raw - surviving;
  pairCount = raw (raw - 1);
  {StringJoin[ToString /@ bits], probability, raw, surviving, loss, pairCount}
];

overlapBits95 = Tuples[{0, 1}, 6];
overlapRows95 = overlapRow95 /@ overlapBits95;
overlapRaw95 = Total[(#[[2]] #[[3]]) & /@ overlapRows95];
overlapSurv95 = Total[(#[[2]] #[[4]]) & /@ overlapRows95];
overlapLoss95 = Total[(#[[2]] #[[5]]) & /@ overlapRows95];
overlapPair95 = Total[(#[[2]] #[[6]]) & /@ overlapRows95];
overlapHash95 = Hash[StringRiffle[rowText95 /@ overlapRows95, "\n"], "SHA256", "HexString"];
overlapPass95 = Length[overlapRows95] == 64 &&
  overlapRaw95 == 6/31 && overlapSurv95 == 148230000/887503681 &&
  overlapLoss95 == 23544906/887503681 && overlapPair95 == 30/961 &&
  overlapLoss95 <= overlapPair95 && AllTrue[overlapRows95, #[[5]] <= #[[6]] &];
overlapData95 = <|
  "patterns_enumerated" -> Length[overlapRows95],
  "pattern_rows_sha256" -> overlapHash95,
  "raw_mass_expectation" -> str95[overlapRaw95],
  "surviving_mass_expectation" -> str95[overlapSurv95],
  "lost_mass_expectation" -> str95[overlapLoss95],
  "expected_ordered_pair_bound" -> str95[overlapPair95],
  "exact_loss_le_pair_bound" -> TrueQ[overlapLoss95 <= overlapPair95],
  "pointwise_loss_le_l_times_l_minus_1" -> TrueQ[AllTrue[overlapRows95, #[[5]] <= #[[6]] &]]|>;

c95 = 1/1024; eta95 = 1/128; delta95 = 1/1024;
lambda95 = 5/4 c95 delta95;
good95 = 5/8 (1 - 5/4 c95);
packLower95 = (good95 - 3/10)/(7/10) - 4 lambda95;
guardLower95 = 1 - 38 eta95;
longEdge95 = 2 eta95 - eta95/2;
constantPass95 = good95 == 20455/32768 && packLower95 == 485691/1048576 &&
  guardLower95 == 45/64 && longEdge95 == 3/256 && longEdge95 > eta95;
constantData95 = <|
  "c" -> str95[c95], "eta" -> str95[eta95], "delta" -> str95[delta95],
  "lambda_upper" -> str95[lambda95],
  "minimum_expected_success_fraction" -> str95[good95],
  "packing_simultaneous_success_probability_lower" -> str95[packLower95],
  "guard_simultaneous_success_probability_lower" -> str95[guardLower95],
  "long_edge_deleted_fraction_lower" -> str95[longEdge95],
  "strictly_exceeds_eta" -> TrueQ[longEdge95 > eta95],
  "derived_with_exact_fraction_arithmetic" -> True|>;

generalRow95[{m_, k_, p_}] := Module[{raw, surv, loss, bound, ratio},
  raw = m k p;
  surv = m (1 - (1 - p)^k) (1 - p)^(k (m - 1));
  loss = raw - surv;
  bound = m k (m k - 1) p^2;
  ratio = If[bound == 0, 0, loss/bound];
  {m, k, p, raw, surv, loss, bound, ratio}
];
generalRows95 = generalRow95 /@ Tuples[{Range[16], Range[8], {1/31, 1/127}}];
generalCanonicalRows95 = Take[#, 7] & /@ generalRows95;
generalHash95 = Hash[StringRiffle[rowText95 /@ generalCanonicalRows95, "\n"], "SHA256", "HexString"];
minSurvRow95 = First[MinimalBy[generalRows95, #[[5]] &]];
maxSurvRow95 = First[MaximalBy[generalRows95, #[[5]] &]];
minPairMarginRow95 = First[MinimalBy[generalRows95, #[[7]] - #[[6]] &]];
maxLossRatioRow95 = First[MaximalBy[generalRows95, #[[8]] &]];
generalPass95 = Length[generalRows95] == 256 && AllTrue[generalRows95, 0 <= #[[6]] <= #[[7]] &];
generalData95 = <|
  "case_count" -> Length[generalRows95],
  "canonical_case_rows_sha256" -> generalHash95,
  "loss_le_pair_bound_all_cases" -> TrueQ[AllTrue[generalRows95, 0 <= #[[6]] <= #[[7]] &]],
  "minimum_surviving_mass" -> summaryRow95[Take[minSurvRow95, 7]],
  "maximum_surviving_mass" -> summaryRow95[Take[maxSurvRow95, 7]],
  "minimum_pair_bound_minus_loss" -> <|"value" -> str95[minPairMarginRow95[[7]] - minPairMarginRow95[[6]]],
    "m" -> minPairMarginRow95[[1]], "k" -> minPairMarginRow95[[2]], "p" -> str95[minPairMarginRow95[[3]]]|>,
  "maximum_loss_over_pair_bound" -> <|"value" -> str95[maxLossRatioRow95[[8]]],
    "m" -> maxLossRatioRow95[[1]], "k" -> maxLossRatioRow95[[2]], "p" -> str95[maxLossRatioRow95[[3]]]|>|>;

theta95 = Unique["theta95"]; h95 = Unique["h95"];
gap95 = -Log[1 - theta95] - theta95;
logFacts95 = <|
  "gap_at_zero" -> TrueQ[(gap95 /. theta95 -> 0) == 0],
  "gap_derivative_identity" -> TrueQ[FullSimplify[D[gap95, theta95] == theta95/(1 - theta95),
      Assumptions -> theta95 < 1]],
  "gap_derivative_positive_on_domain" -> TrueQ[Resolve[ForAll[theta95,
      Implies[0 < theta95 < 1, theta95/(1 - theta95) > 0]], Reals]],
  "positive_H_preserves_order" -> TrueQ[Resolve[ForAll[h95,
      Implies[h95 > 0, 1/h95 > 0]], Reals]]|>;
logContract95 = <|
  "statement" -> "theta/H_N(p) <= d_N(p,theta) <= -log(1-theta)/H_N(p)",
  "domain" -> "0 < theta < 1; H_N(p) > 0",
  "exact_ordering_facts" -> logFacts95,
  "numeric_log_evaluation_used" -> False,
  "scope" -> "checks only the exact logarithmic boundary ordering; the d_N theorem is supplied by core95 (6)"|>;

If[!DirectoryQ[FileNameJoin[{out95, "wolfram"}]],
  CreateDirectory[FileNameJoin[{out95, "wolfram"}], CreateIntermediateDirectories -> True]];
pass95 = TrueQ[overlapPass95 && constantPass95 && generalPass95 && And @@ Values[logFacts95]];
verification95 = <|
  "pass" -> pass95,
  "configured_workers" -> 1,
  "actual_kernel_ids" -> {},
  "main_kernel_computation" -> True,
  "kernels_closed" -> True,
  "overlap" -> overlapData95,
  "constants" -> constantData95,
  "generalization" -> generalData95,
  "log_contract" -> logContract95,
  "scope" -> "Independent exact checks of finite Bernoulli point models, rational margins, and a symbolic logarithmic order contract; not a physical onset or high-dimensional source construction."|>;
Export[FileNameJoin[{out95, "wolfram", "verification95.json"}], verification95, "RawJSON"];
CloseKernels[];
Print[If[pass95, "PASS ", "FAIL "], "overlap patterns=", Length[overlapRows95],
  "; generalized cases=", Length[generalRows95], "; main kernel only"];
If[!pass95, Exit[3]];
