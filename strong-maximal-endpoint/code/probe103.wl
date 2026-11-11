(* High-precision finite-parameter diagnostics; no asymptotic theorem is claimed. *)
ClearAll[numString, stableGMSResidual, computeCase, guardedCase];
$HistoryLength = 0;
wp = 300;
taskMemoryBytes = 160*2^20;
taskSeconds = 90;
reserveBytes = 2^30;

numString[x_, prec_: 32] := ToString[FortranForm[N[x, prec]], OutputForm];

(* Divide the defining equation by no small cancellation-prone quantity:
   after b=P Log[c], the original residual is exactly this stable form. *)
stableGMSResidual[b_?NumericQ, p_?NumericQ] := Module[{q = 1/p},
  Exp[b]*(-1 - (Exp[-q*b] - 1)/q) - 1
];

computeCase[{l_Integer, kName_String, kSpec_}] := Module[
  {k, p, w, delta, upperRate, lowerRate, logTailUpper, logMoment,
   stirlingApprox, stirlingResidual, root, b, c, stableResidual,
   directResidual, prefixCost, lowTailActive, memoryFloor},
  k = If[kName === "exp_sqrt_L", N[Exp[Sqrt[l]], wp], N[kSpec, wp]];
  p = Max[N[l^2, wp], k];
  w = 8*(Sqrt[p*l] + l);
  delta = w/p;

  upperRate = p*(delta - Log[1 + delta]);
  lowTailActive = delta < 1;
  If[lowTailActive,
    lowerRate = p*(-delta - Log[1 - delta]);
    logTailUpper = -Min[upperRate, lowerRate] +
      Log[1 + Exp[-Abs[upperRate - lowerRate]]],
    lowerRate = Missing["LowerEndpointNonpositive"];
    logTailUpper = -upperRate
  ];

  logMoment = p + LogGamma[p + 1] - p*Log[p];
  stirlingApprox = Log[2*Pi*p]/2 + 1/(12*p) - 1/(360*p^3) +
    1/(1260*p^5) - 1/(1680*p^7);
  stirlingResidual = logMoment - stirlingApprox;

  root = FindRoot[stableGMSResidual[b, p] == 0,
    {b, N[13/10, wp]}, WorkingPrecision -> wp,
    AccuracyGoal -> 220, PrecisionGoal -> 220, MaxIterations -> 80];
  b = b /. root;
  c = Exp[b/p];
  stableResidual = stableGMSResidual[b, p];
  directResidual = (p - 1)*c^p - p*c^(p - 1) - 1;
  prefixCost = (k/p + Log[p])/l;
  memoryFloor = MemoryAvailable[];

  <|
    "status" -> "pass",
    "L" -> l,
    "k_family" -> kName,
    "k_definition" -> Switch[kName, "L", "L", "L^2", "L^2", "L^4", "L^4", "exp_sqrt_L", "exp(sqrt(L))"],
    "k_value_numeric" -> numString[k, 40],
    "P_definition" -> "max(L^2,k)",
    "P_value_numeric" -> numString[p, 40],
    "W_definition" -> "8*(sqrt(P*L)+L)",
    "W_value_numeric" -> numString[w, 40],
    "truncated_moment_formula" -> "exp(P)*Gamma(P+1)/P^P",
    "log_truncated_moment_high_precision_numeric" -> numString[logMoment, 60],
    "stirling_log_approximation_order" -> "log(2*pi*P)/2 + 1/(12P) - 1/(360P^3) + 1/(1260P^5) - 1/(1680P^7)",
    "stirling_residual_high_precision_numeric" -> numString[stirlingResidual, 60],
    "center_window" -> "[max(0,P-W),P+W] for X~Gamma(shape=P, rate=1)",
    "center_window_lower_tail_active" -> lowTailActive,
    "center_window_tail_log_bound_type" -> "strict Chernoff upper; not an evaluated exact tail",
    "center_window_tail_log_bound_formula" -> If[lowTailActive,
      "log(exp(-P*(delta-log(1+delta))) + exp(-P*(-delta-log(1-delta)))) where delta=W/P",
      "-P*(delta-log(1+delta)) where delta=W/P; lower tail is zero"],
    "center_window_tail_log_upper_numeric" -> numString[logTailUpper, 60],
    "center_window_upper_rate_numeric" -> numString[upperRate, 40],
    "center_window_lower_rate_numeric" -> If[lowTailActive, numString[lowerRate, 40], "not_applicable"],
    "gms_root_equation" -> "(P-1)c^P - P*c^(P-1) - 1 = 0, c>1",
    "gms_b_definition" -> "b=P*log(c)",
    "gms_stable_residual_formula" -> "exp(b)*(-1-(exp(-b/P)-1)*P)-1",
    "gms_b_root_numeric" -> numString[b, 60],
    "gms_c_root_numeric" -> numString[c, 60],
    "gms_stable_equation_residual_numeric" -> numString[stableResidual, 60],
    "gms_direct_c_equation_residual_numeric" -> numString[directResidual, 60],
    "multilog_prefix_cost_formula" -> "(k/P + log(P))/L",
    "multilog_prefix_cost_high_precision_numeric" -> numString[prefixCost, 60],
    "worker_memory_available_bytes_sample" -> memoryFloor
  |>
];

guardedCase[spec_] := MemoryConstrained[
  TimeConstrained[computeCase[spec], taskSeconds,
    <|"status" -> "timeout", "L" -> First[spec], "k_family" -> spec[[2]]|>],
  taskMemoryBytes,
  <|"status" -> "memory_limit", "L" -> First[spec], "k_family" -> spec[[2]]|>
];

lValues = {16, 32, 64, 128, 256, 512, 1024, 2048, 4096};
tasks = Flatten[Table[
  {{l, "L", l}, {l, "L^2", l^2}, {l, "L^4", l^4},
   {l, "exp_sqrt_L", Exp[Sqrt[l]]}},
  {l, lValues}], 1];

If[MemoryAvailable[] <= reserveBytes,
  Print["Insufficient available memory before kernel launch; refusing to run."];
  Exit[2]
];

launch = LaunchKernels[12];
kernelIds = ParallelEvaluate[$KernelID];
If[Length[kernelIds] != 12,
  CloseKernels[];
  Print["Expected 12 parallel Wolfram kernels."];
  Exit[3]
];
DistributeDefinitions[numString, stableGMSResidual, computeCase, guardedCase,
  wp, taskMemoryBytes, taskSeconds];
rows = ParallelMap[guardedCase, tasks, Method -> "CoarsestGrained"];
allPass = Length[rows] == 36 && AllTrue[rows, AssociationQ[#] && # ["status"] === "pass" &];
CloseKernels[];
kernelsClosed = Length[Kernels[]] == 0;
result = <|
  "schema" -> "probe103-finite-high-precision-v1",
  "pass_all" -> (allPass && kernelsClosed),
  "case_count" -> Length[rows],
  "expected_case_count" -> 36,
  "configured_workers" -> 12,
  "worker_kernel_ids" -> kernelIds,
  "kernels_closed" -> kernelsClosed,
  "working_precision_digits" -> wp,
  "per_task_memory_limit_MiB" -> taskMemoryBytes/2^20,
  "per_task_time_limit_seconds" -> taskSeconds,
  "system_available_memory_reserve_GiB" -> 1,
  "moment_evaluation_kind" -> "high-precision numerical evaluation of the exact symbolic Gamma formula",
  "tail_bound_kind" -> "strict analytic Chernoff bound; reported logarithm is its high-precision numerical evaluation",
  "gms_root_source" -> "r100/strip100.md, equation (P-1)c^P-Pc^(P-1)-1=0",
  "scope" -> "Finite parameter diagnostics only; no general maximal-operator bound is claimed.",
  "rows" -> rows
|>;
outputPath = FileNameJoin[{DirectoryName[$InputFileName], "probe103-results.json"}];
Export[outputPath, result, "JSON"];
Print[ExportString[<|"pass_all" -> result["pass_all"],
  "case_count" -> result["case_count"], "configured_workers" -> 12,
  "kernels_closed" -> kernelsClosed, "first_rows" -> Take[rows, UpTo[4]]|>,
  "JSON", "Compact" -> True]];
If[!TrueQ[result["pass_all"]], Exit[4]];
Exit[0];
