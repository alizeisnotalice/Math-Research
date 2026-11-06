(* r032: continuous, split-branch audit and two low-dimensional multiplier probes.
   Floating point evidence only.  Run: wolframscript -file receiver_continuous32.wl *)

base = DirectoryName[$InputFileName];
seedPath = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]],
   FileNameJoin[{base, "variable29_archived_theta0.json"}]];
outPath = If[Length[$ScriptCommandLine] >= 3, $ScriptCommandLine[[3]],
   FileNameJoin[{base, "receiver_continuous32.json"}]];
seed = Import[seedPath, "RawJSON"];
n = seed["n"]; b = n/2; cap = seed["height_cap"]; c = seed["c"]; tau = seed["tau_ref"];
w0 = N[seed["weights"]]; vv = N[seed["v_nonnegative"]]; dd = N[seed["receiver_d"]];
ng = Length[w0]; coeff = Join[{1.}, ConstantArray[2., ng - 1]];
u = Join[{{0., 0.}}, Transpose[{Rest[dd] + Rest[vv], Rest[dd] - Rest[vv]}],
   Transpose[{Rest[dd] - Rest[vv], Rest[dd] + Rest[vv]}]];
p = 1./(1. + Exp[c]); mz = 1. + c p; vz = 1. + c^2 p (1. - p);
lo = Min[u]; hi = Max[u] + b mz + 16. Sqrt[b vz];
baseScale = 1.0486657924183962; target = 1.00001;
ks = Range[0, b]; bin = N[Binomial[b, #] & /@ ks]; polynorm = 1./(Gamma[b] (1. + Exp[-c])^b);
gbw = bin Exp[-c ks]/(1. + Exp[-c])^b;

points[j_Integer] := If[j == 1, {{0., 0.}}, {{u[[j, 1]], u[[j, 2]]}, {u[[j + ng - 1, 1]], u[[j + ng - 1, 2]]}}];
poly[t_?NumericQ] := If[t <= 0., 0., polynorm Total[bin Clip[t - c ks, {0., Infinity}]^(b - 1)]];
Get[FileNameJoin[{base, "r032_compiled_kernel.wl"}]];
pc = buildPolyCV[b, c];
polyVec[tt_List] := pc[Developer`ToPackedArray[N[tt]]];
gbcdf[z_?NumericQ] := Total[gbw GammaRegularized[b, 0., Clip[z - c ks, {0., Infinity}]]];
phi[h_?NumericQ] := h (1. + Max[Log[Max[h, 1.]], 0.])^(n - 1);
phiPrime = (1. + Log[cap])^(n - 2) (n + Log[cap]);

(* One receiver-density vector and likelihood vector.  The latter is the
   likelihood relative to exp(-t1-t2) after the tau_ref scaling. *)
kernels[t1_?NumericQ, t2_?NumericQ] := Module[{z, rr, pp},
  z = polyVec[t1 - u[[All, 1]]] polyVec[t2 - u[[All, 2]]] Exp[Total[u, {2}]];
  rr = Join[{z[[1]]}, z[[2 ;; ng]] + z[[ng + 1 ;; 2 ng - 1]]] Exp[-tau];
  pp = Join[{z[[1]]}, z[[2 ;; ng]] + z[[ng + 1 ;; 2 ng - 1]]] Exp[-t1 - t2]; {rr, pp}];
like[t1_?NumericQ, t2_?NumericQ, ww_] := Dot[First[kernels[t1, t2]], ww];
extendWeights[ww_] := Join[{ww[[1]]}, Rest[ww], Rest[ww]];
sourceAmp[t1_?NumericQ, ww_] := extendWeights[ww] polyVec[t1 - u[[All, 1]]] Exp[Total[u, {2}] - tau];
likeAmp[t2_?NumericQ, amp_] := Dot[amp, polyVec[t2 - u[[All, 2]]]];
probVec[t1_?NumericQ, t2_?NumericQ] := Module[{z},
  z = polyVec[t1 - u[[All, 1]]] polyVec[t2 - u[[All, 2]]] Exp[Total[u, {2}] - t1 - t2];
  Join[{z[[1]]}, z[[2 ;; ng]] + z[[ng + 1 ;; 2 ng - 1]]]];
profile[s_?NumericQ] := Module[{z, r, d}, If[s < 1., Return[0.]]; If[s <= n, Return[1.]];
  d = Log[s]; z = Max[Exp[(d - Log[n])/(n - 1)] - 1., 0.];
  Do[r = (n - 2) Log[1. + z] + Log[n + z] - d; z = Max[0., z - r/((n - 2)/(1. + z) + 1./(n + z))], {10}]; Min[Exp[z], cap]];

quadCache = <||>;
stdQuad[ord_Integer] := If[KeyExistsQ[quadCache, ord], quadCache[ord], quadCache[ord] = Module[{j, ev, vec, ix, x, wt},
  j = N[DiagonalMatrix[Range[1, ord - 1]/Sqrt[4 Range[1, ord - 1]^2 - 1], 1] +
      DiagonalMatrix[Range[1, ord - 1]/Sqrt[4 Range[1, ord - 1]^2 - 1], -1]];
  {ev, vec} = Eigensystem[j]; ix = Ordering[ev]; x = ev[[ix]]; wt = 2. vec[[ix, 1]]^2;
  {x, wt}]];
quad[ord_Integer, aa_?NumericQ, bb_?NumericQ] := Module[{q = stdQuad[ord]}, {((bb - aa)/2.) q[[1]] + (aa + bb)/2., ((bb - aa)/2.) q[[2]]}];
outerQuad[nOuter_Integer] := Module[{pan = Ceiling[nOuter/16.], edge, x, wt, outx = {}, outw = {}},
  edge = Subdivide[lo, hi, pan]; Do[{x, wt} = quad[16, edge[[i]], edge[[i + 1]]]; outx = Join[outx, x]; outw = Join[outw, wt], {i, pan}]; {outx, outw}];

boundary[lev_?NumericQ, amp_] := Module[{l = lo, r = hi, mid, fl, fr},
  fl = likeAmp[l, amp] - lev; fr = likeAmp[r, amp] - lev;
  If[fl >= 0., Return[l]]; If[fr <= 0., Return[r]];
  Do[mid = (l + r)/2.; If[likeAmp[mid, amp] < lev, l = mid, r = mid], {45}]; (l + r)/2.];

constantPart[t1_?NumericQ, left_?NumericQ, right_?NumericQ, h_?NumericQ] := Module[{cost = 0., av = ConstantArray[0., ng], pt, z, fac},
  If[right <= left, Return[{cost, av}]]; cost = Exp[tau - t1] (Exp[-left] - Exp[-right]) phi[h];
  Do[pt = points[j]; Do[fac = poly[t1 - z[[1]]] Exp[z[[1]] - t1];
      If[fac > 0., av[[j]] += fac (gbcdf[right - z[[2]]] - gbcdf[left - z[[2]]])], {z, pt}], {j, ng}]; {cost, av}];
softPart[t1_?NumericQ, left_?NumericQ, right_?NumericQ, ord_Integer, amp_] := Module[{x, wt, cost = 0., av = ConstantArray[0., ng], h},
  If[right <= left, Return[{cost, av}]]; {x, wt} = quad[ord, left, right];
  Do[h = profile[likeAmp[x[[i]], amp]];
    cost += wt[[i]] Exp[tau - t1 - x[[i]]] phi[h]; av += wt[[i]] h probVec[t1, x[[i]]], {i, Length[x]}]; {cost, av}];

evaluate[ww0_, nOuter_Integer, inner_Integer] := Module[{ww = baseScale ww0, tq, tx, tw, cost = 0., av = ConstantArray[0., ng], roots, q1, q2, q3, a, amp, z, start = AbsoluteTime[]},
  {tx, tw} = outerQuad[nOuter];
  Do[amp = sourceAmp[tx[[i]], ww]; roots = boundary[#, amp] & /@ {1., N[n], phiPrime}; {q1, a} = constantPart[tx[[i]], roots[[1]], roots[[2]], 1.]; cost += tw[[i]] q1; av += tw[[i]] a;
    {q1, a} = constantPart[tx[[i]], roots[[3]], hi, cap]; cost += tw[[i]] q1; av += tw[[i]] cap a;
    {q1, a} = softPart[tx[[i]], roots[[2]], roots[[3]], inner, amp]; cost += tw[[i]] q1; av += tw[[i]] a, {i, Length[tx]}];
  <|"outer_nodes" -> Length[tx], "soft_inner_order" -> inner, "scaled_cost" -> cost,
   "actual_logJ" -> Log[cost] - tau, "averages_by_group" -> (av/coeff), "min_average" -> Min[av/coeff],
   "seconds" -> (AbsoluteTime[] - start)|>];

(* Receiver-dependent c_j check.  The source likelihood and H remain fixed
   at c; only the probability law used for each receiver average uses cj.
   Thus an observed cost reduction after rescaling is not a fixed-c dual move. *)
polyC[t_?NumericQ, cj_?NumericQ] := If[t <= 0., 0., Total[bin Clip[t - cj ks, {0., Infinity}]^(b - 1)]/(Gamma[b] (1. + Exp[-cj])^b)];
gbcdfC[z_?NumericQ, cj_?NumericQ] := Module[{ww = bin Exp[-cj ks]/(1. + Exp[-cj])^b}, Total[ww GammaRegularized[b, 0., Clip[z - cj ks, {0., Infinity}]]]];
probVecC[t1_?NumericQ, t2_?NumericQ, cc_List] := Module[{ans = ConstantArray[0., ng], pt, z, q1, q2},
 Do[pt = points[j]; Do[q1 = polyC[t1 - z[[1]], cc[[j]]]; q2 = polyC[t2 - z[[2]], cc[[j]]];
   ans[[j]] += q1 q2 Exp[z[[1]] + z[[2]] - t1 - t2], {z, pt}], {j, ng}]; ans];
constantPartC[t1_?NumericQ, left_?NumericQ, right_?NumericQ, h_?NumericQ, cc_List] := Module[{cost = 0., av = ConstantArray[0., ng], pt, z, fac},
 If[right <= left, Return[{cost, av}]]; cost = Exp[tau - t1] (Exp[-left] - Exp[-right]) phi[h];
 Do[pt = points[j]; Do[fac = polyC[t1 - z[[1]], cc[[j]]] Exp[z[[1]] - t1];
   If[fac > 0., av[[j]] += fac (gbcdfC[right - z[[2]], cc[[j]]] - gbcdfC[left - z[[2]], cc[[j]]])], {z, pt}], {j, ng}]; {cost, av}];
softPartC[t1_?NumericQ, left_?NumericQ, right_?NumericQ, ord_Integer, amp_, cc_List] := Module[{x, wt, cost = 0., av = ConstantArray[0., ng], h},
 If[right <= left, Return[{cost, av}]]; {x, wt} = quad[ord, left, right];
 Do[h = profile[likeAmp[x[[i]], amp]]; cost += wt[[i]] Exp[tau - t1 - x[[i]]] phi[h];
   av += wt[[i]] h probVecC[t1, x[[i]], cc], {i, Length[x]}]; {cost, av}];
evaluateReceiverC[cc_List, srcScale_?NumericQ, nOuter_Integer, inner_Integer] := Module[{ww = srcScale w0, tx, tw, cost = 0., av = ConstantArray[0., ng], roots, q1, a, amp, start = AbsoluteTime[]},
 {tx, tw} = outerQuad[nOuter];
 Do[amp = sourceAmp[tx[[i]], ww]; roots = boundary[#, amp] & /@ {1., N[n], phiPrime};
   {q1, a} = constantPartC[tx[[i]], roots[[1]], roots[[2]], 1., cc]; cost += tw[[i]] q1; av += tw[[i]] a;
   {q1, a} = constantPartC[tx[[i]], roots[[3]], hi, cap, cc]; cost += tw[[i]] q1; av += tw[[i]] cap a;
   {q1, a} = softPartC[tx[[i]], roots[[2]], roots[[3]], inner, amp, cc]; cost += tw[[i]] q1; av += tw[[i]] a, {i, Length[tx]}];
 <|"source_scale" -> srcScale, "receiver_c" -> cc, "outer_nodes" -> Length[tx], "soft_inner_order" -> inner,
  "scaled_cost" -> cost, "actual_logJ" -> Log[cost] - tau, "averages_by_group" -> av/coeff,
  "min_average" -> Min[av/coeff], "seconds" -> AbsoluteTime[] - start|>];

mode = If[TrueQ[$r032LibraryMode], "library", If[Length[$ScriptCommandLine] >= 4, $ScriptCommandLine[[4]], "full"]];
If[mode === "quick", quick = evaluate[w0, 16, 8]; Export[outPath, <|"status" -> "quick timing only", "quick" -> quick|>, "RawJSON"]; Print[quick]; Exit[]];
If[mode === "allccompare",
  localC = ConstantArray[c, ng]; localC[[{23, 24, 25}]] += .05;
  allC = ConstantArray[c + .05, ng];
  local256 = evaluateReceiverC[localC, baseScale, 256, 32];
  all256 = evaluateReceiverC[allC, baseScale, 256, 32];
  compare = <|"status" -> "same-source-H receiver-c comparison; floating numerical evidence only",
    "source_c" -> c, "source_scale" -> baseScale, "outer_nodes" -> 256, "soft_inner_order" -> 32,
    "local_groups" -> {23, 24, 25}, "local_group_receiver_points" -> "for j>1, both mirror points u[j] and u[j+32] use the same c_j",
    "local_c" -> localC, "all_c" -> allC, "local_only" -> local256, "all_receivers" -> all256|>;
  Export[outPath, compare, "RawJSON"]; Print[<|"local_min" -> local256["min_average"], "all_min" -> all256["min_average"]|>]; Exit[]];

If[mode =!= "library",
(* Compute the continuous q-Avg residual at the safe r031 scale.  Directions
   preserve total multiplier mass: move mass from the two most surplus active
   groups to the most deficient group.  This is precisely a projected q-Avg
   direction, not a new full-dimensional optimization. *)
base128 = evaluate[w0, 128, 32]; base256x32 = evaluate[w0, 256, 32]; base256x64 = evaluate[w0, 256, 64];
av = base256x64["averages_by_group"]; deficit = First@Ordering[av, 1]; active = Flatten@Position[w0, _?(# > 10^-12 &)];
surplus = Take[Reverse@SortBy[active, av[[#]] &], UpTo[2]];
dirs = Table[ReplacePart[ConstantArray[0., ng], {deficit -> 1., surplus[[k]] -> -1.}], {k, Length[surplus]}];
epss = {0.002, 0.004};
probes = Flatten[Table[With[{wc = w0 + e dirs[[k]]},
   If[Min[wc] < -10^-12, <|"direction" -> k, "epsilon" -> e, "invalid" -> True|>,
    Join[<|"direction" -> k, "epsilon" -> e, "from_group" -> surplus[[k]], "to_group" -> deficit|>, evaluate[wc, 128, 32]]]], {k, Length[dirs]}, {e, epss}], 1];
eligible = Select[probes, KeyExistsQ[#, "min_average"] && #["min_average"] >= target &];
chosen = If[Length[eligible] == 0, Missing["No128Candidate"], First@SortBy[eligible, #["scaled_cost"] &]];
final = If[AssociationQ[chosen], With[{wc = w0 + chosen["epsilon"] dirs[[chosen["direction"]]]},
   <|"choice" -> chosen, "weights" -> wc, "N256_inner32" -> evaluate[wc, 256, 32], "N256_inner64" -> evaluate[wc, 256, 64]|>], <|"choice" -> chosen|>];
bestImproves = AssociationQ[chosen] && KeyExistsQ[final, "N256_inner64"] && final["N256_inner64"]["min_average"] >= target && final["N256_inner64"]["scaled_cost"] < base256x64["scaled_cost"];
(* The three binding groups are 23--25 at the corrected continuous baseline.
   This is a direct numerical perturbation only: c changes both the Bernoulli
   success probability and its shift, so there is no parameterwise FOSD claim. *)
receiverGroups = {23, 24, 25}; receiverC = ConstantArray[c, ng]; receiverC[[receiverGroups]] += .05;
receiverAtBase128 = evaluateReceiverC[receiverC, baseScale, 128, 32];
scaleDrop = Min[.02, Max[0., .8 (receiverAtBase128["min_average"] - target)/.04]];
receiverTrialScale = baseScale (1. - scaleDrop);
receiverTrial128 = evaluateReceiverC[receiverC, receiverTrialScale, 128, 32];
receiverFinal = If[receiverTrial128["min_average"] >= target,
 <|"N256_inner32" -> evaluateReceiverC[receiverC, receiverTrialScale, 256, 32],
   "N256_inner64" -> evaluateReceiverC[receiverC, receiverTrialScale, 256, 64]|>, <|"not_validated_at_256" -> True|>];
result = <|"status" -> "r032 split-boundary continuous quadrature; floating evidence only, not interval certified", "seed" -> seedPath,
 "base_uniform_scale" -> baseScale, "finite_domain" -> {lo, hi}, "target_min_average" -> target,
 "baseline_N128_inner32" -> base128, "baseline_N256_inner32" -> base256x32, "baseline_N256_inner64" -> base256x64,
 "q_minus_average" -> (1. - av), "deficit_group" -> deficit, "surplus_active_groups" -> surplus,
 "directions" -> dirs, "probes_N128_inner32" -> probes, "final" -> final,
 "receiver_dependent_c_mechanism" -> <|"source_c" -> c, "raised_receiver_groups" -> receiverGroups, "receiver_c_increment" -> .05,
   "interpretation" -> "direct floating-point recomputation only; Gamma_b+c Binomial_b is not parameterwise first-order stochastically ordered in c",
   "at_base_scale_N128_inner32" -> receiverAtBase128, "trial_source_scale" -> receiverTrialScale,
   "trial_N128_inner32" -> receiverTrial128, "final" -> receiverFinal|>,
 "strict_improvement_over_r031_at_N256_inner64" -> bestImproves|>;
Export[outPath, result, "RawJSON"];
Print["r032 finished; baseline scaled cost=", base256x64["scaled_cost"], " final=", If[KeyExistsQ[final, "N256_inner64"], final["N256_inner64"]["scaled_cost"], "none"]];
];
