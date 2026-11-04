(* Exact scalar certificates used in the independent review. No time limits. *)
ClearAll["Global`*"];
SetAttributes[record, HoldRest];
failures = 0;
record[label_String, expr_] := Module[{res = expr},
  Print[label, " = ", ToString[res, InputForm]];
  If[res =!= True, failures++]
];
noCounterexample[vars_List, claim_] := Resolve[Exists[vars, claim], Reals] === False;

Print["Kernel: ", $Version];
eta = 1/64;
hUpper = eta/(16 9^2);
record["Local clipping and the deep-cover slack, all n >= 2",
  3/8 + hUpper < (2/5) (1 - eta/4) && 1/2 - eta/2 + hUpper < 1/2 - eta/8];
record["The thin strip fits outside any half-box meeting a private center",
  noCounterexample[{d, y1, xk, yk},
    0 < d <= eta && 1/2 < y1 < 1/2 + d/4 && -d <= xk <= d && -1/8 < yk < 1/8 &&
    (Abs[d - y1] > 1/2 - d/4 || Abs[xk - yk] > 1/2 - d/4)]];
record["Individual large witness captures all of its deep cell",
  noCounterexample[{w, delta}, w >= 0 && delta >= 0 && w + delta < 1/10 && 9/10 + w >= 1 - delta]];
record["Threshold integral",
  FullSimplify[Integrate[1/l, {l, T, T^2}, Assumptions -> T > 2] == Log[T], T > 2]];

(* Normalize s = r_i + r_j = 1; w = c delta, u = r_i/s. *)
range = 0 < delta <= w < 1/64 && 1/2 - w <= u <= 1/2 + w;
record["Neighbor radii and half-box packing radius",
  noCounterexample[{delta, w, u}, range &&
    (u < (1-u)/2 || 1-u < u/2 || 1/2 + w + (1-u)/2 > 3 u)]];
record["Critical-coordinate intersection width",
  noCounterexample[{delta, w, u, g}, range && g >= 1/2 - 3 w - 3 delta &&
    (1/2+w) u - (g-(1/2+w)(1-u)) > 4 w + 3 delta]];
record["Every noncritical separation fits both witness lower radii, the critical sign is strict",
  noCounterexample[{delta, w, u}, range &&
    ((1-delta) Min[u, 1-u] <= 1/2-2 w-3 delta || 1/2-4 w-3 delta <= 0)]];
record["Captured-neighbor normal interval for a nonattained private endpoint",
  noCounterexample[{r, d, eps, sep}, r > 0 && d < r && 0 < eps < r-d && 0 < sep <= d+eps && Abs[sep] >= r]];
record["A tail beats the slack uniformly after the radius comparison (q=2^-n, ratio=delta/b)",
  noCounterexample[{q, delta, ratio}, 0 < q <= 1/4 && 0 < delta < 1/2 && 0 < ratio <= q/2 && q (1-delta) <= ratio]];
record["Bilateral separation implies zero normal overlap",
  noCounterexample[{d, ri, rj}, ri > 0 && rj > 0 && d >= ri && d >= rj && d < (ri+rj)/2]];

record["Quantile and overlap exponents for T=v^-1/8",
  FullSimplify[(v^(-1/8))^2 v == v^(3/4) && v^(3/4)/Sqrt[v] == v^(1/4) &&
    (v^(3/4))^2/v == Sqrt[v] && 1/Log[v^(-1/8)] == 8/Log[1/v], 0 < v < 1]];
record["Periodic lower envelope separates its first-order binomial and quadratic reciprocal errors",
  Factor[(a p/(1-a p)-a)/(a v) - (p-1)/v - a p^2/(v(1-a p))] === 0];
record["First binomial coefficient for arbitrary n",
  FullSimplify[(D[(1+v)^(-n), v] /. v -> 0) == -n, n >= 2]];
record["Uniform second derivative bound on 0<=v<=1, fixed real n>=2",
  FullSimplify[D[(1+v)^(-n), {v, 2}] == n (n+1) (1+v)^(-n-2), n >= 2 && 0 <= v <= 1] &&
  FullSimplify[Log[1+v] >= 0, 0 <= v <= 1] &&
  noCounterexample[{n, logBase}, n >= 2 && logBase >= 0 && (-n-2) logBase > 0] &&
  noCounterexample[{exponent}, exponent <= 0 && (Exp[exponent] <= 0 || Exp[exponent] > 1)]];
record["Reciprocal error is at most of order a v^2 for n>=2 (a=v^n)",
  FullSimplify[Log[v] < 0, 0 < v < 1] &&
  noCounterexample[{n, logBase}, n >= 2 && logBase < 0 && (2 n - (n+2)) logBase > 0] &&
  noCounterexample[{exponent}, exponent <= 0 && Exp[exponent] > 1]];
record["Return to the original fractional coefficient",
  FullSimplify[2^n 2 == 2^(n+1), n >= 2]];
Print["Failures: ", failures];
Exit[If[failures === 0, 0, 1]];
