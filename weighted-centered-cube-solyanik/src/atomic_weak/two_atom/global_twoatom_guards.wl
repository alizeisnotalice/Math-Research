(* Exact symbolic guards for the all-tau two-atom diagnostic.
   They prove identities only; positivity and root uniqueness remain analytic. *)
$MaxExtraPrecision = 2000;
ClearAll[guard, r, t, d, n, ff, h, j, rr];
guard[name_, q_] := If[TrueQ[q], Print["PASS ", name], Print["FAIL ", name, ": ", q]; Exit[1]];

n = (r - 1)^t + (r + 1)^t;
ff = r D[n, r] - n;
guard["F_4over3_factorization",
 FullSimplify[(ff /. r -> 4/3)/3^-t == 4 t - 1 + (4 t - 7) 7^(t - 1),
  Assumptions -> 1/4 < t < 1]];

h = Log[(7 - 4 t)/(4 t - 1)] + (t - 1) Log[7];
guard["H_prime",
 FullSimplify[D[h, t] == Log[7] - 24/((7 - 4 t) (4 t - 1)),
  Assumptions -> 1/4 < t < 1]];
guard["H_at_one", FullSimplify[Limit[h, t -> 1, Direction -> "FromBelow"] == 0]];

j = (d^t + (2 + d)^t)/(2 t (1 + d));
rr = ((1 + 2 d)^t - (1 - 2 d)^t + d^t + (2 - d)^t)/(2 t (1 + d));
guard["R_minus_J_integral_numerator",
 FullSimplify[rr - j == ((1 + 2 d)^t - (1 - 2 d)^t - (2 + d)^t + (2 - d)^t)/(2 t (1 + d)),
  Assumptions -> 0 < d < 1/3 && 0 < t < 1]];

(* Numerical high-precision spot diagnostics, independent of Python roots.
   This is a bracketed bisection in u=(r-1)/(r+1), avoiding FindRoot's
   loss of scale for t extremely close to zero or one. *)
taus = {10^-12, 10^-6, 1/10, 1/2, 9/10, 1 - 10^-12};
kk[u_, tt_] := -(1 - tt) Log[u] + Log[tt - (2 - tt) u] - Log[(2 - tt) - tt u];
root[t0_] := Module[{tt = N[t0, 160], lo, hi, mid, u},
 lo = N[10^-140, 160]; hi = tt/(2 - tt) (1 - N[10^-140, 160]);
 guard["numeric_bracket", TrueQ[kk[lo, tt] > 0 && kk[hi, tt] < 0]];
 Do[mid = (lo + hi)/2; If[kk[mid, tt] > 0, lo = mid, hi = mid], {650}];
 u = (lo + hi)/2; (1 + u)/(1 - u)];
rows = Table[
  rt = root[t0]; dd = rt - 1;
  {N[t0, 20], N[rt, 40], N[t0 (rr - j) /. {t -> t0, d -> dd}, 30],
   TrueQ[1 < rt < 4/3]},
  {t0, taus}];
guard["numeric_root_interval_and_gap", And @@ (#[[3]] > 0 && #[[4]] & /@ rows)];
Export[FileNameJoin[{DirectoryName[$InputFileName], "global_twoatom_wolfram_rows.json"}],
 Association["version" -> "global-twoatom-wl-v1-20261105", "rows" -> (AssociationThread[{"tau", "r", "tau_gap", "in_interval"}, #] & /@ rows)], "RawJSON"];
Print["PASS all guards"];
