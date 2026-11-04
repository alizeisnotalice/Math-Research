(* Independent symbolic check of the strict margin used in (GC.10). *)
Clear[d, f, eps, q, m, c];
Print[{"reference-slack", FullSimplify[
   (1 - d) q + d f - (1 - d) q,
   Assumptions -> 0 < d < 1 && f > 0 && q > 0]}];
Print[{"perturbation-positive", FullSimplify[
   d f - (2 - d) eps > 0,
   Assumptions -> 0 < d <= 1 && f > 0 && 0 < eps < d f/4]}];
Print[{"equal-exponent-gap", FullSimplify[
   (c + c f/q)^m - c^m,
   Assumptions -> Element[m, Integers] && m >= 2 && c > 0 && f > 0 && q > 0]}];
Print[{"equal-exponent-gap-positive-m=2,3,4,8", Table[FullSimplify[
   (c + c f/q)^mm - c^mm > 0,
   Assumptions -> c > 0 && f > 0 && q > 0], {mm, {2, 3, 4, 8}}]}];
