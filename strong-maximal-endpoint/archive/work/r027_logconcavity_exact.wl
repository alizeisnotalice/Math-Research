(* Exact piecewise log-concavity audit for R_b at c=1. *)
Clear[x, r, ncurv, dom, piece, pieces, b3r, b3n];
piece[b_, j_] := Module[{rr, nn, dd},
  rr = Sum[Binomial[b, k] (x - k)^(b - 1), {k, 0, j}];
  nn = Factor[rr D[rr, {x, 2}] - D[rr, x]^2];
  dd = If[j < b, j <= x <= j + 1, x >= b];
  <|"piece" -> j, "curvature_numerator" -> ToString[nn, InputForm],
    "forall_nonpositive" -> Resolve[ForAll[x, Implies[dd, nn <= 0]], Reals]|>
];
pieces = Association@Table[ToString[b] -> Table[piece[b, j], {j, 0, b}], {b, 4, 7}];
b3r = Sum[Binomial[3, k] (x - k)^2, {k, 0, 1}];
b3n = Factor[b3r D[b3r, {x, 2}] - D[b3r, x]^2];
result = <|"wolfram_version" -> $Version, "bases" -> pieces,
  "b3_piece_1_numerator" -> ToString[b3n, InputForm],
  "b3_at_x_1" -> (b3n /. x -> 1)|>;
Export["/Users/zhengzhihao/Desktop/T/work/r027_logconcavity_exact_output.json", result, "RawJSON"];
Print[result];
