(* Exact Wolfram Language check for the integer-mass critical identity. *)
Print["Wolfram version: ", $Version];
cases = {{2, 1/7, 17}, {3, 2/9, 29}, {4, 1/3, 41}, {8, 1/5, 97}};
Do[
  {m, f, n} = oneCase;
  c = (m - 1)/m;
  p = Join[{c + f/n}, ConstantArray[c, m - 1]];
  excess = Together[Times @@ p - c^m];
  ok = FullSimplify[excess == c^(m - 1) f/n && excess > 0];
  Print["m=", m, ", f=", f, ", N=", n, ", excess=", excess,
    ", identity-and-positivity=", ok],
  {oneCase, cases}
]
