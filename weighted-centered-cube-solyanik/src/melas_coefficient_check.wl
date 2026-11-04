(* Exact algebraic check of the published one-dimensional weak-type value.
   This does not reprove Melas's theorem. *)
Print["Wolfram version: ", $Version];
weakMelas = (11 + Sqrt[61])/12;
Print["published quadratic identity: ", FullSimplify[12 weakMelas^2 - 22 weakMelas + 5 == 0]];
Print["m, fixed-Lebesgue residue, whole-weight-family residue (n=1)"];
Do[Print[{m, N[weakMelas (1 - 1/m)^(m-1), 40], N[2 (1 - 1/m)^(m-1), 40]}], {m, {2,3,4,8}}];
