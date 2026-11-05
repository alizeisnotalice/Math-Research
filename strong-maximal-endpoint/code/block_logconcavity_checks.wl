(* Exact piecewise log-concavity audit for R_b at c=1. *)
Clear[x, r, ncurv, dom, piece, pieces, b3r, b3n];
piece[b_, j_] := Module[{rr, nn, dd},
  rr = Sum[Binomial[b, k] (x - k)^(b - 1), {k, 0, j}];
  nn = Factor[rr D[rr, {x, 2}] - D[rr, x]^2];
  dd = If[j < b, j <= x <= j + 1, x >= b];
  <|"piece" -> j, "forall_nonpositive" -> Resolve[ForAll[x, Implies[dd, nn <= 0]], Reals]|>
];
pieces = Association@Table[ToString[b] -> Table[piece[b, j], {j, 0, b}], {b, 4, 7}];
b3r = Sum[Binomial[3, k] (x - k)^2, {k, 0, 1}];
b3n = Factor[b3r D[b3r, {x, 2}] - D[b3r, x]^2];
result = <|"wolfram_version" -> $Version, "bases" -> pieces,
  "b3_piece_1_numerator" -> ToString[b3n, InputForm],
  "b3_at_x_1" -> (b3n /. x -> 1), "b3_at_x_101_over_100" -> (b3n /. x -> 101/100)|>;
Clear[cc,zz]; pp=1/(1+Exp[cc]);
lap=(1+Exp[-(1+zz)cc])/((1+zz)(1+Exp[-cc]));
res1=FullSimplify[D[lap,cc]-pp*(lap-Exp[-zz*cc]),Assumptions->{cc>0,zz>0}];
res2=FullSimplify[D[lap,{cc,2}]-pp*(2pp-1)*(lap-Exp[-zz*cc])-pp*zz*Exp[-zz*cc],Assumptions->{cc>0,zz>0}];
res3=FullSimplify[D[lap,{cc,3}]-pp*(1-6pp+6pp^2)*(lap-Exp[-zz*cc])-pp*(3pp-2)*zz*Exp[-zz*cc]+pp*zz^2*Exp[-zz*cc],Assumptions->{cc>0,zz>0}];
result=Join[result,<|"Laplace_kernel_derivative_residuals"->ToString[{res1,res2,res3},InputForm],"quadratic_surrogate_volume_constant"->ToString[N[Sqrt[Pi*(1+ProductLog[Exp[-1]])/2],80],InputForm]|>];
base=DirectoryName[DirectoryName[$InputFileName]];
Export[FileNameJoin[{base,"rounds","r027","symbolic_checks.json"}],result,"RawJSON"];
Print[<|"base_orders"->Keys[pieces],"all_piece_checks"->And@@Flatten[Table[Lookup[pieces[ToString[bb]],"forall_nonpositive"],{bb,4,7}]],"b3_at_x_1"->(b3n/.x->1),"kernel_residuals"->{res1,res2,res3}|>];
