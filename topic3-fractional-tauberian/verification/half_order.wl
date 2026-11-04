(* Half-order symmetry, classical normalization, and a nonsymmetric test. *)
ClearAll[a, b, z, d];
Print["Wolfram: ", $Version];
identity = Resolve[ForAll[{a, b},
  Abs[b + a] + Abs[-b + a] == 2 Max[Abs[b], Abs[a]]], Reals];
Print["Paired-face identity: ", identity];
If[identity =!= True, Exit[1]];
omega[d_] := Pi^(d/2)/Gamma[1 + d/2];
dimensions = {2, 3, 4, 8, 16, 32, 64};
rows = Table[{d, ToString[FullSimplify[(omega[d]/omega[d-1])^d], InputForm],
  ToString[N[(omega[d]/omega[d-1])^d, 75], InputForm],
  ToString[N[2^d/d!, 75], InputForm]}, {d, dimensions}];
Print["Dimension / ellipsoid coefficient / cube coefficient:"];
Scan[Print, rows];
(* For T=conv{(-1,-1),(2,-1),(-1,2)}, v=(1,1), V=9/2. *)
trianglePhi[a_] := a + (6 Max[1-a,0] + 3 Max[-2-a,0])/(9/2);
values = {trianglePhi[0], trianglePhi[1]};
Print["Nonsymmetric triangle: fixed scale, improved scale = ", values];
If[values =!= {4/3, 1}, Exit[1]];
result = <|"version" -> $Version, "pairedFaceIdentity" -> identity,
  "columns" -> {"n", "ellipsoidExact", "ellipsoid75digits", "cube75digits"},
  "rows" -> rows, "triangleFixedScale" -> "4/3", "triangleImprovedScale" -> "1",
  "scope" -> "Algebra and constants; Petty's inequality is a cited classical theorem."|>;
Export[FileNameJoin[{DirectoryName[$InputFileName], "..", "results", "half_order_wolfram.json"}],
  result, "RawJSON"];
