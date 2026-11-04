(* Reproducible exact-arithmetic checks for the polytope tangent duality. *)

ClearAll[cross2, hull2, polygonArea, qCandidates, tangentPolygon, polarPolygon,
  reportCase, formatPoint];

cross2[a_, b_, c_] := Det[{b - a, c - a}];

hull2[points_] := Module[{s, append, lower, upper},
  s = SortBy[DeleteDuplicates[points], {#[[1]], #[[2]]} &];
  If[Length[s] <= 2, Return[s]];
  append[h_, p_] := Module[{hh = h},
    While[Length[hh] >= 2 && cross2[hh[[-2]], hh[[-1]], p] <= 0,
      hh = Most[hh]];
    Append[hh, p]
  ];
  lower = Fold[append, {}, s];
  upper = Fold[append, {}, Reverse[s]];
  Most[lower] ~Join~ Most[upper]
];

polygonArea[p_] := Abs[Total[Det /@ Partition[Append[p, First[p]], 2, 1]]]/2;

qCandidates[areas_, supports_, volume_, beta_] := Module[
  {m = Length[areas], target, candidates = {}, indices, rest, bits, known, qi, q},
  target = (2 - beta) volume;
  Do[
    rest = Delete[Range[m], i];
    Do[
      known = Total[MapThread[#1 #2 #3 &, {areas[[rest]], supports[[rest]], bits}]];
      qi = (target - known)/(areas[[i]] supports[[i]]);
      If[0 <= qi <= 1,
        q = ConstantArray[0, m]; q[[rest]] = bits; q[[i]] = qi;
        AppendTo[candidates, q]
      ],
      {bits, Tuples[{0, 1}, m - 1]}
    ],
    {i, Range[m]}
  ];
  DeleteDuplicates[candidates]
];

tangentPolygon[areas_, normals_, supports_, volume_, beta_] := Module[{q},
  q = qCandidates[areas, supports, volume, beta];
  hull2[(-Transpose[normals].(areas #))/volume & /@ q]
];

polarPolygon[p_] := Module[{pairs, raw},
  pairs = Partition[Append[p, First[p]], 2, 1];
  raw = LinearSolve[#, {1, 1}] & /@ pairs;
  hull2[raw]
];

formatPoint[p_] := ToString[p, InputForm];

reportCase[name_, areas_, normals_, supports_, volume_, beta_] := Module[
  {p, s, parea, sarea, q, moment, testVectors, primal, reflected, dual, maxGap, reflectionGap},
  q = qCandidates[areas, supports, volume, beta];
  p = tangentPolygon[areas, normals, supports, volume, beta];
  s = polarPolygon[p];
  parea = polygonArea[p]; sarea = polygonArea[s];
  moment = Table[Total[areas supports q[[j]]], {j, Length[q]}];
  testVectors = {{1, 0}, {0, 1}, {1, 1}, {2, -1}, {-1, 2}};
  primal[v_] := MinValue[
    {(2 - beta) a + Total[areas MapThread[Max[-#1 - a #2, 0] &, {normals.v, supports}]]/volume,
      a \[Element] Reals}, a
  ];
  reflected[v_] := MinValue[
    {-beta a + Total[areas MapThread[Max[#1 + a #2, 0] &, {normals.v, supports}]]/volume,
      a \[Element] Reals}, a
  ];
  dual[v_] := Max[(# . v) & /@ p];
  maxGap = Max[Abs[primal[#] - dual[#]] & /@ testVectors];
  reflectionGap = Max[Abs[primal[#] - reflected[#]] & /@ testVectors];
  {
    "case=" <> name <> ", beta=" <> ToString[beta, InputForm],
    "volume=" <> ToString[volume, InputForm],
    "targetMoment=" <> ToString[(2 - beta) volume, InputForm],
    "candidateMoments=" <> ToString[DeleteDuplicates[moment], InputForm],
    "PVertices=" <> StringRiffle[formatPoint /@ p, "; "],
    "areaP=" <> ToString[parea, InputForm],
    "SVertices=" <> StringRiffle[formatPoint /@ s, "; "],
    "areaS=" <> ToString[sarea, InputForm],
    "areaSOverVolume=" <> ToString[Together[sarea/volume], InputForm],
    "maxPrimalDualGapFiveExactDirections=" <> ToString[maxGap, InputForm],
    "maxReflectionGapFiveExactDirections=" <> ToString[reflectionGap, InputForm],
    ""
  }
];

scriptDirectory = DirectoryName[$InputFileName];
projectDirectory = DirectoryName[scriptDirectory];
resultDirectory = FileNameJoin[{projectDirectory, "results"}];
If[!DirectoryQ[resultDirectory], CreateDirectory[resultDirectory, CreateIntermediateDirectories -> True]];

squareAreas = {1, 1, 1, 1};
squareNormals = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};
squareSupports = {1/2, 1/2, 1/2, 1/2};
squareVolume = 1;

triangleAreas = {3, 3 Sqrt[2], 3};
triangleNormals = {{0, -1}, {1/Sqrt[2], 1/Sqrt[2]}, {-1, 0}};
triangleSupports = {1, 1/Sqrt[2], 1};
triangleVolume = 9/2;

lines = Join[
  {"Wolfram Language version=" <> $Version,
   "All polytope calculations below use exact rational/algebraic arithmetic.", ""},
  Flatten[reportCase["square[-1/2,1/2]^2", squareAreas, squareNormals, squareSupports, squareVolume, #] & /@ {1/2, 1, 3/2}],
  Flatten[reportCase["triangle[(-1,-1),(2,-1),(-1,2)]", triangleAreas, triangleNormals, triangleSupports, triangleVolume, #] & /@ {1/2, 1, 3/2}]
];

diskDefect = Integrate[-Cos[z], {z, Pi/2, 3 Pi/2}];
If[diskDefect =!= 2, Exit[1]];
lines = Join[lines, {"", "unitDisk beta=1: directionalDefect=" <> ToString[diskDefect, InputForm],
  "unitDisk beta=1: areaSOverVolume=" <> ToString[Pi^2/4, InputForm]}];

Export[FileNameJoin[{resultDirectory, "polytope_wolfram_results.txt"}], StringRiffle[lines, "\n"], "Text"];
Print[StringRiffle[lines, "\n"]];
