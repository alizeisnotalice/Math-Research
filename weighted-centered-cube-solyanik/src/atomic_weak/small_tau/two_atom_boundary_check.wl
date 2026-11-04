Print["Wolfram version: ", $Version];
Do[
  delta[u_] := Exp[-u/tau];
  f[u_] := ((1 + 2 delta[u])^tau - (1 - 2 delta[u])^tau +
      delta[u]^tau + (2 - delta[u])^tau)/(2 + 2 delta[u]);
  starts = {3 tau, 5 tau, 8 tau, 12 tau, 20 tau};
  values = Table[Quiet[FindMaximum[{f[u], tau Log[3] < u < 2},
      {{u, start}}, WorkingPrecision -> 80, AccuracyGoal -> 45, PrecisionGoal -> 45]],
      {start, starts}];
  sol = First@MaximalBy[values, First];
  Print["tau=", tau, ", max tau*R(two atoms)=", N[sol[[1]],30],
    ", delta=", N[delta[u /. sol[[2]]],30], ", u=", N[u /. sol[[2]],30]],
  {tau, {1/64, 1/128, 1/256}}
]
