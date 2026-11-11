(* Independent Wolfram Language integer evaluation of every weight interval. *)
m = 256; d = 10^9;
v = Table[Floor[d (101^k - 100^k)/(101^k + 100^k)], {k, 0, 1024}];
mass = Map[Max[0, # - 1] &, Differences[v]];
sums = Table[Module[{ks},
    ks = Table[Quotient[128 (m + i) - (i + 1) j, m], {j, 1024}];
    mass . (v[[1 + Clip[ks, {0, 1024}]]])
    ], {i, 0, m - 1}];
proof = {128 (1/100 - 1/20000 + 1/3000000) < 637/500,
  Sum[(32/25)^k/k!, {k, 0, 8}] + (32/25)^9/9!/(1 - (32/25)/10) < 18/5,
  1 + 5/18 > 637/500,
  1/(3 Sqrt[3]) < 1/5,
  Min[sums] == 512686147030756979, Min[sums]/d^2 - 1/90 > 1/2};
Print[InputForm[{Min[sums], First[Ordering[sums]] - 1, proof}]];
If[!And @@ proof, Exit[1]];
