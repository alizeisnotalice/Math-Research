(* R078 narrow exact scalar checks. No actual input geometry encoded. *)
ClearAll[b, v, c, f, aa, a, qq, vv, mm, a0, n, rr, q];
Print["Version: ", $Version];
theta = b v/c;
weighted = Sum[theta^(9 - h) (c f[h] - b v f[h - 1]), {h, 1, 9}];
telescoping = Factor[weighted - (c f[9] - b v theta^8 f[0])];
rq = aa/(aa - a (1 + qq));
reverseRatio = Factor[(1 + qq) (vv + mm) - aa mm/a -
   (aa - a (1 + qq))/a vv (rq - (1 + mm/vv))];
endpoint = Factor[b v theta^8 (a0 v a)/a];
cdfNormalization = FullSimplify[v (v rr)^(n - 1)/v^n - rr^(n - 1),
  Assumptions -> {v > 0, rr > 0, Element[n, Integers], n >= 2}];
Print["Telescoping residual: ", telescoping];
Print["Reverse-ratio residual: ", reverseRatio];
Print["CDF endpoint-normalization residual: ", cdfNormalization];
Print["Endpoint using f0 <= A0 v a S: ", endpoint, " S"];
Print["Prefix exponents: terminal9, h8=", 1 + 8, ", h7=", 2 + 8,
  ", h6 little-o=", 3 + 6, ", h<=5 at least=", 4 + 6];
beta[q_] := q (2 q - 1)/(q - 1);
Print["Bad short q3: ", 7 + 3];
Print["Bad long q3 Young remainder: ", beta[3] + 3];
Print["Bad long q4 Young remainder: ", beta[4] + 3];
Print["Good q3 intrinsic exponent: ", 3*3,
  "; improved gap exponent: ", 2 + 9];
Print["Anchor/face9 times v: ", 9 + 1];
passed = TrueQ[telescoping == 0] && TrueQ[reverseRatio == 0] &&
  TrueQ[cdfNormalization == 0] &&
  TrueQ[Factor[endpoint - a0 b^9 v^10/c^8] == 0];
Print["All exact assertions passed: ", passed];
If[! passed, Exit[1]];
