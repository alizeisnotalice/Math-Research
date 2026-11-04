(* Round 26 uniform two-atom versus single-atom symbolic checks.
   Run: wolframscript -file uniform_two_atom_verify.wl *)
Clear[t, d, r];
n = (r - 1)^t + (r + 1)^t;
f = FullSimplify[r D[n, r] - n];
f43 = FullSimplify[f /. r -> 4/3];
jExpected = 3^(-t) (4 t - 1 + (4 t - 7) 7^(t - 1));
j = (d^t + (2 + d)^t)/(2 t (1 + d));
rr = ((1 + 2 d)^t - (1 - 2 d)^t + d^t + (2 - d)^t)/(2 t (1 + d));
diff = FullSimplify[rr - j];
integralDiff = ((1 + 2 d)^t - (1 - 2 d)^t - (2 + d)^t + (2 - d)^t)/(2 t (1 + d));
h = Log[(7 - 4 t)/(4 t - 1)] + (t - 1) Log[7];
hExpected = Log[7] - 24/((7 - 4 t) (4 t - 1));
assF = 0 < t < 1 && r > 1;
assD = 0 < t < 1 && 0 < d < 1/3;
checks = {
  TrueQ[FullSimplify[f43 == jExpected]],
  TrueQ[FullSimplify[diff == integralDiff, assD]],
  TrueQ[FullSimplify[D[f, r] < 0, assF]],
  TrueQ[FullSimplify[D[h, t] == hExpected, 1/4 < t < 1]],
  TrueQ[FullSimplify[Log[7] < 8/3]]
  };
Print["Wolfram version: ", $Version];
Print["F(4/3) / 3^(-t) = ", FullSimplify[f43/3^(-t)]];
Print["R_delta-J(1+delta) = ", diff];
Print["Integral form difference = ", FullSimplify[diff - integralDiff]];
Print["F'(r) = ", FullSimplify[D[f, r]]];
Print["Exact guard vector = ", checks];
If[!And @@ checks, Print["UNIFORM CERTIFICATE FAILED"]; Abort[]];
Do[
  root = r /. FindRoot[(f /. t -> tt) == 0, {r, 1.05, 1.32}];
  delta = root - 1;
  Print["tau=", tt, "; r*=", N[root, 20], "; tau*R-J=", N[tt (rr-j) /. {t -> tt, d -> delta}, 20]],
  {tt, {1/10, 1/4, 1/2, 3/4, 9/10}}];
Print["UNIFORM CERTIFICATE: PASS"];
