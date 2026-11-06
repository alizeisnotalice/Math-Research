(* r055: symbolic audit of the capped Phi_N dual density q. *)
out = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[-1]],
  "/tmp/dilation55.json"];
Clear[n,h,a];
phi[n_, h_] := h (1 + Log[h])^(n - 1);
ap[n_, h_] := (1 + Log[h])^(n - 2) (n + Log[h]);
phiDerivativeCheck = FullSimplify[D[phi[n, h], h] == ap[n, h], h > 0 && n > 1];
qpar[n_, h_] := (1 + Log[h])^(n - 2) (h + n - h n + Log[h]);
qfactor = qpar[n, h];
derivCheck = FullSimplify[D[qpar[n, h], h] == (1 - h) D[ap[n, h], h], h > 1 && n > 1];
qOverA = FullSimplify[qpar[n, h]/ap[n, h], h > 1 && n > 1];
n0 = 64;
hStar = h /. FindRoot[(n0 - 1) h == n0 + Log[h], {h, 1 + 1/n0}, WorkingPrecision -> 100];
aStar = N[ap[n0, hStar], 90];
phi2 = N[phi[n0, 2], 90];
aCap = N[ap[n0, 2], 90];
qCap = N[phi2 - aCap, 90];
rootChecks = {
  N[qpar[n0, hStar], 70],
  N[(n0 - 1) hStar - n0 - Log[hStar], 70],
  N[qpar[n0, 1], 70],
  N[qCap, 70]
};
(* Exact elementary moments under two independent unit exponentials. *)
t1mom = Integrate[(x - y)^2 Exp[-x - y], {x, 0, Infinity}, {y, 0, Infinity},
  Assumptions -> True];
linmom = Integrate[(x + y) Exp[-x - y], {x, 0, Infinity}, {y, 0, Infinity},
  Assumptions -> True];
tailD = FullSimplify[2 Integrate[d^2 Exp[-d]/2, {d, ell, Infinity}], ell > 0];
tailS = FullSimplify[Integrate[s^2 Exp[-s], {s, ell, Infinity}], ell > 0];
res = <|
  "wolfram" -> $Version,
  "symbolic_phi_derivative_check" -> phiDerivativeCheck,
  "symbolic_phi_prime" -> ToString[ap[n, h], InputForm],
  "symbolic_q_parametric" -> ToString[qfactor, InputForm],
  "q_parametric_derivative_check" -> derivCheck,
  "q_over_a_parametric" -> ToString[qOverA, InputForm],
  "q_sign_rule" -> "q=a on [0,1]; q=1 on [1,N]; q(phi'(h)) has the sign of N+log(h)-(N-1)h for 1<h<2; q=Phi(2)-a for a>=Phi'(2).",
  "N64_h_zero" -> ToString[N[hStar, 80], InputForm],
  "N64_a_zero" -> ToString[aStar, InputForm],
  "N64_phi2" -> ToString[phi2, InputForm],
  "N64_a_cap" -> ToString[aCap, InputForm],
  "N64_q_at_cap" -> ToString[qCap, InputForm],
  "N64_root_residuals" -> (ToString[#, InputForm]& /@ rootChecks),
  "unit_exp_D_squared" -> ToString[t1mom, InputForm],
  "unit_exp_tsum" -> ToString[linmom, InputForm],
  "unit_exp_D_squared_tail" -> ToString[tailD, InputForm],
  "unit_exp_tsum_tail" -> ToString[tailS, InputForm]
|>;

(* Independent pullback differentiation; N=2b. *)
Clear[c,theta,x1,x2,bb];
lx=Log[c^2-theta^2]-(c+theta)x1-(c-theta)x2;
lnu=bb(Log[c+theta]+Log[c-theta])-bb(Log[1+Exp[-c-theta]]+Log[1+Exp[-c+theta]])-(c+theta)x1-(c-theta)x2;
pp=1/(1+Exp[c]);vv=1+c^2 pp(1-pp);
checksPull=<|
"first_variations_equal"->TrueQ[FullSimplify[(D[lx,theta]-D[lnu,theta])/.theta->0,Assumptions->c>0]===0],
"xi_second"->TrueQ[FullSimplify[(D[lx,{theta,2}]+D[lx,theta]^2)/.theta->0,Assumptions->c>0]==(x1-x2)^2-2/c^2],
"nu_second"->TrueQ[FullSimplify[((D[lnu,{theta,2}]+D[lnu,theta]^2)/.theta->0)-((x1-x2)^2-2bb vv/c^2),Assumptions->c>0]===0]|>;
If[!TrueQ[phiDerivativeCheck]||!TrueQ[derivCheck]||!And@@Values[checksPull],Print[checksPull];Exit[2]];
res["pullback_checks"]=checksPull;

Export[out, res, "RawJSON"]; Print[res];
