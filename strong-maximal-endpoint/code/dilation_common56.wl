(* r056: common-parameter pullback algebra, no numerical optimization. *)
out = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[-1]], "/tmp/dilation_common56.json"];
Clear[c,n,y,p,m,z,lambda,A,B,rho,ellx,ellnu];
p[c_] := 1/(1 + Exp[c]);
m[c_] := 1 + c p[c];
z[c_] := 1 + Exp[-c];
logK = (n - 2) Log[c] - n Log[z[c]];
rhoExpr = (n m[c] - 2)/c;
logKCheck = FullSimplify[D[logK, c] == rhoExpr, c > 0 && n > 1];
ellx = 2/c - y;
ellnu = n m[c]/c - y;
ellDiffCheck = FullSimplify[ellnu - ellx == rhoExpr, c > 0 && n > 1];
ellxSecond = FullSimplify[D[ellx, c], c > 0];
ellnuSecond = FullSimplify[D[ellnu, c], c > 0 && n > 1];
vExpr = 1 + c^2 p[c] (1 - p[c]);
ellnuSecondCheck = FullSimplify[ellnuSecond == -n vExpr/c^2, c > 0 && n > 1];
formalC = FullSimplify[(0 - B) - rhoExpr (1 - 1) == -B];
schurCheck = FullSimplify[(-B/A) (-lambda B) == lambda B^2/A, A > 0];
res = <|
 "wolfram" -> $Version,
 "logK_derivative_check" -> logKCheck,
 "logR_common_derivative" -> ToString[rhoExpr, InputForm],
 "score_difference_check" -> ellDiffCheck,
 "xi_log_derivative" -> ToString[ellx, InputForm],
 "nu_log_derivative" -> ToString[ellnu, InputForm],
 "xi_log_second" -> ToString[ellxSecond, InputForm],
 "nu_log_second_check" -> ellnuSecondCheck,
 "formal_C_over_lambda_equals_minus_B" -> formalC,
 "schur_sign_check" -> schurCheck,
 "common_formula" -> "P'=Integral[q Lxi dxi]; P''=Integral[q(Lxi^2-2/c^2)dxi]+lambda B^2/A, eta=-B/A, A=E_nu[a h'(a)], B=E_nu[H Lnu]."
|>;
If[!And@@{logKCheck,ellDiffCheck,ellnuSecondCheck,formalC,schurCheck},Print[res];Exit[2]];
Export[out,res,"RawJSON"]; Print[res];
