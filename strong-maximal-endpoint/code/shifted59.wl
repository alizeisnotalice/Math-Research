(* Algebra/numerical receipt for r059.  This is corroboration, not a
   replacement for the analytic Gamma and coarea estimates in shifted59.md. *)
Clear[b, j, c, t, q, alpha, coefficient, supportLower]

q = Exp[c];
alpha = q/(1 + q);
(* Binomial weight for J=j times e^t Gamma_b(t-c j). *)
coefficient = Binomial[b, j] (1/(1 + q))^j (q/(1 + q))^(b - j) Exp[c j];
coefficientCheck = FullSimplify[coefficient == Binomial[b, j] alpha^b,
  Assumptions -> Element[{b, j}, Integers] && 0 <= j <= b && c > 0];

(* If lambda R_1 R_2 >= 1 and lambda <= 2 N=4b, while
   R_i <= 2^b t_i^(b-1)/Gamma[b], then t_1 t_2 is at least this value. *)
supportLower = (Gamma[b]^2/(4 b 4^b))^(1/(b - 1));
elementaryLower = ((b - 1)/E)^2 (4 b)^(-1/(b - 1)) 4^(-b/(b - 1));
constantLower = b^2/(128 E^2);

table = Table[
  <|"b" -> bb, "gamma_bound_holds" -> N[Gamma[bb] >= ((bb - 1)/E)^(bb - 1), 30],
    "elementary_over_constant" -> N[(elementaryLower/constantLower) /. b -> bb, 30],
    "exact_support_over_constant" -> N[(supportLower/constantLower) /. b -> bb, 30]|>,
  {bb, {6, 8, 16, 64, 256, 1024}}];

(* Exact implicit-differentiation audit for d(a): G(d(a),a)=1. *)
Clear[Gd, Ga, Gdd, Gad, Gaa, y, beta, gamma]
riccatiCheck = FullSimplify[
  -(Gaa + 2 Gad y + Gdd y^2)/Gd == -Gaa/Gd - 2 (Gad/Gd) y - (Gdd/Gd) y^2];

result = <|
  "wolfram_version" -> $Version,
  "coefficient_cancellation" -> coefficientCheck,
  "source_support_formula" -> ToString[InputForm[supportLower]],
  "claimed_elementary_constant" -> "b^2/(128 E^2), b>=6",
  "implicit_riccati_algebra" -> riccatiCheck,
  "table" -> table,
  "scope" -> "symbolic coefficient cancellation and numerical checks of the elementary support constant only; no sign or tail claim"
  |>;
out59=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/shifted59.json"];
If[!TrueQ[coefficientCheck && riccatiCheck],Exit[2]];
Export[out59, result, "JSON"];
Print[ExportString[result, "JSON", "Compact" -> False]];
