$MaxExtraPrecision = 1000;
Print["Wolfram version: ", $Version];
(* Formal Taylor audit after eliminating d from G(d,v,a)=1 at (0,0,c).
   gd>0, gaa<0; K is the mixed coefficient after radial elimination. *)
Clear[gd, gdd, gdvv, gvvvv, gaa, gad, gavv, d2, d4f, a2, kappa, delta];
d2 = -gvv/(2 gd);
kappa = FullSimplify[gavv + 2 gad d2];
a2 = FullSimplify[-kappa/(2 gaa)];
d4f = FullSimplify[-(gdd d2^2/2 + gdvv d2/2 + gvvvv/24)/gd];
delta = FullSimplify[-(gaa a2^2/2 + gad d2 a2 + gavv a2/2)/gd];
Print[<|
 "K" -> kappa,
 "a2" -> a2,
 "d2" -> d2,
 "d4Fixed" -> d4f,
 "d4OptMinusFixed" -> delta,
 "deltaEqualsKFormula" -> TrueQ[FullSimplify[delta == kappa^2/(8 gaa gd)]],
 "signUnderGaaNegativeGdPositive" -> "K^2/(8 gaa gd) <= 0"
 |>];

checks=<|"radially_corrected_K"->TrueQ[FullSimplify[kappa==gavv-gad gvv/gd,Assumptions->gd>0]],"radial_delta"->TrueQ[FullSimplify[delta==kappa^2/(8gaa gd),Assumptions->gd>0&&gaa<0]],"negative_quartic_coefficient"->TrueQ[FullSimplify[kappa^2/(8gaa gd)<=0,Assumptions->gd>0&&gaa<0&&Element[{gavv,gad,gvv},Reals]]]|>;
If[!And@@Values[checks],Print[checks];Exit[2]];
out=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/geometry57.json"];
Export[out,<|"wolfram"->$Version,"checks"->checks,"all_pass"->True,"scope"->"Exact Taylor algebra only; analytic hypotheses in geometry57.md"|>,"RawJSON"];
