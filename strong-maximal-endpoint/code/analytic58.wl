(* r058: symbolic checks for the radial-curvature and density identities.
   WolframScript input; no numerical/asymptotic inference is made here. *)
Clear["Global`*"];

(* A generic positive one-block density is f=Exp[-t+r(t,a)], so
   k=r_t and s=r_a are the likelihood scores of R=Exp[r]. *)
r1 = r[t1, a]; r2 = r[t2, a];
p = Exp[-t1 - t2 + r1 + r2];
lv[x_] := -D[x, t1] + D[x, t2];
q = FullSimplify[lv[lv[p]]/p];
qa = FullSimplify[D[lv[lv[p]], a]/p];
k1 = D[r1, t1]; k2 = D[r2, t2];
s1 = D[r1, a]; s2 = D[r2, a];
qExpected = (k1 - k2)^2 + D[k1, t1] + D[k2, t2];
qaExpected = FullSimplify[(s1 + s2) qExpected +
   2 (k1 - k2) (D[s1, t1] - D[s2, t2]) +
   D[s1, {t1, 2}] + D[s2, {t2, 2}]];

(* Smooth-region version of the one-integration-by-parts response formula.
   At a lower KKT level this is read distributionally; B includes that
   coarea contribution. *)
vv[x_] := D[x, t1] - D[x, t2];
z = lam Exp[r1 + r2]; h0 = h[z];
hv = FullSimplify[vv[h0]];
kminus = k1 - k2; splus = s1 + s2;
tdiffs = D[s1, t1] - D[s2, t2];
hvExpected = z h'[z] kminus;
ibpPointwise = FullSimplify[vv[hv splus p] == vv[hv] splus p + hv vv[splus p]];
responseCentering = FullSimplify[
  -hv (tdiffs + splus kminus) p + beta hv kminus p ==
  -hv (tdiffs + (splus - beta) kminus) p];

(* Along G(rho(a),0,a)=1, rho'(c)=0 when G_a=0.  The coefficient
   kappa(a)=-G_vv/(2G_d) has the claimed derivative at c. *)
rhoPrime = -ga/gd;
kappaPrimeAtBase = FullSimplify[-((gavv + gdvv rhoPrime) gd -
      gvv (gad + gdd rhoPrime))/(2 gd^2) /. ga -> 0];
kappaExpected = FullSimplify[-(gavv - gad gvv/gd)/(2 gd)];

checks = <|
  "density_angular_second" -> TrueQ[FullSimplify[q == qExpected]],
  "density_action_angular_second" -> TrueQ[FullSimplify[qa == qaExpected]],
  "smooth_response_first_variation" -> TrueQ[FullSimplify[hv == hvExpected]],
  "one_ibp_product_identity" -> TrueQ[ibpPointwise],
  "centered_response_algebra" -> TrueQ[responseCentering],
  "radial_curvature_derivative" -> TrueQ[FullSimplify[kappaPrimeAtBase == kappaExpected]],
  "q_formula" -> ToString[InputForm[qExpected]],
  "qa_formula" -> ToString[InputForm[qaExpected]],
  "kappa_prime_formula" -> ToString[InputForm[kappaExpected]],
  "wolfram_version" -> $Version
|>;
Print[checks];
If[!And @@ Values[KeyTake[checks, {"density_angular_second", "density_action_angular_second", "smooth_response_first_variation", "one_ibp_product_identity", "centered_response_algebra", "radial_curvature_derivative"}]], Exit[1]];

out=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/analytic58.json"];
Export[out,checks,"RawJSON"];
