Print["VERSION=", $Version];
Print["TYPE=exact symbolic audit; not numerical proof of all-scale geometry"];
Print["RATIONAL_WIDTH_MARGIN=", FullSimplify[
 (e + 4 a^2) - 4 (a + e/(16 a))^2]];
Print["MARGIN_IDENTITY=", FullSimplify[
 (e + 4 a^2) - 4 (a + e/(16 a))^2 == e/2-e^2/(64 a^2)]];
Print["MARGIN_POSITIVE=", FullSimplify[
 (e + 4 a^2) - 4 (a + e/(16 a))^2 > 0,
 Assumptions -> a > 0 && 0 < e <= a^2]];
Print["ANNULAR_LAYER_AREA=", FullSimplify[2 (a-a/2) e/(16 a)]];
Print["FIRST_LAYER_AREA=", FullSimplify[2 a e/(16 a)]];
Print["TOTAL_LAYER_AREA=", FullSimplify[e/8+(n-1)e/16]];
Print["PREFIX_TELESCOPE=", FullSimplify[
 Sum[4 (2^(2 (j-n))-2^(2 (j-1-n))), {j, 2, k}]+4*2^(2 (1-n)),
 Assumptions -> Element[{n,k},Integers] && n >= k >= 1]];
Print["NORMALIZED_LOSS_RATIO=", FullSimplify[(area/total)/(e/total), Assumptions->total>0 && e>0]];
Print["L1_FIXED_RADIUS_MARGIN=",FullSimplify[4 a^2+e-4(a+e/(16a))^2]];
Print["L1_LAYER_AREA=",FullSimplify[a*(e/(16a))/2]];
Print["L1_SUPPORT_FITS_UNIT_BOX=", FullSimplify[(1+2*(e/64))^2<4+e,Assumptions->0<e<=1]];
Print["PROBABILITY_FAMILY_EPSILON_RATIO=",Factor[(1-(1+q)^(-2))/(16 q/(1+q)^2)]];
Print["PROBABILITY_FAMILY_EPSILON_RATIO_LT_ONE=",FullSimplify[(2+q)/16<1,Assumptions->0<q<=1/64]];
Print["PROBABILITY_FAMILY_RIGHT_SUPPORT_LT_HALF=",FullSimplify[
 1/(2(1+q))+(1-(1+q)^(-2))/(64/(2(1+q)))<1/2,
 Assumptions->0<q<=1/64]];
Print["EPSILON_SANDWICH=",FullSimplify[q<=1-(1+q)^(-2)<=2q,Assumptions->0<q<=1/64]];
Print["THRESHOLD_MINIMUM_IDENTITY=",FullSimplify[
 A/(1-th)+B/th-(Sqrt[A]+Sqrt[B])^2==
 (Sqrt[A]*th-Sqrt[B]*(1-th))^2/(th*(1-th)),
 Assumptions->A>0&&B>0&&0<th<1]];
Print["EXCESS_IDENTITY=",FullSimplify[
 (Sqrt[K*(1-e)]+Sqrt[e])^2-K==2*Sqrt[K*e*(1-e)]-(K-1)*e,
 Assumptions->K>=1&&0<e<1]];
Print["SUBCRITICAL_BRANCH_MONOTONICITY=",FullSimplify[
 D[(1-d/r)^2,r]>0 && D[m/(4r^2),r]<0,
 Assumptions->m>0&&d>0&&r>d]];
Print["SUBCRITICAL_CROSSING_VALUE=",FullSimplify[
 {(1-d/r)^2,m/(4r^2)}/.r->d+Sqrt[m]/2,
 Assumptions->m>0&&d>0]];
