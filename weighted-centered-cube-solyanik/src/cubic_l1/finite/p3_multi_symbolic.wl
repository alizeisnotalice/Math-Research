(* Exact independent symbolic checks used by the finite-vector audit. *)
Print["wolfram_version=", $Version];
Clear[q, a, b, delta, x, y, qs, endpoint, c3];
qs=(3-Sqrt[5])/2; c3=(11+5 Sqrt[5])/8;
delta=q (1-a); b=(1-q a)/(1-q);
r=FullSimplify[(q a^3+(1-q)b^3-1)/(4 delta^2),Assumptions->0<q<1&&0<a<1];
Print["binary_ratio=",r];
endpoint=FullSimplify[r/.a->0,Assumptions->0<q<1/2];
Print["endpoint=",endpoint];
Print["derivative_factor=",Factor[D[endpoint,q]]];
Print["critical_value=",FullSimplify[endpoint/.q->qs]];
Print["critical_matches_C3=",FullSimplify[(endpoint/.q->qs)==c3]];
Print["endpoint_minus_C3_factor=",Factor[Together[endpoint-c3],Extension->Sqrt[5]]];
Print["two_point_jensen_gap=",Factor[x^3+y^3-2((x+y)/2)^3]];
Print["jensen_gap_nonnegative=",FullSimplify[x^3+y^3>=2((x+y)/2)^3,Assumptions->x>0&&y>0]];
ratioExpected=(2+a-q-2*a*q)/(4*(1-q)^2*q);
derivativeExpected=(1-3*q+q^2)/(2*(-1+q)^3*q^2);
dIdentity=FullSimplify[q*(1-a)+(1-q)*(b-1)==2*q*(1-a),Assumptions->0<q<1&&0<a<1];
jensenIdentity=FullSimplify[x^3+y^3-2((x+y)/2)^3==3*(x-y)^2*(x+y)/4];
guards={
  TrueQ[FullSimplify[r==ratioExpected,Assumptions->0<q<1&&0<a<1]],
  TrueQ[FullSimplify[D[endpoint,q]==derivativeExpected,Assumptions->0<q<1/2]],
  TrueQ[FullSimplify[(endpoint/.q->qs)==c3]],
  TrueQ[dIdentity],
  TrueQ[jensenIdentity],
  TrueQ[FullSimplify[((1+x)^3+(1-x)^3)/2-1==3*x^2]],
  TrueQ[FullSimplify[x^3+y^3>=2((x+y)/2)^3,Assumptions->x>0&&y>0]]
  };
Print["guard_vector=",guards];
If[!And@@guards,Print["symbolic_guard_failure"];Exit[1]];
Print["all_symbolic_guards=True"];
