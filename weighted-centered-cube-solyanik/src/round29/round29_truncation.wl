(* Wolfram independent check: density integration, exact moment and split guards. *)
$MaxExtraPrecision=10000; ClearAll[guard,mom,le,ta,u,n,g,q,a];
guard[name_,v_]:=If[TrueQ[v],Print["PASS ",name],Print["FAIL ",name];Exit[1]];
mom[n_,g_,q_]:=Module[{u},1/n+(n-1)/(n g) Integrate[u^(q-1/g-1),{u,1,Infinity},Assumptions->0<g<1&&q<1/g]];
le[n_,g_,q_,a_]:=Module[{u},1/n+(n-1)/(n g) Integrate[u^(q-1/g-1),{u,1,a},Assumptions->a>=1&&0<g<1]];
ta[n_,g_,q_,a_]:=Module[{u},(n-1)/(n g) Integrate[u^(q-1/g-1),{u,a,Infinity},Assumptions->a>=1&&0<g<1&&q<1/g]];
Do[guard["moment",FullSimplify[mom[n,g,q]==(n-g q)/(n(1-g q)),0<q<1/g]],{n,{2,3}},{g,{1/4,1/2,3/4}},{q,{1,6/5}}];
split=Flatten[Table[Module[{p=1/(1-(1-g)c)},FullSimplify[le[n,g,q,a]+ta[n,g,q,a]==mom[n,g,q]]],
 {n,{2,3}},{g,{1/4,1/2,3/4}},{c,{1/4,1/2,3/4,19/20}},{q,{0,1,1/(1-(1-g)c)}},{a,{1,2,257}}]];
guard["partial_plus_tail",And@@split];
(* Two exact rational diagnostics protecting the former gamma*q-1 exponent bug. *)
guard["rational_q0_A2",le[2,1/2,0,2]+ta[2,1/2,0,2]==mom[2,1/2,0]];
guard["rational_q1_A4",le[2,1/2,1,4]+ta[2,1/2,1,4]==mom[2,1/2,1]];
guard["log_partial_q2_A2",le[2,1/2,2,2]==1/2+Log[2]];
Print["PASS round29 Wolfram checks"];
