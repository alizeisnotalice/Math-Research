(* Exact Wolfram guards for the round-24 fixed-d cubic envelope identities. *)
Print["wolfram_version=",$Version];
Clear[q,d,del];
gdel=q (1-del/q)^3+(1-q) (1+del/(1-q))^3-1;
g=(d^2/4)*(3 q (1-q)+(d/2)(2q-1))/(q^2 (1-q)^2);
r=1-3q+3q^2;dq=3q(1-q)(1-2q)/r;
stationary=6q^3-3(3+d)q^2+3(1+d)q-d;
paramH=27q(1-q)(1-2q)^2(1-2q+2q^2)/(8r^3);
guards={
 TrueQ[FullSimplify[gdel==g/.d->2del,Assumptions->0<del<q<1]],
 TrueQ[FullSimplify[D[gdel,{q,2}]==6 del^2*((1-del/q)/q^3+(1+del/(1-q))/(1-q)^3),Assumptions->0<del<q<1]],
 TrueQ[FullSimplify[(D[gdel,q]/.q->del)==(3del-1)/(1-del)^3,Assumptions->0<del<1]],
 TrueQ[FullSimplify[(D[g,q]/.q->1/3)==81d^2(3d-2)/32]],
 TrueQ[FullSimplify[(D[g,q]/.q->1/2)==4d^3]],
 TrueQ[FullSimplify[D[g,d]==3d*(4q(1-q)+d(2q-1))/(8q^2(1-q)^2)]],
 TrueQ[FullSimplify[stationary==0/.d->dq]],
 TrueQ[FullSimplify[(g/.d->dq)==paramH]],
 TrueQ[FullSimplify[(g/.q->d/2)==d(4-d)/(2-d)^2]],
 TrueQ[FullSimplify[{dq/.q->1/3,paramH/.q->1/3}=={2/3,5/4}]]
 };
Print["guard_vector=",guards];
If[!And@@guards,Print["guard_failure"];Exit[1]];
Print["all_guards=True"];
