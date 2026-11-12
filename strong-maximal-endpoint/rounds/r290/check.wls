(* Integrate the missing-corner Laplace kernel in one variable. *)
Clear[a,b,l];
r=Integrate[Exp[-l a b],{a,0,1},Assumptions->l>0 && b>0];
If[!TrueQ[FullSimplify[r==(1-Exp[-l b])/(l b),Assumptions->l>0 && b>0]],Exit[1]];
Do[
 val=NIntegrate[(1-Exp[-v])/v,{v,0,1,L},WorkingPrecision->40];
 If[!(val<1+Log[L]),Exit[2]],
 {L,{2,8,32,128,512}}];
Print["symbolic integral and five corner bounds passed"];
