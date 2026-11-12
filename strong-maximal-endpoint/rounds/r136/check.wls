(* Both log inequalities start at equality, with nonnegative derivative differences. *)
checks={FullSimplify[D[r-Log[1+r],r]>=0,r>=0],FullSimplify[D[Log[1+r]-r/(1+r),r]>=0,r>=0],FullSimplify[l/(2 Sqrt[l/(2a)])==Sqrt[a l/2],a>0&&l>0],Limit[Sqrt[l] Exp[-Sqrt[l]/2],l->Infinity]==0};Print[checks];If[!TrueQ[And@@checks],Exit[1]];
