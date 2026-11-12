Clear[s,r,v,z,u];
checks={FullSimplify[D[r^(s-1)+(s-1)v^s/r-s v^(s-1),v]==s(s-1)v^(s-2)(v/r-1),Assumptions->r>0&&s>2&&v>0],FullSimplify[Limit[Exp[(1-u)z]+(1/u-1)Exp[-u z]-1/u,u->0]==Exp[z]-1-z,Assumptions->Element[z,Reals]]};Print[checks];If[!And@@checks,Exit[1]];
