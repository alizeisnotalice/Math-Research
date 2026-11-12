checks={Limit[Exp[4Log[l]^2/l],l->Infinity]==1,Limit[l^(1/4)/Log[l],l->Infinity]==Infinity};Print[checks];If[!TrueQ[And@@checks],Exit[1]];
