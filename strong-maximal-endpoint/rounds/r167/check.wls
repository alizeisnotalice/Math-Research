tests={FullSimplify[Integrate[g-2h u/b,{u,0,b g/(2h)}]==b g^2/(4h),g>0&&b>0&&h>0],FullSimplify[2Integrate[x,{x,0,d/(2h)}]h==d^2/(4h),d>0&&h>0]};Print[tests];If[!TrueQ[And@@tests],Exit[1]];
