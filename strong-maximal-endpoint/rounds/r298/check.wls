Clear[a,b,t];
checks={FullSimplify[a/b+(1-a/b)/(1+b t)==(1+a t)/(1+b t),Assumptions->0<=a<b<1&&t>=0],FullSimplify[a/b+(1-a/b)/(1-b)==(1-a)/(1-b),Assumptions->0<=a<b<1]};
Do[AppendTo[checks,Sum[T^j/j!,{j,0,d-1}]/(2 Sum[T^j/j!,{j,0,d-2}])>=T/(4(d-1))],{d,{2,4,8,16,32}},{T,{2d,d^2}}];
Print[{Length[checks],And@@checks}];If[!And@@checks,Exit[1]];
