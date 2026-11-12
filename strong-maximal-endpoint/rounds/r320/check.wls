Clear[t,n];checks={Resolve[ForAll[{t,n},Implies[n>=1&&t>=2n,-1<=-1+(n-1)/t<=-1/2]],Reals]};
Do[AppendTo[checks,Sum[Binomial[n,b] 3^b/b!,{b,0,n}]<=Sum[(3n)^b/(b!)^2,{b,0,n}]],{n,{2,8,24,64}}];Print[checks];If[!And@@checks,Exit[1]];
