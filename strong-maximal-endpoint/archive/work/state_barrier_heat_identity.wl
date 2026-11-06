Clear[x,t,m,k,t0];f=Exp[-x^2/(2m t)]/Sqrt[2Pi m t];h=k(x^2/(m t)+2Log[t/t0]);
res=FullSimplify[(D[h f,t]-m D[h f,{x,2}]/2)/f,Assumptions->{m>0,t>0,t0>0}];
Print[InputForm[res]];
Print[FullSimplify[res-k(1+x^2/(m t))/t,Assumptions->{m>0,t>0,t0>0}]];
