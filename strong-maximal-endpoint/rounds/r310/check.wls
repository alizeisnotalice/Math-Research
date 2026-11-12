Clear[v,m,x];checks={FullSimplify[2/(1-v)-1==(1+v)/(1-v),Assumptions->0<v<1],FullSimplify[Log[(1+Tanh[m/2])/(1-Tanh[m/2])]==m,Assumptions->m>0]};Print[checks];If[!And@@checks,Exit[1]];
