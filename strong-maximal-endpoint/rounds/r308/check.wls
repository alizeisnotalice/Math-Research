Clear[u,m,c,p,v];
checks={FullSimplify[Integrate[u(1+u^2)^(-m/2),{u,0,Infinity},Assumptions->m>2]==1/(m-2),Assumptions->m>2],FullSimplify[Integrate[Exp[-c p Sqrt[m] v],{v,0,Infinity},Assumptions->c>0&&p>0&&m>0]==1/(c p Sqrt[m]),Assumptions->c>0&&p>0&&m>0]};
Print[checks];If[!And@@checks,Exit[1]];
