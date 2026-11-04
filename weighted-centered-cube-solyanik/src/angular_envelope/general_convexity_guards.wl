Print[$Version];
Clear[q24,delta24,p24];
assume24=0<delta24<q24<1&&p24>1;
a24=1-delta24/q24;b24=1+delta24/(1-q24);
g24=q24*a24^p24+(1-q24)*b24^p24-1;
f24[x_]:=p24*x^(p24-1)-(p24-1)*x^p24;
guards24={TrueQ[FullSimplify[D[g24,{q24,2}]==p24*(p24-1)*delta24^2*(a24^(p24-2)/q24^3+b24^(p24-2)/(1-q24)^3),Assumptions->assume24]],TrueQ[FullSimplify[D[g24,delta24]==p24*(b24^(p24-1)-a24^(p24-1)),Assumptions->assume24]],TrueQ[FullSimplify[D[g24,q24]==f24[a24]-f24[b24],Assumptions->assume24]],TrueQ[FullSimplify[-f24[1/(1-delta24)]==(p24*delta24-1)/(1-delta24)^p24,Assumptions->0<delta24<1&&p24>1]]};
Print["all_p_guards=",guards24];If[!And@@guards24,Exit[1]];
