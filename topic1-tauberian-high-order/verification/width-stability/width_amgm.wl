ClearAll[s,kap,v,v0,clo,cup,cylo,cyup,z,lo,xx,up];
Print["Wolfram ",$Version,"; PID ",$ProcessID];
Print["Domain positive constants, c<=C, cy<=Cy, 0<v<=v0, cv<=kap<=Cv, cy<=z<=Cy."];
h=Log[1+kap Exp[s]];
Print["hessian residual: ",FullSimplify[D[h,{s,2}]-kap Exp[s]/(1+kap Exp[s])^2,kap>0&&Element[s,Reals]]];
Print["denominator-cross inequality violation 0<=lo<=xx<=up: ",Reduce[0<=lo<=xx<=up&&xx(1+up)^2<lo(1+xx)^2,{lo,xx,up},Reals]];
Print["sum of positive terms residual: ",Expand[xx(1+up)^2-lo(1+xx)^2-((xx-lo)(1+xx)^2+xx(up-xx)(2+up+xx))]];
Print["n2 gap residual: ",FullSimplify[(1+kap Exp[s])(1+kap Exp[-s])-(1+kap)^2-2kap(Cosh[s]-1),Element[s,Reals]&&kap>0]];
Print["complete"];
