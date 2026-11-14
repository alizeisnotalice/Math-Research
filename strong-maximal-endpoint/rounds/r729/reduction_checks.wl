Print[FullSimplify[D[(1+b)/(x+b),b]-(x-1)/(x+b)^2]];
Print[Integrate[Exp[-aa Cos[t]]aa Sin[t],{t,0,2Pi},Assumptions->Element[aa,Reals]]];
Print[TimeConstrained[Assuming[n>0&&Element[qq,Reals],FullSimplify[Limit[Gamma[n,t-qq]/Gamma[n,t],t->Infinity]]],15,"not-completed"]];
