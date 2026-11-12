Clear[w,r];m=1-w/(Exp[w]-1);v=1-w^2 Exp[w]/(Exp[w]-1)^2;
checks={Limit[m,w->0]==0,Limit[v,w->0]==0,Limit[m,w->Infinity]==1,Limit[v,w->Infinity]==1,FullSimplify[D[1-r Log[r]/(1+r),r]==(-1-r-Log[r])/(1+r)^2,Assumptions->r>0]};
mu=1+ProductLog[1/E];critical=2 Sqrt[mu] Sqrt[2] InverseErf[1/2];
Print[{checks,N[critical,24],Abs[N[critical,48]-N[critical,24]]<10^-22}];If[!And@@checks,Exit[1]];
