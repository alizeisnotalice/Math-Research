Clear[v,r,s,x];checks={FullSimplify[Integrate[v^(-1/s),{v,0,r},Assumptions->s>1&&r>0]==r^(1-1/s)/(1-1/s),Assumptions->s>1&&r>0],FullSimplify[D[-r Log[r]/(1+r),r]==(-1-r-Log[r])/(1+r)^2,Assumptions->r>0],FullSimplify[(D[s (Log[1+x]-x),x]/.x->0)==0]};
Print[checks];If[!And@@checks,Exit[1]];
