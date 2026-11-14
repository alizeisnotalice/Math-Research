base=DirectoryName[$InputFileName];p[t_]:=-5t^19-18t^15+45t^9+30t^4+36;a[t_]:=(t^15+t^10+2t^6)/(t^15+3);
root=t/.FindRoot[p[t]==0,{t,11/10},WorkingPrecision->120,AccuracyGoal->100,PrecisionGoal->100];
root80=t/.FindRoot[p[t]==0,{t,11/10},WorkingPrecision->80,AccuracyGoal->65,PrecisionGoal->65];
mu=N[1+ProductLog[Exp[-1]],100];
lo=Floor[root*10^60]/10^60;hi=lo+10^-60;
res=<|"polynomial"->ToString[p[t],InputForm],"derivative_residue"->FullSimplify[D[a[t],t]-t^5*p[t]/(t^15+3)^2],"sign_changes_descending"->1,"p_at_1"->p[1],"rational_bracket"->(ToString[#,InputForm]&/@{lo,hi}),"bracket_signs_exact"->{p[lo]>0,p[hi]<0},"root"->ToString[N[root,70],InputForm],"precision_doubling_difference"->ToString[Abs[root-root80],InputForm],"grid_exponential_rate"->ToString[N[Log[a[root]],70],InputForm],"optimal_alpha"->ToString[N[root^15/(3+root^15),70],InputForm],"mu_minus_grid_rate"->ToString[N[mu-Log[a[root]],70],InputForm]|>;
Export[FileNameJoin[{base,"grid_rate_wolfram.json"}],res,"RawJSON"];Print[res];Quit[];
