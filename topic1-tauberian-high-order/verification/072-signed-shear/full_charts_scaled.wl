ClearAll[e,f,w,d,lam];
data=Get[FileNameJoin[{DirectoryName[$InputFileName],"full_charts_results.wl"}]];
expr=FullSimplify[(data["total"]/e^4)/.{f->e/d,w->e/lam},e>0&&d>0&&lam>0];
Print["Wolfram ",$Version," PID ",$ProcessID];
Print["exact dimensionless full single-f total: ",Factor[expr]];
Print["expanded delta: ",Collect[Expand[expr],d,Factor]];
Print["leading 7/4 in lambda->0 residual: ",FullSimplify[Limit[Coefficient[expr,d,1],lam->0]-7/4]];
