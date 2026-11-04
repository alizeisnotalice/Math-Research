ClearAll[eps,del];
Print["Mathematica: ",$Version];
Print["domain: 0<eps<1/100, 0<=del<=eps/100"];
a=(4eps-del)(eps-2del); b=(eps-3del)(4eps-2del);
Print["residuals: ",{Expand[a-2eps^2],Expand[b-2eps^2]}];
Print["either strict violation: ",Reduce[0<eps<1/100&&0<=del<=eps/100&&(a<2eps^2||b<2eps^2),{eps,del},Reals]];
Print["equality or worse: ",Reduce[0<eps<1/100&&0<=del<=eps/100&&(a<=2eps^2||b<=2eps^2),{eps,del},Reals]];
