v=-Log[m]/2+x^2(1-1/m)/2;
q=Exp[-x^2/2]/Sqrt[2 Pi];p=Exp[-x^2/(2m)]/Sqrt[2 Pi m];
out=<|"gaussian_budget"->ToString[FullSimplify[Integrate[q Exp[v],{x,-Infinity,Infinity},Assumptions->m>0]],InputForm],
"gaussian_gain"->ToString[FullSimplify[Integrate[p v,{x,-Infinity,Infinity},Assumptions->m>0]],InputForm],
"modular_ratio_limit"->ToString[Limit[Exp[u/n](1+u/(n(1+d/Sqrt[n])))^(n-1),n->Infinity,Assumptions->{Element[{u,d},Reals],d>0}],InputForm]|>;
Print[out];Export[FileNameJoin[{DirectoryName[$InputFileName],"wolfram.json"}],out,"RawJSON"];
