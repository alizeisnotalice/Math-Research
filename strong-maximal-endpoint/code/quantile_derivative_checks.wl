Clear[c,l,n];
m[c_,l_]:=(1+Exp[-(1+l)c])/((1+l)(1+Exp[-c]));
p[c_]:=1/(1+Exp[c]);
res=FullSimplify[D[m[c,l],c]-p[c](m[c,l]-Exp[-l c]),Assumptions->{c>0,l>0}];
mu=1+ProductLog[Exp[-1]];sigma=Sqrt[mu];
g[n_,y_]:=If[y<=0,0,GammaRegularized[n,0,y]];
f[n_,y_]:=If[y<=0,0,Exp[(n-1)Log[y]-y-LogGamma[n]]];
rows=Table[Module[{cc,pp,q,weights,fn,previous,direct,identity},
 cc=N[mu,80]; pp=p[cc];q=Floor[N[nn mu-sigma Sqrt[nn Log[nn]],40]];
 weights=Table[Binomial[nn,k]pp^k(1-pp)^(nn-k),{k,0,nn}];
 fn=Sum[weights[[k+1]]g[nn,q-cc k],{k,0,nn}];
 previous=Sum[Binomial[nn-1,k]pp^k(1-pp)^(nn-1-k)g[nn-1,q-cc-cc k],{k,0,nn-1}];
 direct=Sum[weights[[k+1]]((nn pp-k)g[nn,q-cc k]-k f[nn,q-cc k]),{k,0,nn}];
 identity=nn pp(fn-previous);
 <|"n"->nn,"q"->q,"CDF_derivative_residual"->ToString[N[direct-identity,20],InputForm],"identity_agrees_50_digits"->TrueQ[Abs[direct-identity]<10^-50]|>],{nn,{8,32,64}}];
result=<|"wolfram_version"->$Version,"Laplace_derivative_symbolic_residual"->ToString[res,InputForm],"working_precision"->80,"rows"->rows|>;
root=DirectoryName[DirectoryName[$InputFileName]];
Export[FileNameJoin[{root,"rounds","r016","symbolic_checks.json"}],result,"RawJSON"];
Print[ExportString[result,"RawJSON"]];
