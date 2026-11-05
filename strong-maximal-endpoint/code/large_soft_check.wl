root=DirectoryName[DirectoryName[$InputFileName]];
data=Import[FileNameJoin[{root,"rounds","r008","large_soft.json"}],"RawJSON"];
row=First[data["rows"]]; n=row["n"]; r=N[Rationalize[data["r"],0],60];
q=N[Rationalize[row["q"],0],60]; theta=1/n;
p=r^(1-theta)/(1+r^(1-theta)); mode=Floor[(n+1)p]; sd=Sqrt[n p(1-p)];
lower=Max[0,Floor[mode-12sd]]; upper=Min[n,Ceiling[mode+12sd]];
tail=Sum[x=(1-theta)(q+k Log[r]);
 Binomial[n,k]p^k(1-p)^(n-k)If[x<=0,1,GammaRegularized[n,x]],{k,lower,upper}];
omitted=1-Sum[Binomial[n,k]p^k(1-p)^(n-k),{k,lower,upper}];
mgf=n(-Log[1-theta]+Log[1+r^(1-theta)]-Log[1+r]);
integral=Sum[Factorial[n-1]/Factorial[n-1-k](1/n)^k/(1-1/n)^(k+1),{k,0,n-1}];
score=q-Log[integral];
result=<|"wolfram_version"->$Version,"n"->n,
 "truncated_expectation_margin60"->ToString[N[Exp[mgf-theta q]tail-1,25],InputForm],
 "omitted_binomial_mass60"->ToString[N[omitted,20],InputForm],
 "matlab_log_score_difference"->N[Abs[score-Rationalize[row["log_lower_candidate"],0]]],
 "status"->"60-digit diagnostic at n=10000; not directed certificate"|>;
Export[FileNameJoin[{root,"rounds","r008","crosscheck.json"}],result,"RawJSON"];
Print[ExportString[result,"RawJSON"]];
