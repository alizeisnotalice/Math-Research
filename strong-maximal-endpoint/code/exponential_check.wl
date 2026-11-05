root=DirectoryName[DirectoryName[$InputFileName]];
data=Import[FileNameJoin[{root,"rounds","r003","exponential.json"}],"RawJSON"];
checks=Map[Function[row,Module[{n=row["n"],r,a,q,theta,p,mgf,tail,x,int,score},
 r=N[Rationalize[row["r"],0],60]; a=N[Rationalize[row["a"],0],60];
 q=N[Rationalize[row["q"],0],60]; theta=a/n;
 p=r^(1-theta)/(1+r^(1-theta));
 mgf=n(-Log[1-theta]+Log[1+r^(1-theta)]-Log[1+r]);
 tail=Sum[x=(1-theta)(q+k Log[r]);
 Binomial[n,k]p^k(1-p)^(n-k)If[x<=0,1,GammaRegularized[n,x]],{k,0,n}];
 int=Sum[Factorial[n-1]/Factorial[n-1-k](a/n)^k/(1-a/n)^(k+1),{k,0,n-1}];
 score=q-Log[int];
 <|"n"->n,"expectation_margin60"->ToString[N[Exp[mgf-theta q]tail-1,25],InputForm],
 "log_score60"->ToString[N[score,30],InputForm],
 "matlab_score_difference"->N[Abs[score-Rationalize[row["log_lower_candidate"],0]]]|>
]],Select[data["rows"],MemberQ[{16,256,1024,4096},#["n"]]&]];
result=<|"wolfram_version"->$Version,"status"->"60-digit numerical checks, not directed certificate","checks"->checks|>;
Export[FileNameJoin[{root,"rounds","r003","crosscheck.json"}],result,"RawJSON"];
Print[ExportString[result,"RawJSON"]];
