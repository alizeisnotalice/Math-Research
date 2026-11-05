root=DirectoryName[DirectoryName[$InputFileName]];
data=Import[FileNameJoin[{root,"rounds","r002","multilayer.json"}],"RawJSON"];
exact[x_]:=Rationalize[x,0];
tail[n_,r_,q_]:=Module[{p=r/(1+r),x},
 Sum[x=q+k Log[r]; Binomial[n,k]p^k(1-p)^(n-k)If[x<=0,1,GammaRegularized[n,x]],{k,0,n}]];
checks=Map[Function[row,Module[{n=row["n"],r,t,h,a,p,w,e,j,score},
 r=N[exact[row["r"]],80]; t=N[exact/@row["thresholds"],80];
 h=N[exact/@row["heights"],80]; a=tail[n,r,#]&/@t;
 p=a-Append[Rest[a],0]; w=Exp[-t]-Append[Rest[Exp[-t]],0];
 e=p.h;
 j=Sum[w[[k]]If[h[[k]]==0,0,h[[k]](1+Max[0,Log[h[[k]]]])^(n-1)],{k,Length[h]}];
 score=-Log[j];
 <|"n"->n,"expectation80"->ToString[N[e,35],InputForm],
 "strict_margin80"->ToString[N[e-1,25],InputForm],
 "log_score80"->ToString[N[score,35],InputForm],
 "matlab_score_difference"->N[Abs[score-exact[row["log_lower_candidate"]]]]|>
]],data["rows"]];
result=<|"wolfram_version"->$Version,"status"->"80-digit numerical re-evaluation; no directed certificate", "checks"->checks|>;
Export[FileNameJoin[{root,"rounds","r002","crosscheck.json"}],result,"RawJSON"];
Print[ExportString[result,"RawJSON"]];
