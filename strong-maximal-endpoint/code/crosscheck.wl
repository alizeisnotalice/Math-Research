(* Independent arbitrary-precision Gamma/binomial evaluation at MATLAB candidates. *)
root = DirectoryName[DirectoryName[$InputFileName]];
data = Import[FileNameJoin[{root,"rounds","r001","baseline.json"}],"RawJSON"];
tail[n_, r_, q_, prec_] := Module[{rr=N[Rationalize[r,0],prec],qq=N[Rationalize[q,0],prec],p,x},
 p=rr/(1+rr);
 Sum[x=qq+k Log[rr]; Binomial[n,k] p^k (1-p)^(n-k) If[x<=0,1,GammaRegularized[n,x]],{k,0,n}]
];
checks = Map[Function[row, Module[{n=row["n"],a40,a80,score},
 a40=tail[n,row["r"],row["q"],40]; a80=tail[n,row["r"],row["q"],80];
 score=Log[a80]+Rationalize[row["q"],0]-(n-1) Log[1-Log[a80]];
 <|"n"->n,"alpha80"->ToString[N[a80,35],InputForm],
 "precision40_80_difference"->ToString[N[Abs[a40-a80],15],InputForm],
 "matlab_alpha_difference"->N[Abs[a80-Rationalize[row["alpha"],0]]],
 "log_score80"->ToString[N[score,35],InputForm],
 "matlab_score_difference"->N[Abs[score-Rationalize[row["log_lower_candidate"],0]]]|>
 ]],data["rows"]];
result=<|"wolfram_version"->$Version,"input_policy"->"Exact binary MATLAB candidate coordinates via Rationalize[x,0]; recomputed at 40 and 80 digits", "checks"->checks|>;
Export[FileNameJoin[{root,"rounds","r001","crosscheck.json"}],result,"RawJSON"];
Print[ExportString[result,"RawJSON"]];
