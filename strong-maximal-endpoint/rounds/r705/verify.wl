pc=1/(1+Exp[c]);
out=<|"bernoulli_exponential_moment"->ToString[FullSimplify[1-pc+pc Exp[c]],InputForm],
"moment_maxima"->Table[<|"order"->j,"critical_point"->2j,"maximum"->ToString[(2j/E)^j,InputForm],
"stationary_check"->ToString[FullSimplify[D[x^j Exp[-x/2],x]/.x->2j],InputForm]|>,{j,1,8}],
"zero_ratio_increment_mean"->1,"infinite_c_increment_mean"->1|>;
Print[out];Export[FileNameJoin[{DirectoryName[$InputFileName],"wolfram.json"}],out,"RawJSON"];
