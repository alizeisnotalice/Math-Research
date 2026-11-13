out=Table[Module[{x=N[n(1+1/128),80],v,exact,stable},
v=x/n-1;exact=(n-1)Log[x]-x-LogGamma[n];
stable=-Log[2 Pi n]/2-Log[1+v]+n(Log[1+v]-v)-1/(12n)+1/(360n^3)-1/(1260n^5);
<|"n"->n,"formula_difference"->ToString[N[stable-exact,25],InputForm],
"stirling_remainder_bound"->ToString[N[1/(1680n^7),25],InputForm]|>],
{n,{262144,1048576,4194304}}];
Print[out];Export[FileNameJoin[{DirectoryName[$InputFileName],"wolfram.json"}],out,"RawJSON"];
