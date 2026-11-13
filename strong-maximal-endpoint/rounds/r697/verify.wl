mu=1+ProductLog[1/E];data={{64,2871751305155139/2500000000000000,3812474239382891/50000000000000},{65536,126280439373593/125000000000000,8312739222316067/100000000000},{4194304,1001625583303663/1000000000000000,1338868075066401/250000000}};
out=Table[Module[{n=row[[1]],beta=N[row[[2]],80],t1=N[row[[3]],80],rate,z,q,m,v},
rate=1-beta/n;z=rate n/beta;q=GammaRegularized[n,z];
m=n/rate GammaRegularized[n+1,z]/q;
v=n(n+1)/rate^2 GammaRegularized[n+2,z]/q-m^2;
<|"n"->n,"budget_center"->ToString[N[(t1-n/beta+m-n mu)/Sqrt[n],30],InputForm],
"budget_variance"->ToString[N[v/n,30],InputForm]|>],{row,data}];
Print[out];Export[FileNameJoin[{DirectoryName[$InputFileName],"wolfram.json"}],out,"RawJSON"];
