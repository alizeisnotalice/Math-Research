mu=1+ProductLog[1/E];data={{64,722564147337939/10000000000000,3487320373735237/5000000000000000},{256,7584309221102491/25000000000000,327386929368049/625000000000000},{1024,782416052813487/625000000000,1382327117667989/4000000000000000},{4096,2551801248715943/500000000000,2044102285903243/10000000000000000},{16384,1032253566213799/50000000000,224552970388227/2000000000000000},{65536,8311573032712653/100000000000,737710782594071/12500000000000000}};
out=Table[Module[{n=row[[1]],t0=N[row[[2]],80],th=N[row[[3]],80],c=N[mu,80],p,j,x,lf,dd,t,mg,b},
p=1/(1+Exp[c]);j=Floor[n p];x=t0-c j;
lf=LogGamma[n+1]-LogGamma[j+1]-LogGamma[n-j+1]+j Log[p]+(n-j)Log[1-p]-x+(n-1)Log[x]-LogGamma[n];
dd=Exp[(-t0-lf+n Log[2(1-p)]-LogGamma[n])/(n-1)];
t=n c+18 Sqrt[n c];mg=-n Log[1-th]+n Log[1-p+p Exp[c th]];
b=mg-(th-dd)t;
<|"n"->n,"log_mean"->ToString[N[b,30],InputForm],"log_scaled_budget"->ToString[N[b-lf,30],InputForm]|>],{row,data}];
Print[out];Export[FileNameJoin[{DirectoryName[$InputFileName],"wolfram.json"}],out,"RawJSON"];
