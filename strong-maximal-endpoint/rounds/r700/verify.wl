data={{64,143587540956379/125000000000000},{4194304,400650229594289/400000000000000}};out=Table[Module[{n=row[[1]],beta=N[row[[2]],80],z,rate,tau,egx,der},
z=n/beta-1;rate=1-beta/n;tau=Exp[n Log[z]-z-LogGamma[n]-Log[GammaRegularized[n,z]]];
egx=n/rate-n/beta+tau/rate;
der=beta (D[Log[n/b]+n/b-1+Log[Gamma[n,n/b-1]]-n Log[n/b-1],b]/.b->beta);
<|"n"->n,"derivative_identity_error"->ToString[N[der-(egx-1),30],InputForm]|>],{row,data}];
Print[out];Export[FileNameJoin[{DirectoryName[$InputFileName],"wolfram.json"}],out,"RawJSON"];
