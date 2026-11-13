out=Flatten[Table[Module[{yr=N[6/(5 Sqrt[n]),80],ys,hr,hs,lr,g,lim},
ys=yr-Log[N[r,80]]/n;hr=Exp[yr];hs=Exp[ys];
lr=(n-2)(Log[1+yr]-Log[1+ys])+Log[n+yr]-Log[n+ys];
g=n(Exp[lr]hr(1+yr)/(n+yr)-hs(1+ys)/(n+ys)-(hr-hs));
lim=r-1-Log[r];
<|"n"->n,"ratio"->ToString[r,InputForm],"normalized_gap"->ToString[N[g,25],InputForm],
"limit"->ToString[N[lim,25],InputForm],"difference"->ToString[N[g-lim,15],InputForm]|>],
{n,{10^6,10^10,10^14}},{r,{1/5,1/2,1,2,5}}]];
Export[FileNameJoin[{DirectoryName[$InputFileName],"wolfram.json"}],out,"RawJSON"];
Print[Take[out,-5]];
