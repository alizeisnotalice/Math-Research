mu=1+ProductLog[1/E];
score[n_,a_,bb_,prec_]:=Module[{c=N[mu,prec],beta,w,r,k,p,weights,avg,q,root,logK},
 beta=N[1+a Sqrt[Log[n]/n],prec];w=N[(1+a^2)Log[n]/2+Log[Log[n]]+bb,prec];
 r=n/beta-1;p=1/(1+Exp[c]);k=Range[0,n];weights=Binomial[n,k]p^k(1-p)^(n-k);
 avg[q_?NumericQ]:=Module[{x=Map[Max[0,#]&,q+w-c k],y=Map[Max[0,#]&,q-c k]},
 weights.(Exp[-beta(q+w-c k)/n](1-beta/n)^(-n) GammaRegularized[n,(1-beta/n)x]+GammaRegularized[n,0,x]-GammaRegularized[n,0,y])-1];
 root=q/.FindRoot[avg[q]==0,{q,N[n mu-Sqrt[n Log[n]],prec]},WorkingPrecision->prec];
 logK=Log[1-Exp[-w]+Exp[-w] n/beta Exp[r] Gamma[n,r]/r^n];root-logK];
rows=Import[FileNameJoin[{DirectoryName[$InputFileName],"results.json"}],"RawJSON"];
out=Flatten[Table[Module[{row=SelectFirst[rows,#["n"]==n&],a,b,g},
 a=Rationalize[row["a"],0];b=Rationalize[row["b"],0];
 g=score[n,a,b,prec]-score[n,0,0,prec];
 <|"n"->n,"precision"->prec,"gain"->ToString[N[g,25],InputForm]|>],{n,{64,256}},{prec,{40,80}}]];
Export[FileNameJoin[{DirectoryName[$InputFileName],"wolfram.json"}],out,"RawJSON"];Print[out];
