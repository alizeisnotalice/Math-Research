mu=1+ProductLog[1/E];
score[n_,a_,bb_,prec_]:=Module[{c=N[mu,prec],beta,w,r,k,p,weights,avg,q,root,logK},
 beta=N[1+a Sqrt[Log[n]/n],prec];w=N[(1+a^2)Log[n]/2+Log[Log[n]]+bb,prec];
 r=n/beta-1;p=1/(1+Exp[c]);k=Range[0,n];weights=Binomial[n,k]p^k(1-p)^(n-k);
 avg[q_?NumericQ]:=Module[{x=Map[Max[0,#]&,q+w-c k],y=Map[Max[0,#]&,q-c k]},
 weights.(Exp[-beta(q+w-c k)/n](1-beta/n)^(-n) GammaRegularized[n,(1-beta/n)x]+GammaRegularized[n,0,x]-GammaRegularized[n,0,y])-1];
 root=q/.FindRoot[avg[q]==0,{q,N[n mu-Sqrt[n Log[n]],prec]},WorkingPrecision->prec];
 logK=Log[1-Exp[-w]+Exp[-w] n/beta Exp[r] Gamma[n,r]/r^n];root-logK];

n=64;t0=722564147337939/10^13;
out=Table[Module[{cc=N[mu,prec+20],p,kk,wt,lf,f,t1,ys,limit,extra,cost},
p=1/(1+Exp[cc]);kk=Range[0,n];wt=Binomial[n,kk]p^kk(1-p)^(n-kk);
f[t_?NumericQ]:=Sum[If[t>cc k,wt[[k+1]] Exp[-(t-cc k)](t-cc k)^(n-1)/Gamma[n],0],{k,0,n}];
lf=Log[f[N[t0,prec+20]]];
t1=x/.FindRoot[x-t0+Log[f[x]]-lf==Log[n],{x,N[t0+4,prec]},WorkingPrecision->prec];
ys[t_?NumericQ]:=Module[{target=t-t0+Log[f[t]]-lf,y=(t-t0+Log[f[t]]-lf-Log[n])/(n-1)},
Do[y=y-((n-2)Log[1+y]+Log[n+y]-target)/((n-2)/(1+y)+1/(n+y)),{12}];y];
limit=N[n mu+18 Sqrt[n mu],prec];
extra=NIntegrate[(Exp[ys[x]]-1)f[x],{x,t1,limit},WorkingPrecision->prec,AccuracyGoal->22,PrecisionGoal->22];
cost=1-Exp[-(t1-t0)]+NIntegrate[Exp[ys[x]](1+ys[x])/(n+ys[x]) f[x],{x,t1,limit},WorkingPrecision->prec,AccuracyGoal->22,PrecisionGoal->22]/Exp[lf];
<|"precision"->prec,"gain"->ToString[N[t0-Log[cost]-score[n,0,0,prec],25],InputForm],
"mean_residual"->ToString[N[extra-Sum[wt[[k+1]]If[t0>cc k,GammaRegularized[n,0,t0-cc k],0],{k,0,n}],20],InputForm]|>],{prec,{40,60}}];
Print[out];Export[FileNameJoin[{DirectoryName[$InputFileName],"wolfram.json"}],out,"RawJSON"];
