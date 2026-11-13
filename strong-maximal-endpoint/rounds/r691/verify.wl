mu=1+ProductLog[1/E];
score[n_,a_,old_,prec_]:=Module[{c=N[mu,prec],b,w,r,k,p,weights,x,avg,q,res,logK},
 b=N[1+a Sqrt[Log[n]/n],prec]; w=N[(1+a^2)Log[n]/2+If[old,1,2]Log[Log[n]],prec];
 r=n/b-1;p=1/(1+Exp[c]);k=Range[0,n];weights=Binomial[n,k]p^k(1-p)^(n-k);
 avg[q_?NumericQ]:=Module[{x=Map[Max[0,#]&,q+w-c k],y=Map[Max[0,#]&,q-c k]},
 weights.(Exp[-b(q+w-c k)/n](1-b/n)^(-n) GammaRegularized[n,(1-b/n)x]+GammaRegularized[n,0,x]-GammaRegularized[n,0,y])-1];
 res=q/.FindRoot[avg[q]==0,{q,N[n mu-Sqrt[n Log[n]],prec]},WorkingPrecision->prec];
 logK=Log[1-Exp[-w]+Exp[-w] n/b Exp[r] Gamma[n,r]/r^n];
 res-logK];
Do[Print[{n,prec,N[score[n,Sqrt[mu],False,prec]-score[n,0,True,prec],22]}],{n,{16,64}},{prec,{40,80}}];
