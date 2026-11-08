(* Explicit inputs.json, output directory, optional worker count. *)
args=Rest[$ScriptCommandLine];If[Length[args]<2,Print["Need inputs.json and output directory"];Exit[2]];
in66=Import[args[[1]],"RawJSON"];out66=args[[2]];
If[!DirectoryQ[out66],CreateDirectory[out66,CreateIntermediateDirectories->True]];
workers66=If[Length[args]>=3,ToExpression[args[[3]]],12];
(* Independent finite Gamma--binomial / beta characteristic-function calculation. *)
ClearAll[mix66];
mix66[spec_]:=Module[{b=spec[[1]],kind=spec[[2]],kap=spec[[3]],wp=spec[[4]],mu,c,z,s,k,t1,t2,om,kk,gg,logs,ww,keep,d,dm,Aj,Ap,Bj,Bp,base,cf,der,lc1,lc2,phase,jj,probs,A,AA,pref,start,xv,h0,h1,ibp,pairs},
start=AbsoluteTime[];mu=cc/.FindRoot[cc==1+Exp[-cc],{cc,N[5/4,wp]},WorkingPrecision->wp];c=mu;z=zz/.FindRoot[Erfc[zz/Sqrt[2]]/2==Sqrt[mu/b]zz,{zz,N[Sqrt[Log[b]],wp]},WorkingPrecision->wp];s=b mu-Sqrt[mu b]z;k=If[kind==1,1,b/2];t1=N[kap,wp]Sqrt[(b-k)/(b k)];t2=-N[kap,wp]Sqrt[k/(b(b-k))];om=t1-t2;
kk=Range[0,Min[b,Ceiling[s/c]-1]];gg=s-c kk;logs=LogGamma[N[b+1,wp]]-LogGamma[N[kk+1,wp]]-LogGamma[N[b-kk+1,wp]]+(b-1)Log[gg];ww=Exp[logs-Max[logs]];ww=ww/Total[ww];keep=Flatten[Position[ww,x_/;x>10^(-wp+20)]];kk=kk[[keep]];gg=gg[[keep]];ww=ww[[keep]];d=-(b-1)kk/gg;dm=ww.d;
If[kind==1,Aj=1-kk/b+kk/b Exp[I om c];Ap=I om kk/b Exp[I om c],
lc1=LogGamma[N[k+1,wp]]-LogGamma[N[Range[0,k]+1,wp]]-LogGamma[N[k-Range[0,k]+1,wp]];lc2=LogGamma[N[b-k+1,wp]]-LogGamma[N[Range[0,b-k]+1,wp]]-LogGamma[N[b-k-Range[0,b-k]+1,wp]];phase=Exp[I om c Range[0,k]];
A=Table[jj=Range[Max[0,K-(b-k)],Min[k,K]];probs=Exp[lc1[[jj+1]]+lc2[[K-jj+1]]-LogGamma[N[b+1,wp]]+LogGamma[N[K+1,wp]]+LogGamma[N[b-K+1,wp]]];{probs.phase[[jj+1]],probs.(I om jj phase[[jj+1]])},{K,kk}];Aj=A[[All,1]];Ap=A[[All,2]]];
If[kind==4,
xv=om gg;h0=Hypergeometric0F1[N[k+1/2,wp],#]&/@(-xv^2/16);h1=Hypergeometric0F1[N[k+3/2,wp],#]&/@(-xv^2/16);Bj=Exp[I xv/2]h0;Bp=-om kk Exp[I xv/2](I h0/2-xv h1/(8(k+1/2))),
(* For Beta(1,b-1), integration by parts has remainder at most
   FallingFactorial[b-1,m]/x^m.  Its derivative remainder is at most
   (1+m/x) times that bound. This converges rapidly for x>b. *)
ibp[xx_,KK_]:=If[xx>6(b-1)/5,
 Module[{term=I(b-1)/xx,su,du,j=1,ep},su=term;du=term;
 While[Abs[term]>10^(-wp-20)&&j<b-1,term=term(b-1-j)/(I xx);su+=term;du+=(j+1)term;j++];If[j==b-1,ep=-Exp[I xx]term;su+=ep;du+=(b-1-I xx)ep];{su,om KK du/xx}],
 {Hypergeometric1F1[N[1,wp],N[b,wp],I xx],-I om KK/b Hypergeometric1F1[N[2,wp],N[b+1,wp],I xx]}];
pairs=MapThread[ibp,{om gg,kk}];Bj=pairs[[All,1]];Bp=pairs[[All,2]]];
pref=Exp[I t2 s];cf=pref ww.(Aj Bj);der=pref ww.((d-dm)Aj Bj+Ap Bj+Aj Bp);
<|"b"->b,"kind"->kind,"kappa"->ToString[kap,InputForm],"precision"->wp,"retained_counts"->Length[kk],"omitted_posterior_mass"->ToString[N[1-Total[ww],20],InputForm],"CF"->ToString[N[cf,35],InputForm],"derivative"->ToString[N[der,35],InputForm],"derivative_abs"->ToString[N[Abs[der],35],InputForm],"seconds"->N[AbsoluteTime[]-start]|>];
specs=({ToExpression[#["b"]],#["kind"],ToExpression[#["kappa"]],#["precision"]}&/@in66["mixture_specs"]);
LaunchKernels[workers66];ParallelEvaluate[$HistoryLength=0];DistributeDefinitions[mix66];start=AbsoluteTime[];rows=ParallelMap[mix66,specs,Method->"FinestGrained"];
Export[FileNameJoin[{out66,"mixture66.json"}],<|"version"->$Version,"workers"->$KernelCount,"seconds"->N[AbsoluteTime[]-start],"rows"->rows,"scope"->"Finite conditional mixture with posterior cutoff; high precision, no directed rounding"|>,"RawJSON"];CloseKernels[];Print["Completed conditional mixture checks"];
