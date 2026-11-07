ClearAll[c,d,u,s,a,x,alpha,beta,T,M,A,p,m];alpha=Sqrt[2]-1;beta=2-Sqrt[2];p[u_]:=1/(1+Exp[u]);m[u_]:=1+u p[u];T=1+alpha Exp[-c]+beta Exp[-d];M=1+(alpha c Exp[-c]+beta d Exp[-d])/T;A=T/Sqrt[(1+Exp[-c])(1+Exp[-d])];
checks=<|"hellinger_piecewise_integral"->TrueQ[FullSimplify[(Integrate[Exp[-x],{x,0,c}]+Sqrt[2]Integrate[Exp[-x],{x,c,d}]+2Integrate[Exp[-x],{x,d,Infinity}])-T,Assumptions->0<c<d]==0],"hellinger_diagonal_slope"->TrueQ[FullSimplify[(D[Log[A],d]/.d->c)+(3/2-Sqrt[2])p[c]]==0],"geometric_mean_partial_c"->TrueQ[FullSimplify[(D[M,c]/.{c->u,d->u})-alpha D[m[u],u]]==0],"geometric_mean_partial_d"->TrueQ[FullSimplify[(D[M,d]/.{c->u,d->u})-beta D[m[u],u]]==0]|>;
b2=2(a/(2(s-d))+Integrate[1/(s-x),{x,a,d},Assumptions->0<a<d<s/2]+(s/2-d)/(s-d));checks["exact_b2_envelope"]=TrueQ[FullSimplify[b2-(1-(d-a)/(s-d)+2Log[(s-a)/(s-d)]),Assumptions->0<a<d<s/2]==0];
f2=4Exp[-s](s-c)/(1+Exp[-c])^2;checks["b2_flux_identity"]=TrueQ[FullSimplify[2p[c]-D[Log[f2],c]-1/(s-c),Assumptions->0<c<s/2]==0];
checks["affine_Dirichlet_Jacobian"]=And@@Table[TrueQ[Det[(V/(V+s))IdentityMatrix[k]]==(V/(V+s))^k],{k,1,5}];
joint00=s/(4(s-d));joint10=Log[(s-a)/(s-d)]/4+1/4;
joint11=Log[(s-a)/(s-d)]/2+(s/2-d)/(2(s-d));
checks["exact_b2_joint_label_envelope"]=TrueQ[FullSimplify[joint00+2joint10+joint11-(1+Log[(s-a)/(s-d)]),Assumptions->0<a<d<s/2]==0];
out=<|"version"->$Version,"checks"->checks,"all_pass"->And@@Values[checks],"scope"->"Exact symbolic identities; classifier/asymptotic proof separately written"|>;Export[If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],FileNameJoin[{$TemporaryDirectory,"verify64.json"}]],out,"RawJSON"];Print[out];If[!out["all_pass"],Exit[1]];
