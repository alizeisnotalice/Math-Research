(* Explicit inputs.json, output directory, optional worker count. *)
args=Rest[$ScriptCommandLine];If[Length[args]<2,Print["Need inputs.json and output directory"];Exit[2]];
in66=Import[args[[1]],"RawJSON"];out66=args[[2]];
If[!DirectoryQ[out66],CreateDirectory[out66,CreateIntermediateDirectories->True]];
workers66=If[Length[args]>=3,ToExpression[args[[3]]],12];
ClearAll[conditional66];
conditional66[spec_]:=Module[{b=spec[[1]],kind=spec[[2]],kap0=spec[[3]],ref=spec[[4]],wp,mu,c,z,s,theta,a,p,v,m,k,t1,t2,h,cut,js,xx,lg,sc,part,out,I0,I1,D0,D1,result,start},
start=AbsoluteTime[];wp=2 Ceiling[Log[10,b]]+80+If[ref==1,40,0];mu=c0/.FindRoot[c0==1+Exp[-c0],{c0,N[5/4,wp]},WorkingPrecision->wp];c=mu;
z=z0/.FindRoot[Erfc[z0/Sqrt[2]]/2==Sqrt[mu/b] z0,{z0,N[Sqrt[Log[b]],wp]},WorkingPrecision->wp];s=b mu-Sqrt[mu b] z;
theta=th/.FindRoot[1/(1-th)+c/(1+Exp[c(1-th)])==s/b,{th,-N[z/Sqrt[mu b],wp]},WorkingPrecision->wp];a=1-theta;p=1/(1+Exp[c a]);m=1/a+c p;v=1/a^2+c^2 p(1-p);
k=If[kind==1,1,b/2];t1=N[kap0,wp]Sqrt[N[(b-k)/(b k),wp]];t2=-N[kap0,wp]Sqrt[N[k/(b(b-k)),wp]];h=If[ref==1,1/16,1/8];cut=If[ref==1,40,32];js=N[Range[-cut/h,cut/h]h,wp];
lg[x_]:=-Log[1-I x/a]+Log[1+Exp[-c a+I c x]]-Log[1+Exp[-c a]]-I m x;
sc[x_]:=a p+Exp[-c a+I c x](-a+I x)/(1+Exp[-c a+I c x]);
part[u_]:=Module[{tau=u/Sqrt[b v],base,num,d0,d1},base=Exp[b lg[tau]];num=Exp[k lg[t1+tau]+(b-k)lg[t2+tau]];{base,num,base b sc[tau],num(k sc[t1+tau]+(b-k)sc[t2+tau])}];
out=Total[part/@js]h;{I0,I1,D0,D1}=out;result=(D1 I0-I1 D0)/I0^2;
<|"b"->ToString[b],"kind"->kind,"kappa"->ToString[kap0,InputForm],"refined"->ref,"precision"->wp,"step"->ToString[h,InputForm],"cutoff"->cut,"theta"->ToString[N[theta,30],InputForm],"CF"->ToString[N[I1/I0,30],InputForm],"derivative"->ToString[N[result,30],InputForm],"derivative_abs"->ToString[N[Abs[result],30],InputForm],"normalizer_imaginary_residual"->ToString[N[Abs[Im[I0]],20],InputForm],"normalizer"->ToString[N[Re[I0],30],InputForm],"tilt_mean_residual"->ToString[N[Abs[m-s/b],20],InputForm],"seconds"->N[AbsoluteTime[]-start]|>];
specs=({ToExpression[#["b"]],#["kind"],ToExpression[#["kappa"]],#["refined"]}&/@in66["fourier_specs"]);
LaunchKernels[workers66];ParallelEvaluate[$HistoryLength=0];DistributeDefinitions[conditional66];start=AbsoluteTime[];rows=ParallelMap[conditional66,specs,Method->"FinestGrained"];
Export[FileNameJoin[{out66,"conditional66.json"}],<|"version"->$Version,"workers"->$KernelCount,"seconds"->N[AbsoluteTime[]-start],"rows"->rows,"scope"->"Conditional Fourier ratios; finite patterns and quadrature diagnostics, not a uniform operator theorem"|>,"RawJSON"];CloseKernels[];Print["Completed conditional Fourier checks"];
