(* Usage: wolframscript -file actions65.wl inputs.json output-directory precision [workers] *)
args=Rest[$ScriptCommandLine];
If[Length[args]<2,Print["Need inputs.json and output directory"];Exit[2]];
inputPath=args[[1]];outDir=args[[2]];precision65=If[Length[args]>=3,ToExpression[args[[3]]],80];
If[!IntegerQ[precision65]||precision65<40,Exit[2]];
If[!DirectoryQ[outDir],CreateDirectory[outDir,CreateIntermediateDirectories->True]];
ClearAll[case65];
case65[r_]:=Module[{wp=precision65,n,nn,c,q,al,p,k,w,lg,g,g2,F,Q,Q2,pos,pos2,hm1,pdf,pdf2,dp,aa,ac,acc,grad,curv,start},
start=AbsoluteTime[];n=Round[r["N"]];nn=N[n,wp];c=ToExpression[r["c_decimal"]<>"`"<>ToString[wp]];q=ToExpression[r["q_decimal"]<>"`"<>ToString[wp]];al=ToExpression[r["slope_decimal"]<>"`"<>ToString[wp]];p=1/(1+Exp[c]);k=Range[Max[0,Floor[n p-Sqrt[110n]]],Min[n,Ceiling[n p+Sqrt[110n]]]];lg=LogGamma[nn];w=Exp[LogGamma[nn+1]-LogGamma[N[k+1,wp]]-LogGamma[N[n-k+1,wp]]-c k-n Log[1+Exp[-c]]];g=q-c k;g2=g+n/al;
(* Gamma(n+1,x)=Gamma(n,x)+x^n exp(-x)/Gamma(n+1).
   Two tail calls suffice; the earlier implementation used six. *)
Q=GammaRegularized[nn,#]&/@g;Q2=GammaRegularized[nn,#]&/@g2;F=1-Q;
pdf=Exp[(n-1)Log[g]-g-lg];pdf2=Exp[(n-1)Log[g2]-g2-lg];
pos=(n-g)Q+g pdf;pos2=(n-g2)Q2+g2 pdf2;hm1=-F+al/n (pos-pos2);
dp=(n-1)/g-1;aa=n p-k;ac=k pdf+al k/n (Q-Q2);acc=-k^2 pdf dp+al k^2/n (pdf-pdf2);grad=w.(aa hm1+ac);curv=w.((aa^2-n p(1-p))hm1+2aa ac+acc);
<|"N"->n,"profile"->r["profile"],"q_decimal"->r["q_decimal"],"c_decimal"->r["c_decimal"],"slope_decimal"->r["slope_decimal"],"precision"->wp,"mean_minus_one"->ToString[N[w.hm1,35],InputForm],"gradient"->ToString[N[grad,35],InputForm],"curvature"->ToString[N[curv,35],InputForm],"seconds"->N[AbsoluteTime[]-start]|>];
inputs=Import[inputPath,"RawJSON"];LaunchKernels[If[Length[args]>=4,ToExpression[args[[4]]],12]];ParallelEvaluate[$HistoryLength=0];DistributeDefinitions[case65,precision65,outDir];start=AbsoluteTime[];rows=ParallelMap[Function[r,Module[{x=case65[r]},Export[outDir<>"/wl_"<>ToString[x["N"]]<>"_"<>ToString[x["profile"]]<>".json",x,"RawJSON"];x]],inputs,Method->"FinestGrained"];
Export[outDir<>"/actions65_wolfram.json",<|"status"->"complete","version"->$Version,"workers"->$KernelCount,"seconds"->N[AbsoluteTime[]-start],"rows"->rows,"scope"->"High-precision response at fixed decimal profiles/actions, not certified global optima"|>,"RawJSON"];CloseKernels[];Print["Completed actions65 high-precision checks"];
