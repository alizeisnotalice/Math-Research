(* Run: wolframscript -file verify.wl PATH_TO_CERTIFICATES_JSON OUTPUT_TXT *)
args=Take[$ScriptCommandLine,-2]; input=args[[1]]; output=args[[2]];
Print["VERSION ",$Version," MASTER ",$ProcessID];
Print["implicit leading equation: ",InputForm[Limit[4(1+2u)^(n-2) Sum[Binomial[n,j](1/(1+2u))^(n-j) H^j u^(6j-6),{j,1,n}]-k^2,u->0,Assumptions->Element[n,Integers]&&n>=2]]];
ass=Element[n,Integers]&&n>=2&&k>0; hlead=k^2/(4n);aLead=2(k hlead-(2n/k)hlead^2);gammaLead=FullSimplify[(2^n k^2-2aLead 2^(n-2))/2^(n+8),ass];target=(k^2-k^3/(8n))/256;
Print["symbolic ",InputForm[{"areaLead"->FullSimplify[aLead,ass],"gammaIdentity"->FullSimplify[gammaLead==target,ass],"criticalK"->FullSimplify[Solve[D[target,k]==0&&k>0,k,Reals],ass],"criticalValue"->FullSimplify[target/.k->16n/3,ass],"edgeExtraPower"->3}]];
primitive=((r+h)^(n+1)-r^(n+1))/(n+1)-r^n h;
Print["areaPrimitiveIdentity ",InputForm[FullSimplify[{D[primitive,h]==(r+h)^n-r^n,(primitive/.h->0)==0},Element[n,Integers]&&n>=2&&r>0&&h>0]]];
data=Import[input,"RawJSON"];records=data["cases"];Print["INPUT_COUNT ",Length[records]];
Clear[decExact];decExact[s_String]:=Module[{parts=StringSplit[s,RegularExpression["[eE]"]],base,ex=0,sgn=1,dp,nums},base=First[parts];If[Length[parts]>1,ex=ToExpression[parts[[2]]]];If[StringStartsQ[base,"-"],sgn=-1;base=StringDrop[base,1]];dp=If[StringContainsQ[base,"."],StringLength[Last[StringSplit[base,"."]]],0];nums=StringReplace[base,"."->""];sgn FromDigits[nums]10^(ex-dp)];
Clear[calc];calc[rec_]:=Module[{nn=rec["n"],kap=ToExpression[rec["kappa"]],u=ToExpression[rec["u"]],r,del,e,K,m,det,h0,hs,h,A,b,g,lo,hi,res,midErr,inside,sign},
 r=1/(1+2u);del=u/(1+2u);e=kap u^4;K=4u^2 r^(2-nn);m=(2del)^nn;det=1-e^2;h0=kap^2 u^6/(4nn);
 h=hs/.FindRoot[K Sum[Binomial[nn,j]r^(nn-j)hs^j,{j,1,nn}]==e^2,{hs,N[h0,150]},WorkingPrecision->150,AccuracyGoal->125,PrecisionGoal->125];
 A=2h NIntegrate[e-(K/e)Sum[Binomial[nn,j]r^(nn-j)(h zz)^j,{j,1,nn}],{zz,0,1},WorkingPrecision->120,AccuracyGoal->85,PrecisionGoal->85];
 b=2A(2del)^(nn-2)+If[nn>=3,4(nn-2)h(e+h)^2(2del)^(nn-3),0];g=(m e^2-b(1-m))/((det-m)(1-m)(2u)^nn(2u)^8);
 lo=decExact[rec["G"][[1]]];hi=decExact[rec["G"][[2]]];inside=TrueQ[lo<=g<=hi];sign=Which[g>0,"positive",g<0,"negative",True,"unresolved"];
 res=Abs[K Sum[Binomial[nn,j]r^(nn-j)h^j,{j,1,nn}]-e^2];midErr=Abs[g-(lo+hi)/2];
 {nn,rec["kappa"],rec["u"],inside,sign==rec["G_sign"],N[midErr,8],N[res,8],$ProcessID}];
LaunchKernels[15];Print["KERNEL_COUNT ",Length[Kernels[]]];pids=ParallelEvaluate[$ProcessID];DistributeDefinitions[calc,decExact];
results=ParallelMap[calc,records,Method->"FinestGrained"];
summary={"version"->$Version,"master"->$ProcessID,"kernelCount"->Length[pids],"kernelPIDs"->pids,"cases"->Length[results],"intervalPass"->Count[results[[All,4]],True],"signPass"->Count[results[[All,5]],True],"maxMidpointError"->Max[results[[All,6]]],"maxRootResidual"->Max[results[[All,7]]],"numericNotIntervalProof"->True};
Export[output,ToString[InputForm[{summary,results}]],"Text"];Print["SUMMARY ",InputForm[summary]];CloseKernels[];
