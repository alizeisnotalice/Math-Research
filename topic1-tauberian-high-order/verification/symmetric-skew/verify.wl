(* wolframscript -file symmetric_verify.wl CERTIFICATES_JSON OUTPUT_TXT *)
args=Take[$ScriptCommandLine,-2];input=args[[1]];output=args[[2]];
Print["VERSION ",$Version," MASTER ",$ProcessID];data=Import[input,"RawJSON"];records=data["cases"];Print["INPUT_COUNT ",Length[records]];
Clear[decExact];decExact[s_String]:=Module[{parts=StringSplit[s,RegularExpression["[eE]"]],base,ex=0,sgn=1,dp,nums},base=First[parts];If[Length[parts]>1,ex=ToExpression[parts[[2]]]];If[StringStartsQ[base,"-"],sgn=-1;base=StringDrop[base,1]];dp=If[StringContainsQ[base,"."],StringLength[Last[StringSplit[base,"."]]],0];nums=StringReplace[base,"."->""];sgn FromDigits[nums]10^(ex-dp)];
Clear[calc];calc[rec_]:=Module[{nn=rec["n"],kap=ToExpression[rec["kappa"]],u=ToExpression[rec["u"]],r,del,e,K,m,det,h0,hs,hms,h,hm,A,bpair,btotal,g,lo,hi,res1,resM,midErr,inside,sign,M,staridx,Midx,starOK,mOK},
 r=1/(1+2u);del=u/(1+2u);e=kap u^4;K=4u^2 r^(2-nn);m=(2del)^nn;det=(1-e)^(nn-1)(1+(nn-1)e);M=3^nn;h0=kap^2 u^6/(4nn);
 h=hs/.FindRoot[K Sum[Binomial[nn,j]r^(nn-j)hs^j,{j,1,nn}]==e^2,{hs,N[h0,180]},WorkingPrecision->180,AccuracyGoal->150,PrecisionGoal->150];
 hm=hms/.FindRoot[K Sum[Binomial[nn,j]r^(nn-j)hms^j,{j,1,nn}]==M e^2,{hms,N[M h0,180]},WorkingPrecision->180,AccuracyGoal->150,PrecisionGoal->150];
 A=2h NIntegrate[e-(K/e)Sum[Binomial[nn,j]r^(nn-j)(h zz)^j,{j,1,nn}],{zz,0,1},WorkingPrecision->150,AccuracyGoal->110,PrecisionGoal->110];
 bpair=2A(2del)^(nn-2);
 btotal=Binomial[nn,2] bpair+If[nn>=3,8 Binomial[nn,3](e+hm)^3(2del)^(nn-3),0];
 (* det stays exact up to here; in particular 1-det is exact rational. *)
 g=(m(1-det)-btotal(1-m))/((det-m)(1-m)(2u)^nn(2u)^8);lo=decExact[rec["G"][[1]]];hi=decExact[rec["G"][[2]]];inside=TrueQ[lo<=g<=hi];sign=Which[g>0,"positive",g<0,"negative",True,"unresolved"];
 staridx=FromDigits[StringDrop[rec["roots"]["hstar"]["lower_numerator_hex"],2],16];Midx=FromDigits[StringDrop[rec["roots"]["hM"]["lower_numerator_hex"],2],16];
 starOK=TrueQ[(del-e)staridx/2^400<=h<=(del-e)(staridx+1)/2^400];mOK=TrueQ[(del-e)Midx/2^400<=hm<=(del-e)(Midx+1)/2^400];
 res1=Abs[K Sum[Binomial[nn,j]r^(nn-j)h^j,{j,1,nn}]-e^2];resM=Abs[K Sum[Binomial[nn,j]r^(nn-j)hm^j,{j,1,nn}]-M e^2];midErr=Abs[g-(lo+hi)/2];
 {nn,rec["u"],inside,sign==rec["G_sign"],N[midErr,8],N[res1,8],N[resM,8],$ProcessID,starOK,mOK}];
LaunchKernels[15];pids=ParallelEvaluate[$ProcessID];Print["KERNEL_COUNT ",Length[pids]];DistributeDefinitions[calc,decExact];results=ParallelMap[calc,records,Method->"FinestGrained"];
summary={"version"->$Version,"master"->$ProcessID,"kernelCount"->Length[pids],"kernelPIDs"->pids,"cases"->Length[results],"intervalPass"->Count[results[[All,3]],True],"signPass"->Count[results[[All,4]],True],"hstarBracketPass"->Count[results[[All,9]],True],"hMBracketPass"->Count[results[[All,10]],True],"maxMidpointError"->Max[results[[All,5]]],"maxRootResidual"->Max[Join[results[[All,6]],results[[All,7]]]],"numericNotIntervalProof"->True};
Export[output,ToString[InputForm[{summary,results}]],"Text"];Print["SUMMARY ",InputForm[summary]];CloseKernels[];
