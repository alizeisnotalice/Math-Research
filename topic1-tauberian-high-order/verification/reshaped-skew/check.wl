SetDirectory[DirectoryName[$InputFileName]];Print["VERSION ",$Version," MASTER ",$ProcessID];
ass=Element[n,Integers]&&n>=2&&k>0&&lam>0;w=u/(1+2u)+lam u^6;rho=1-2w;rho0=1/(1+2u);det=1-k^2 u^6;
dcoef=FullSimplify[Coefficient[Normal[Series[(w/u)^n-rho^n,{u,0,5}]],u,5],ass];gcoef=FullSimplify[Coefficient[Normal[Series[rho^n-det rho0^n,{u,0,6}]],u,6]/64,ass];
H=k^2/(4n);A0=FullSimplify[2(k H-2n H^2/k),ass];B0=FullSimplify[2 A0 2^(n-2)/2^n,ass];
Print["SYMBOLIC ",InputForm[{"D/a u^5 coefficient"->dcoef,"A/u^7 coefficient"->A0,"B/a u^5 coefficient"->B0,"G limit"->gcoef,"expected G"->(k^2-2n lam)/64,"G identity"->FullSimplify[gcoef==(k^2-2n lam)/64,ass],"edge extra order"->2,"criticalLambda"->k^3/(8n^2),"criticalLimit"->(k^2-k^3/(4n))/64,"opt"->Solve[D[(k^2-k^3/(4n))/64,k]==0&&k>0,k,Reals],"optimalValue"->FullSimplify[((k^2-k^3/(4n))/64)/.k->8n/3,ass]}]];
Clear[decExact];decExact[s_String]:=Module[{p=StringSplit[s,RegularExpression["[eE]"]],b,x=0,sgn=1,dp},b=First[p];If[Length[p]>1,x=ToExpression[p[[2]]]];If[StringStartsQ[b,"-"],sgn=-1;b=StringDrop[b,1]];dp=If[StringContainsQ[b,"."],StringLength[Last[StringSplit[b,"."]]],0];sgn FromDigits[StringReplace[b,"."->""]]10^(x-dp)];
Clear[calc];calc[rec_]:=Module[{nn=rec["n"],kap=ToExpression[rec["kappa"]],ll=ToExpression[rec["lambda"]],u=ToExpression[rec["u"]],w,rr,e,a,det,M,m0,K,h0,hs,h,A,b,D,g,idx,hr,dblo,dbhi,glo,ghi,rootOK,budgetOK,gOK},
 w=u/(1+2u)+ll u^6;rr=1-2w;e=kap u^3;a=(2u)^nn;det=1-e^2;M=a rr^nn;m0=(2u/(1+2u))^nn;D=(2w)^nn-M;K=a/(2w)^(nn-2);h0=kap^2 u^4/(4nn);
 h=hs/.FindRoot[K Sum[Binomial[nn,j]rr^(nn-j)hs^j,{j,1,nn}]==e^2,{hs,N[h0,150]},WorkingPrecision->150,AccuracyGoal->120,PrecisionGoal->120];
 A=2h NIntegrate[e-(K/e)Sum[Binomial[nn,j]rr^(nn-j)(h zz)^j,{j,1,nn}],{zz,0,1},WorkingPrecision->120,AccuracyGoal->85,PrecisionGoal->85];
 b=2A(2w)^(nn-2)+If[nn>=3,4(nn-2)h(e+h)^2(2w)^(nn-3),0];g=(M-det m0)/((det-M)(1-m0)a(2u)^6);
 idx=FromDigits[StringDrop[rec["root"]["lower_hex"],2],16];hr=TrueQ[(w-e)idx/2^320<=h<=(w-e)(idx+1)/2^320];
 dblo=decExact[rec["D_minus_B"][[1]]];dbhi=decExact[rec["D_minus_B"][[2]]];glo=decExact[rec["G"][[1]]];ghi=decExact[rec["G"][[2]]];budgetOK=TrueQ[dblo<=D-b<=dbhi];gOK=TrueQ[glo<=g<=ghi];
 {nn,rec["kappa"],rec["u"],hr,budgetOK,gOK,$ProcessID}];
raw=Import["cases.json","RawJSON"];records=Select[raw["cases"],#["status"]=="certified"&];Print["INPUT ",Length[records]," OF ",Length[raw["cases"]]];
LaunchKernels[15];pids=ParallelEvaluate[$ProcessID];DistributeDefinitions[calc,decExact];results=ParallelMap[calc,records,Method->"FinestGrained"];
summary={"kernelCount"->Length[pids],"pids"->pids,"caseCount"->Length[results],"rootPass"->Count[results[[All,4]],True],"budgetPass"->Count[results[[All,5]],True],"GPass"->Count[results[[All,6]],True],"numericNotProof"->True};
Export["wolfram-results.txt",ToString[InputForm[{summary,results}]],"Text"];Print["SUMMARY ",InputForm[summary]];CloseKernels[];
