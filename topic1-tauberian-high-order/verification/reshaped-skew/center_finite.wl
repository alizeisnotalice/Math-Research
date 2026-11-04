SetDirectory[DirectoryName[$InputFileName]];Print["VERSION ",$Version," MASTER ",$ProcessID];
Clear[decExact];decExact[s_String]:=Module[{p=StringSplit[s,RegularExpression["[eE]"]],b,x=0,sgn=1,dp},b=First[p];If[Length[p]>1,x=ToExpression[p[[2]]]];If[StringStartsQ[b,"-"],sgn=-1;b=StringDrop[b,1]];dp=If[StringContainsQ[b,"."],StringLength[Last[StringSplit[b,"."]]],0];sgn FromDigits[StringReplace[b,"."->""]]10^(x-dp)];
Clear[calc];calc[rec_]:=Module[{nn=rec["n"],kap=ToExpression[rec["kappa"]],ll=ToExpression[rec["lambda"]],u=ToExpression[rec["u"]],w,rr,e,a,det,M,m0,K,h0,hs,h,A,b,D,g,idx,hr,dblo,dbhi,glo,ghi,rootOK,budgetOK,gOK},
 w=u/(1+2u)+ll u^6;rr=1-2w;e=kap u^3;a=(2u)^nn;det=1-e^2;M=a rr^nn;m0=(2u/(1+2u))^nn;D=(2w)^nn-M;K=a/(2w)^(nn-2);h0=kap^2 u^4/(4nn);
 h=hs/.FindRoot[K Sum[Binomial[nn,j]rr^(nn-j)hs^j,{j,1,nn}]==e^2,{hs,N[h0,150]},WorkingPrecision->150,AccuracyGoal->120,PrecisionGoal->120];
 A=2h NIntegrate[e-(K/e)Sum[Binomial[nn,j]rr^(nn-j)(h zz)^j,{j,1,nn}],{zz,0,1},WorkingPrecision->120,AccuracyGoal->85,PrecisionGoal->85];
 b=2A(2w)^(nn-2)+If[nn>=3,4(nn-2)h(e+h)^2(2w)^(nn-3),0];g=(M-det m0)/((det-M)(1-m0)a(2u)^6);
 idx=FromDigits[StringDrop[rec["root"]["lower_hex"],2],16];hr=TrueQ[(w-e)idx/2^320<=h<=(w-e)(idx+1)/2^320];
 dblo=decExact[rec["D_minus_B"][[1]]];dbhi=decExact[rec["D_minus_B"][[2]]];glo=decExact[rec["G"][[1]]];ghi=decExact[rec["G"][[2]]];budgetOK=TrueQ[dblo<=D-b<=dbhi];gOK=TrueQ[glo<=g<=ghi];
 {nn,rec["kappa"],rec["u"],hr,budgetOK,gOK,$ProcessID}];
Clear[centerCheck];centerCheck[rec_]:=Module[{nn=rec["n"],kap=ToExpression[rec["kappa"]],ll=ToExpression[rec["lambda"]],u=ToExpression[rec["u"]],w,rho,e,a,det,mout,M,m0,Nn,radius,upper,margin,fd,basecheck},
 w=u/(1+2u)+ll u^6;rho=1-2w;e=kap u^3;a=(2u)^nn;det=1-e^2;mout=(2w)^nn;M=a rho^nn;m0=(2u/(1+2u))^nn;Nn=Ceiling[u^(-7)];radius=rho+w;upper=(M+2e^2(2w)^(nn-2))/(2radius)^nn;margin=u^nn-upper;fd=((2Nn-5)/(2Nn))^nn M/(det-M)-m0/(1-m0);
 basecheck=calc[rec];{nn,And@@basecheck[[4;;6]],Nn==ToExpression[rec["N"]],2M>mout,margin>0,margin==ToExpression[rec["center"]["strict_margin"]],fd>0,fd==ToExpression[rec["finite_halo_minus_Ln"]["exact_fraction"]],$ProcessID}];
raw=Import["center_finite.json","RawJSON"];records=Map[Last[# ["attempts"]]&,raw["cases"]];Print["INPUT ",Length[records]];
LaunchKernels[15];pids=ParallelEvaluate[$ProcessID];DistributeDefinitions[calc,decExact,centerCheck];results=ParallelMap[centerCheck,records,Method->"FinestGrained"];
summary={"kernelCount"->Length[pids],"pids"->pids,"caseCount"->Length[results],"checksEach"->Length[results[[1,2;;8]]],"allPass"->And@@Flatten[results[[All,2;;8]]],"results"->results};
Export["center_finite_wolfram.txt",ToString[InputForm[summary]],"Text"];Print["SUMMARY ",InputForm[summary]];CloseKernels[];
