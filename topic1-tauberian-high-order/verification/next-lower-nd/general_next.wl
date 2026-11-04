SetDirectory[DirectoryName[$InputFileName]];Print["VERSION ",$Version," MASTER ",$ProcessID];
ass=Element[n,Integers]&&n>=2&&u>0&&lam>0&&kap>0;
g=1/(1+2u);shape=((g+lam u^5)^(2n-4))/(g-2lam u^6)^(n-1);
shapeRelative=(1+lam u^5/g)^(2n-4)(1-2lam u^6/g)^(-(n-1));
Print["pair shape exact relative identity=",InputForm[FullSimplify[shape==g^(n-3)shapeRelative,ass&&g-2lam u^6>0]]];
Print["pair relative series=",InputForm[FullSimplify[Normal[Series[shapeRelative,{u,0,5}]],ass]]];
Dnorm=(g+lam u^5)^n-(g-2lam u^6)^n;Dlead=n lam u^5 g^(n-2);
Print["D-minus-leading through u9=",InputForm[FullSimplify[Normal[Series[Dnorm-Dlead,{u,0,9}]],ass]]];
(* independent ratio expansion in formal X, avoiding algebraic branch changes *)
xs=xx; hs=rho((1+xs)^(1/n)-1);ares=2(hs-(((rho+hs)^(n+1)-rho^(n+1))/(n+1)-rho^n hs)/(xs rho^n));
Print["area normalized=",InputForm[FullSimplify[Series[ares/(xs rho/n),{xs,0,1}],Element[n,Integers]&&n>=2&&rho>0&&xs>0]]];
C2=n(n-1)/2;Kstar=C2 n^2/27;lam0=C2 kap^3/(8n^2 g);leading=C2/64(kap^2 g^n-kap^3/(4n)g^(n-2));
Print["kappa stationary=",InputForm[FullSimplify[Solve[D[leading,kap]==0&&kap>0,kap,Reals],Element[n,Integers]&&n>=2&&u>0]]];
Print["optimized identity=",InputForm[FullSimplify[(leading/.kap->8n g^2/3)==Kstar g^(n+4),Element[n,Integers]&&n>=2&&u>0]]];
Clear[task];task[dim_]:=Module[{gg=1/(1+2u),kap0,lambda0,lambda1,rho0,rho1,aa,m0,M0,M1,d,baseG,G1,Kn,increment,period},
 kap0=8dim gg^2/3;lambda0=Binomial[dim,2]kap0^3/(8dim^2 gg);lambda1=lambda0+c u^2;rho0=gg-2lambda0 u^6;rho1=gg-2lambda1 u^6;aa=(2u)^dim;m0=aa gg^dim;M0=aa rho0^dim;M1=aa rho1^dim;d=(1-kap0 u^3)^(dim-1)(1+(dim-1)kap0 u^3);Kn=Binomial[dim,2]dim^2/27;
 baseG=(M0-d m0)/((d-M0)(1-m0)aa(2u)^6);G1=(M1-d m0)/((d-M1)(1-m0)aa(2u)^6);
 increment=FullSimplify[Coefficient[Normal[Series[((gg+lambda1 u^5)^dim-(gg-2lambda1 u^6)^dim)-((gg+lambda0 u^5)^dim-(gg-2lambda0 u^6)^dim),{u,0,7}]],u,7]];
 period=FullSimplify[Coefficient[Normal[Series[G1-baseG,{u,0,2}]],u,2]];
 {"n"->dim,"GMinusPredThroughU1"->FullSimplify[Normal[Series[baseG-Kn gg^(dim+4),{u,0,1}]]],"LambdaSlackBudgetCoeffU7"->increment,"expectedBudget"->dim c,"periodSlackCoeffU2"->period,"expectedPeriod"->-dim c/32,"budgetCheck"->FullSimplify[increment==dim c],"periodCheck"->FullSimplify[period==-dim c/32]}];
LaunchKernels[4];pids=ParallelEvaluate[$ProcessID];ids=ParallelEvaluate[$KernelID];DistributeDefinitions[task];res=ParallelMap[Function[dim,{dim,$ProcessID,$KernelID,task[dim]}],Range[2,5],Method->"FinestGrained"];
Print["KERNELS ",Length[pids]," PIDS ",InputForm[pids]," IDS ",InputForm[ids]];Export["general_next_results.txt",ToString[InputForm[res]],"Text"];Print["RESULTS ",InputForm[res]];CloseKernels[];
