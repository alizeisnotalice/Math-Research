ClearAll[u,kap,n,j,case];
Print["Wolfram ",$Version];
Print["Full outer-box terms only; fixed n=3..17, balanced even support j. f=4 n u²/(1+2u), e=kap u³, w=u/(1+2u)."];
case[n_]:=Module[{w=u/(1+2u),e=kap u^3,a=(2u)^n,f=4n u^2/(1+2u),B=n(n-1)/2,branches,inside,outside,expin,expout},
 branches=Table[
  inside=n Binomial[n-1,j-1] 2 Binomial[j-1,j/2-1] (2w)^(n-2j+2) e^(2j-1)/(2^(j-1) a f);
  outside=If[j<n,n Binomial[n-1,j] Binomial[j,j/2] 2 (2w)^(n-2j+1) e^(2j)/(2^j a f),0];
  expin=4j-3; expout=4j-1;
  <|"j"->j,"inside_coefficient"->ToString[FullSimplify[Limit[inside/u^expin,u->0,Direction->"FromAbove"]],InputForm],"outside_coefficient"->ToString[FullSimplify[Limit[outside/u^expout,u->0,Direction->"FromAbove"]],InputForm],"outside_exists"->TrueQ[j<n],"inside_expected_order"->expin,"outside_expected_order"->If[j<n,expout,Null],"inside_j2_identity"->If[j==2,TrueQ[FullSimplify[inside==2B(2w)^(n-2)e^3/(a f)]],Null],"outside_j2_identity"->If[j==2,TrueQ[FullSimplify[outside==B(n-2)(2w)^(n-3)e^4/(a f)]],Null]|>,{j,2,n,2}];
 <|"n"->n,"kernel_id"->$KernelID,"pid"->$ProcessID,"branches"->branches|>];
LaunchKernels[15];Print["launched: ",Length[Kernels[]]];DistributeDefinitions[case,u,kap];
rows=ParallelTable[case[n],{n,3,17},Method->"CoarsestGrained"];
Print["actual kernel IDs: ",Sort@DeleteDuplicates@Lookup[rows,"kernel_id"]];
Print["actual PIDs: ",Sort@DeleteDuplicates@Lookup[rows,"pid"]];
Print["j2 inside identities: ",And@@(#[["branches",1,"inside_j2_identity"]]& /@ rows)];
Print["j2 outside identities: ",And@@(#[["branches",1,"outside_j2_identity"]]& /@ rows)];
Print["j>=4 minimum expected exponents: ",{13,15}];
Export[FileNameJoin[{DirectoryName[$InputFileName],"support_branches.json"}],rows,"RawJSON"];
CloseKernels[];Print["complete"];
