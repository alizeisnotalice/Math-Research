folder=DirectoryName[$InputFileName];
verify[n_,a_,position_]:=Module[{star,eps,mass,inc,vol,x,min,lam,z,t,rho,relaxed,formula,om,c},
 star=1/(2+n a(1-a)/(1+a));eps=Switch[position,0,star/2,1,star,2,(star+1/2)/2];
 mass=Join[{eps,1-eps},ConstantArray[a eps,n]];inc=Join[{Join[{1,1},ConstantArray[0,n]]},Table[Join[{1,0},UnitVector[n,i]],{i,n}]];vol=Join[{1},ConstantArray[(1+a)eps,n]];x=Array[xx,n+2];
 min=First[Minimize[{mass.x,And@@Join[Thread[(inc.DiagonalMatrix[mass]).x>=vol/2],Thread[x>=0],Thread[x<=1]]},x]];
 lam=Array[ll,n+1];z=Array[zz,n+2];t=vol.lam;om=1+n a eps;c=1+n a/(1+a);
 relaxed=First[Maximize[{t/2-mass.z,And@@Join[Thread[lam>=0],Thread[z>=0],Thread[z>=Transpose[inc].lam-1],{lam[[1]]<=t/om}]},Join[lam,z]]];
 formula=If[eps<=star,1/2-n a(1-a)eps/(2(1+a)),om/(2c)];
 <|"n"->n,"a"->ToString[a,InputForm],"position"->position,"epsilon"->ToString[eps,InputForm],"true_source"->ToString[min,InputForm],"responsible_relaxation"->ToString[relaxed,InputForm],"formula"->ToString[formula,InputForm],"passes"->(min==1/2&&relaxed==formula)|>];
params=Tuples[{{2,8,32,128},{1/3,1/2,2/3,1},{0,1,2}}];LaunchKernels[12];DistributeDefinitions[verify];res=ParallelMap[Function[p,TimeConstrained[MemoryConstrained[verify@@p,402653184,$Failed],180,$Aborted]],params,Method->"CoarsestGrained"];CloseKernels[];
opt=Maximize[{a(1-a)/(1+a),0<a<1},a];star=1/(2+n a(1-a)/(1+a));om=1+n a star;c=1+n a/(1+a);
checks=<|"optimal_arm_value"->ToString[opt,InputForm],"optimal_value_identity"->FullSimplify[First[opt]-(3-2Sqrt[2])],"core_half_probability"->FullSimplify[star c/om-1/2],"branch_match"->FullSimplify[om/(2c)-star],"rational_ratio"->ToString[FullSimplify[(1/2)/star/.a->1/2],InputForm]|>;
Export[FileNameJoin[{folder,"wolfram.json"}],<|"version"->$Version,"results"->res,"checks"->checks|>,"RawJSON"];Print[<|"cases"->Length[res],"passed"->Count[res,r_Association/;TrueQ[r["passes"]]]|>];Quit[];
