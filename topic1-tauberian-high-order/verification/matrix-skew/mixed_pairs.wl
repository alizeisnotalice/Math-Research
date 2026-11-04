SetDirectory[DirectoryName[$InputFileName]];
Print["VERSION ",$Version," MASTER ",$ProcessID];
(* X=p q/(K rho^n)>0, and h=rho[(1+X)^(1/n)-1]. *)
ass=Element[n,Integers]&&n>=2&&p>0&&K>0&&rho>0&&X>0;
qq=K rho^n X/p;hh=rho((1+X)^(1/n)-1);
prim=K(((rho+hh)^(n+1)-rho^(n+1))/(n+1)-rho^n hh);
area=(p+qq)(hh-prim/(p qq));lead=p qq(p+qq)/(2n K rho^(n-1));
Print["normalized series ",InputForm[FullSimplify[Series[area/lead,{X,0,1}],ass]]];
Print["primitive derivative ",InputForm[FullSimplify[D[K(((rho+h)^(n+1)-rho^(n+1))/(n+1)-rho^n h),h]==K((rho+h)^n-rho^n),ass&&h>0]]];
Print["global exceed set ",InputForm[Reduce[x>=0&&y>=0&&x y (1-(x+y)/8)>64/27,{x,y},Reals]]];
Print["equality set ",InputForm[Reduce[x>=0&&y>=0&&x y (1-(x+y)/8)==64/27,{x,y},Reals]]];
matrices={
{"n2_unequal",{{0,1},{4,0}}},
{"n2_fractional",{{0,4/3},{2/3,0}}},
{"n2_both_negative",{{0,-1},{-2,0}}},
{"n2_opposite",{{0,1},{-3,0}}},
{"n3_cycle",{{0,1,0},{0,0,2},{3,0,0}}},
{"n3_mixed",{{0,1,-2},{3,0,1},{1,-4,0}}},
{"n4_mixed_symmetric",{{0,1,-2,1/2},{1,0,3,-1},{-2,3,0,2},{1/2,-1,2,0}}},
{"n5_equal_opt",Table[If[i==j,0,40/3],{i,5},{j,5}]},
{"n5_asymmetric",{{0,1,0,-1,2},{2,0,3,0,0},{0,-1,0,4,1},{2,0,1,0,-2},{1,0,3,1,0}}}
};
Clear[matrixTask];matrixTask[entry_]:=Module[{id=entry[[1]],B=entry[[2]],dim,S2,pol},dim=Length[B];S2=Sum[B[[i,j]]B[[j,i]],{i,1,dim},{j,i+1,dim}];pol=Expand[Det[IdentityMatrix[dim]+t B]];{id,"matrix"->B,"S2"->S2,"detPolynomial"->pol,"quadraticIdentity"->(Coefficient[pol,t,2]==-S2),"linearCoefficient"->Coefficient[pol,t,1]}];
Clear[integralTask];integralTask[dim_]:=Module[{Fh=K((rho+h)^dim-rho^dim),formula,direct},formula=(p+q)h-(1/p+1/q)K Sum[Binomial[dim,j]rho^(dim-j)h^(j+1)/(j+1),{j,1,dim}];direct=Integrate[q-Fh/p+p-Fh/q,{h,0,H},Assumptions->p>0&&q>0&&K>0&&rho>0&&H>0];{"integral_n"<>ToString[dim],"check"->FullSimplify[direct==(formula/.h->H),p>0&&q>0&&K>0&&rho>0&&H>0]}];
tasks=Join[({"matrix",#}& /@ matrices),({"integral",#}& /@ Range[2,7])];
LaunchKernels[15];pids=ParallelEvaluate[$ProcessID];Print["KERNELS ",Length[pids]," PIDS ",InputForm[pids]," TASKS ",Length[tasks]];DistributeDefinitions[matrixTask,integralTask];
results=ParallelMap[Function[task,Module[{ans,id},ans=If[task[[1]]=="matrix",matrixTask[task[[2]]],integralTask[task[[2]]]];id=First[ans];{id,$ProcessID,ans}]],tasks,Method->"FinestGrained"];
Print["RESULTS ",InputForm[results]];CloseKernels[];
