SetDirectory[DirectoryName[$InputFileName]];Print["VERSION ",$Version," MASTER ",$ProcessID];
Clear[task];task[dim_]:=Module[{K=Binomial[dim,2]dim^2/27,raw,rawFormula,b6,c6,actual,physical,poly,sigma,vroot,comp,paramTarget,q1,old,oldCheck},
 raw=Normal[Series[x^dim/((1+x)^dim-x^dim),{x,0,dim+6}]];
 rawFormula=Sum[(-1)^(6-(k-1)dim) Binomial[dim+5,6-(k-1)dim],{k,1,1+Floor[6/dim]}];
 b6=2^(dim+6)rawFormula;c6=b6+2^(dim+6)K;
 actual=Coefficient[raw,x,dim+6];
 If[dim<=6,
  poly=1+(1-q)x^dim/2+(1-q)(2-q)x^(2dim)/6+(1-q)(2-q)(3-q)x^(3dim)/24;
  sigma=x^dim poly;vroot=2x poly^(1/dim);
  physical=2^dim sigma/((1+vroot)^dim-2^dim sigma)+K vroot^(dim+6);
  comp=Normal[Series[physical,{x,0,dim+6}]];paramTarget=Factor[Coefficient[comp,x,dim+6]];q1=FullSimplify[(paramTarget/.q->1)==c6];
  old=If[dim==3,Factor[Coefficient[comp,x,dim+5]],Missing["NotApplicable"]];oldCheck=If[dim==3,FullSimplify[old==160(1-q)],Missing["NotApplicable"]],
  paramTarget=Missing["NotRequested"];q1=Missing["NotRequested"];old=Missing["NotRequested"];oldCheck=Missing["NotRequested"]];
 {"n"->dim,"K"->K,"rawJ6"->rawFormula,"rawCheck"->(actual===rawFormula),"ordinaryBaseB6"->b6,"ordinaryConditionalC6"->c6,"cZero"->(c6===0),"physicalParameterCoefficient"->paramTarget,"qOneCheck"->q1,"n3OldExponent8over3"->old,"n3OldCheck"->oldCheck}];
LaunchKernels[15];pids=ParallelEvaluate[$ProcessID];ids=ParallelEvaluate[$KernelID];Print["KERNELS ",Length[pids]," PIDS ",InputForm[pids]," IDS ",InputForm[ids]];DistributeDefinitions[task];
res=ParallelMap[Function[dim,{dim,$ProcessID,$KernelID,task[dim]}],Range[2,16],Method->"FinestGrained"];
Export["conditional_transfer_results.txt",ToString[InputForm[res]],"Text"];Print["RESULTS ",InputForm[res]];CloseKernels[];
