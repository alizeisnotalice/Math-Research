ClearAll[y,j,worker];
Print["Wolfram ",$Version];
Print["Conditional scalar certification under cost X>=c min(A/v^(q-1),A^((2q-1)/q)/v^(2q-1)); no geometry is asserted."];
worker[j_]:=Module[{crit,sharp,viol,lo,hi},
 crit=(j/(2j-1))^(1/(j-1));sharp=(j-1)/(2j-1)(j/(2j-1))^(j/(j-1));
 viol=Reduce[y>=0&&y^j-y^(2j-1)>1,y,Reals];
 lo=Reduce[0<=y<=1&&y^j>1,y,Reals];hi=Reduce[y>=1&&y^j>y^(2j-1),y,Reals];
 <|"q"->j,"kernel_id"->$KernelID,"pid"->$ProcessID,"beta"->ToString[j(2j-1)/(j-1),InputForm],"violation"->ToString[viol,InputForm],"branches_pass"->TrueQ[lo===False&&hi===False],"sharp_constant"->ToString[sharp,InputForm],"critical_y"->ToString[crit,InputForm]|>];
LaunchKernels[15];Print["launched: ",Length[Kernels[]]];DistributeDefinitions[worker,y];
rows=ParallelTable[worker[j],{j,3,17},Method->"CoarsestGrained"];
Print["actual IDs: ",Sort[DeleteDuplicates[Lookup[rows,"kernel_id"]]]];Print["actual PIDs: ",Sort[DeleteDuplicates[Lookup[rows,"pid"]]]];
Print["all violations False: ",And@@(#=="False"& /@ Lookup[rows,"violation"])];
Print["all branch checks: ",And@@Lookup[rows,"branches_pass"]];
Export[FileNameJoin[{DirectoryName[$InputFileName],"young_branches.json"}],rows,"RawJSON"];
CloseKernels[];Print["complete"];
